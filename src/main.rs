#![cfg_attr(target_os = "windows", windows_subsystem = "windows")]

use macroquad::prelude::*;
use std::borrow::Cow;
use std::cell::{Cell, RefCell};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
#[cfg(target_os = "windows")]
use std::fs::OpenOptions;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::mpsc::{self, Receiver};

mod acid_base;
mod chemistry;
mod cml;
mod naming;
mod native;
mod pubchem;
mod query;
mod rdkit;
mod reaction;
mod shortcuts;
mod smarts;
mod smiles;

const TOP: f32 = 118.;
const LEFT: f32 = 88.;
const BOTTOM: f32 = 86.;
const RIGHT: f32 = 250.;
const BOND: f32 = 43.;
fn rotate_template_point(point: Vec2, center: Vec2, degrees: u16) -> Vec2 {
    let angle = f32::from(degrees % 360).to_radians();
    let offset = point - center;
    center
        + vec2(
            offset.x * angle.cos() - offset.y * angle.sin(),
            offset.x * angle.sin() + offset.y * angle.cos(),
        )
}
const AMINO_ACID_TEMPLATES: [&str; 20] = [
    "Glycine",
    "Alanine",
    "Valine",
    "Leucine",
    "Isoleucine",
    "Serine",
    "Threonine",
    "Cysteine",
    "Methionine",
    "Aspartic acid",
    "Asparagine",
    "Glutamic acid",
    "Glutamine",
    "Lysine",
    "Arginine",
    "Histidine",
    "Phenylalanine",
    "Tyrosine",
    "Tryptophan",
    "Proline",
];

fn builtin_template_matches_query(name: &str, query: &str) -> bool {
    let query = query.trim().to_lowercase();
    name.to_lowercase().contains(&query)
        || (matches!(
            query.as_str(),
            "amino" | "amino acid" | "amino acids" | "proteinogenic"
        ) && AMINO_ACID_TEMPLATES.contains(&name))
        || (matches!(
            query.as_str(),
            "ring" | "rings" | "fused" | "fused ring" | "fused rings"
        ) && matches!(
            name,
            "Naphthalene" | "Decalin" | "Indole" | "Quinoline" | "Anthracene" | "Phenanthrene"
        ))
        || (matches!(
            query.as_str(),
            "ring" | "rings" | "bridged" | "bridged ring" | "bridged rings"
        ) && matches!(name, "Norbornane"))
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Language {
    English,
    Danish,
}
impl Language {
    fn from_settings(value: &str) -> Self {
        if value.lines().any(|line| line.trim() == "language=da") {
            Self::Danish
        } else {
            Self::English
        }
    }
}
fn settings_path() -> PathBuf {
    let home = std::env::var_os("HOME")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from("."));
    if cfg!(target_os = "macos") {
        home.join("Library/Application Support/Molecula/settings.conf")
    } else {
        std::env::var_os("XDG_CONFIG_HOME")
            .map(PathBuf::from)
            .unwrap_or_else(|| home.join(".config"))
            .join("molecula/settings.conf")
    }
}
fn autosave_path() -> PathBuf {
    settings_path()
        .parent()
        .unwrap_or(Path::new("."))
        .join("workspace.autosave")
}
fn recent_files_path() -> PathBuf {
    settings_path()
        .parent()
        .unwrap_or(Path::new("."))
        .join("recent-files.txt")
}
fn user_templates_path() -> PathBuf {
    settings_path()
        .parent()
        .unwrap_or(Path::new("."))
        .join("templates")
}
#[derive(Clone)]
struct UserTemplate {
    name: String,
    mol: Molecule,
    path: PathBuf,
}
fn template_slug(name: &str) -> Result<String, String> {
    let name = name.trim();
    if name.is_empty()
        || name.chars().count() > 48
        || !name
            .chars()
            .all(|ch| ch.is_ascii_alphanumeric() || ch == ' ' || ch == '-' || ch == '_')
    {
        return Err("Use 1–48 letters, numbers, spaces, hyphens or underscores".into());
    }
    let slug = name
        .chars()
        .filter(|ch| *ch != ' ')
        .map(|ch| match ch {
            '_' => '-',
            value => value.to_ascii_lowercase(),
        })
        .collect::<String>();
    if slug.is_empty() {
        Err("Template name must contain a letter or number".into())
    } else {
        Ok(slug)
    }
}
fn save_user_template_to(
    directory: &Path,
    name: &str,
    mol: &Molecule,
) -> Result<UserTemplate, String> {
    if mol.atoms.is_empty() || mol.atoms.len() > 200 || mol.bonds.len() > 400 {
        return Err("A user template must contain 1–200 atoms and no more than 400 bonds".into());
    }
    if mol
        .bonds
        .iter()
        .any(|bond| bond.a >= mol.atoms.len() || bond.b >= mol.atoms.len() || bond.a == bond.b)
    {
        return Err("A user template contains an invalid bond reference".into());
    }
    let name = name.trim();
    let slug = template_slug(name)?;
    let path = directory.join(format!("{slug}.mols"));
    if path.exists() {
        return Err(format!("A template named ‘{name}’ already exists"));
    }
    fs::create_dir_all(directory).map_err(|error| error.to_string())?;
    let mut stored = mol.clone();
    stored.name = name.into();
    let encoded = native::write(&Snapshot {
        mol: stored.clone(),
        graphics: Vec::new(),
    });
    let temporary = directory.join(format!(".{slug}-{}.tmp", std::process::id()));
    let result = fs::write(&temporary, encoded)
        .and_then(|()| fs::rename(&temporary, &path))
        .map_err(|error| {
            let _ = fs::remove_file(&temporary);
            error.to_string()
        });
    result?;
    Ok(UserTemplate {
        name: name.into(),
        mol: stored,
        path,
    })
}
fn load_user_templates_from(directory: &Path) -> Result<Vec<UserTemplate>, String> {
    let entries = match fs::read_dir(directory) {
        Ok(entries) => entries,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(Vec::new()),
        Err(error) => return Err(error.to_string()),
    };
    let mut templates = Vec::new();
    for entry in entries {
        let entry = entry.map_err(|error| error.to_string())?;
        let path = entry.path();
        if path.extension().and_then(|extension| extension.to_str()) != Some("mols") {
            continue;
        }
        let text = fs::read_to_string(&path).map_err(|error| error.to_string())?;
        let snapshot = native::read(&text)
            .map_err(|error| format!("Could not load user template {}: {error}", path.display()))?;
        if snapshot.mol.name.trim().is_empty()
            || snapshot.mol.atoms.is_empty()
            || snapshot.mol.atoms.len() > 200
            || snapshot.mol.bonds.len() > 400
            || molecule_components(&snapshot.mol).len() != 1
        {
            return Err(format!("Invalid user template file: {}", path.display()));
        }
        templates.push(UserTemplate {
            name: snapshot.mol.name.clone(),
            mol: snapshot.mol,
            path,
        });
        if templates.len() > 200 {
            return Err("The user template library exceeds 200 entries".into());
        }
    }
    templates.sort_by_key(|template| template.name.to_lowercase());
    Ok(templates)
}
fn load_recent_files() -> Vec<String> {
    fs::read_to_string(recent_files_path())
        .map(|s| {
            s.lines()
                .filter(|l| !l.trim().is_empty())
                .take(10)
                .map(str::to_owned)
                .collect()
        })
        .unwrap_or_default()
}
fn save_recent_files(paths: &[String]) -> std::io::Result<()> {
    let path = recent_files_path();
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    fs::write(path, paths.join("\n"))
}
fn load_language() -> Language {
    load_language_from(&settings_path())
}
fn load_language_from(path: &Path) -> Language {
    fs::read_to_string(path)
        .map(|value| Language::from_settings(&value))
        .unwrap_or(Language::English)
}
fn save_language(language: Language) -> std::io::Result<()> {
    save_language_to(&settings_path(), language)
}
fn save_language_to(path: &Path, language: Language) -> std::io::Result<()> {
    save_setting(
        path,
        "language",
        if language == Language::Danish {
            "da"
        } else {
            "en"
        },
    )
}
fn load_dark_theme() -> bool {
    fs::read_to_string(settings_path())
        .map(|s| s.lines().any(|line| line.trim() == "theme=dark"))
        .unwrap_or(false)
}
fn load_implicit_hydrogen_visibility() -> bool {
    load_implicit_hydrogen_visibility_from(&settings_path())
}
fn load_implicit_hydrogen_visibility_from(path: &Path) -> bool {
    fs::read_to_string(path)
        .map(|settings| {
            settings
                .lines()
                .any(|line| line.trim() == "show_implicit_hydrogens=true")
        })
        .unwrap_or(false)
}
fn save_implicit_hydrogen_visibility(show: bool) -> std::io::Result<()> {
    save_implicit_hydrogen_visibility_to(&settings_path(), show)
}
fn save_implicit_hydrogen_visibility_to(path: &Path, show: bool) -> std::io::Result<()> {
    save_setting(
        path,
        "show_implicit_hydrogens",
        if show { "true" } else { "false" },
    )
}
fn load_3d_atom_symbols_visibility() -> bool {
    load_3d_atom_symbols_visibility_from(&settings_path())
}
fn load_3d_atom_symbols_visibility_from(path: &Path) -> bool {
    fs::read_to_string(path)
        .map(|settings| {
            !settings
                .lines()
                .any(|line| line.trim() == "show_3d_atom_symbols=false")
        })
        .unwrap_or(true)
}
fn save_3d_atom_symbols_visibility(show: bool) -> std::io::Result<()> {
    save_3d_atom_symbols_visibility_to(&settings_path(), show)
}
fn save_3d_atom_symbols_visibility_to(path: &Path, show: bool) -> std::io::Result<()> {
    save_setting(
        path,
        "show_3d_atom_symbols",
        if show { "true" } else { "false" },
    )
}
fn load_zoom_preference() -> f32 {
    load_zoom_preference_from(&settings_path())
}
fn load_zoom_preference_from(path: &Path) -> f32 {
    fs::read_to_string(path)
        .ok()
        .and_then(|settings| {
            settings
                .lines()
                .find_map(|line| line.trim().strip_prefix("zoom=")?.parse::<f32>().ok())
        })
        .filter(|zoom| zoom.is_finite())
        .unwrap_or(1.)
        .clamp(0.2, 4.0)
}
fn save_zoom_preference_to(path: &Path, zoom: f32) -> std::io::Result<()> {
    if !zoom.is_finite() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            "Zoom preference must be finite",
        ));
    }
    save_setting(path, "zoom", &format!("{:.4}", zoom.clamp(0.2, 4.0)))
}
fn save_dark_theme(dark: bool) -> std::io::Result<()> {
    save_setting(
        &settings_path(),
        "theme",
        if dark { "dark" } else { "light" },
    )
}
fn save_setting(path: &Path, key: &str, value: &str) -> std::io::Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    let old = fs::read_to_string(path).unwrap_or_default();
    let mut lines: Vec<String> = old
        .lines()
        .filter(|line| !line.trim().starts_with(&format!("{key}=")))
        .map(str::to_string)
        .collect();
    lines.push(format!("{key}={value}"));
    fs::write(path, format!("{}\n", lines.join("\n")))
}
fn localized<'a>(language: Language, source: &'a str) -> Cow<'a, str> {
    if language == Language::English {
        return Cow::Borrowed(source);
    }
    let exact = match source {
        "File" => "Filer",
        "Edit" => "Rediger",
        "View" => "Vis",
        "Insert" => "Indsæt",
        "Atom" => "Atom",
        "Bond" => "Binding",
        "Structure" => "Struktur",
        "Calculations" => "Beregninger",
        "Services" => "Tjenester",
        "Help" => "Hjælp",
        "New" => "Nyt",
        "Open" => "Åbn",
        "Save" => "Gem",
        "Undo" => "Fortryd",
        "Redo" => "Gentag",
        "Delete selection" => "Slet markering",
        "Fit" => "Tilpas visning",
        "Elements" => "Grundstoffer",
        "Templates" => "Skabeloner",
        "Arrow" => "Reaktionspil",
        "Equilibrium arrow" => "Ligevægtspil",
        "Resonance arrow" => "Resonanspil",
        "Retrosynthesis arrow" => "Retrosyntesepil",
        "Electron pair arrow" => "Elektronparspil",
        "Single electron arrow" => "Enkelt-elektronpil",
        "Plus" => "Plus",
        "Text" => "Tekst",
        "Line" => "Linje",
        "Rectangle" => "Rektangel",
        "Ellipse" => "Ellipse",
        "Charge +" => "Ladning +",
        "Charge -" => "Ladning -",
        "Single" => "Enkeltbinding",
        "Double" => "Dobbeltbinding",
        "Triple" => "Tripelbinding",
        "Aromatic" => "Aromatisk",
        "Wedge" => "Fyldt kile",
        "Hash" => "Stiplet kile",
        "Wavy" => "Bølgebinding",
        "Validate" => "Kontrollér struktur",
        "Find substructure matches" => "Find delstrukturmatch",
        "Name structure" => "Navngiv struktur",
        "Auto-map reaction (Complete)" => "Automatisk atomkortlægning (Komplet)",
        "Auto-map reaction (Matching)" => "Automatisk atomkortlægning (Matchende)",
        "Auto-map reaction (Changing)" => "Automatisk atomkortlægning (Ændringer)",
        "Clean 2D" => "Ryd op i 2D",
        "Add explicit hydrogens" => "Tilføj eksplicitte hydrogenatomer",
        "Remove explicit hydrogens" => "Fjern eksplicitte hydrogenatomer",
        "Rotate 90 CW" => "Rotér 90° med uret",
        "Rotate 90 CCW" => "Rotér 90° mod uret",
        "Rotate 180" => "Rotér 180°",
        "Flip left-right" => "Vend venstre-højre",
        "Flip up-down" => "Vend op-ned",
        "Benzene template" => "Benzen-skabelon",
        "Formula & mass" => "Formel og masse",
        "Molecular properties" => "Molekyleegenskaber",
        "Chemistry A acid/base tools" => "Kemi A: syre-baseværktøjer",
        "Chemistry A acid-base workbench" => "Kemi A: syre-baseberegner",
        "Buffer pH from pKa" => "Buffer-pH ud fra pKa",
        "pKa from buffer measurements" => "pKa ud fra buffermåling",
        "Weak-acid pH" => "pH for svag syre",
        "Weak-base pH" => "pH for svag base",
        "Strong acid/base mixture pH" => "pH for blanding af stærk syre/base",
        "Polyprotic speciation" => "Speciering af polyprot syre",
        "Bjerrum diagram data" => "Bjerrumdiagram-data",
        "Polyprotic titration pH" => "pH ved polyprot titrering",
        "Polyprotic titration curve" => "Polyprot titrerkurve",
        "Bjerrum diagram · species distribution" => "Bjerrumdiagram · artsfordeling",
        "Fraction (%)" => "Andel (%)",
        "Added strong base (mL)" => "Tilsat stærk base (mL)",
        "Polyprotic acid · titration curve" => "Polyprot syre · titrerkurve",
        "Copy CSV data" => "Kopiér CSV-data",
        "Choose a calculation, then enter values in the order shown. Use ; between fields and | between successive pKa values." => "Vælg en beregning, og indtast værdierne i den viste rækkefølge. Brug ; mellem felter og | mellem successive pKa-værdier.",
        "Henderson–Hasselbalch: pH = pKa + log₁₀([A⁻]/[HA]). Charts can copy data as CSV." => "Henderson–Hasselbalch: pH = pKa + log₁₀([A⁻]/[HA]). Diagramdata kan kopieres som CSV.",
        "Enter values separated by semicolons. Decimal comma or point is accepted." => "Angiv værdier adskilt af semikolon. Decimal-komma eller punktum accepteres.",
        "Example: 4.76; 0.10; 0.20" => "Eksempel: 4,76; 0,10; 0,20",
        "Example: 5.06; 0.10; 0.20" => "Eksempel: 5,06; 0,10; 0,20",
        "Example: 4.76; 0.10" => "Eksempel: 4,76; 0,10",
        "Example: 9.25; 0.10" => "Eksempel: 9,25; 0,10",
        "Example: 0.10; 0.05" => "Eksempel: 0,10; 0,05",
        "Example: 2.15|7.20|12.35; 7.00" => "Eksempel: 2,15|7,20|12,35; 7,00",
        "Example: 2.15|7.20; 0; 14; 0.5" => "Eksempel: 2,15|7,20; 0; 14; 0,5",
        "Example: 4.76; 25; 0.10; 0.10; 12.5" => "Eksempel: 4,76; 25; 0,10; 0,10; 12,5",
        "Example: 4.76; 25; 0.10; 0.10; 1.5" => "Eksempel: 4,76; 25; 0,10; 0,10; 1,5",
        "pKa values (in increasing order) separated by |; then pH" => "pKa-værdier i stigende rækkefølge adskilt af |; derefter pH",
        "pKa values separated by |; pH minimum; maximum; step" => "pKa-værdier adskilt af |; pH minimum; maksimum; trin",
        "pKa values separated by |; acid volume mL; acid M; base M; added base volume mL" => "pKa-værdier adskilt af |; syrevolumen mL; syre M; base M; tilsat basevolumen mL",
        "pKa values separated by |; acid volume mL; acid M; base M; max equivalents" => "pKa-værdier adskilt af |; syrevolumen mL; syre M; base M; maks. ækvivalenter",
        "Buffer pH (pKa; acid M; conjugate base M)" => "Buffer-pH (pKa; syre M; korresponderende base M)",
        "Measured pH; acid M; conjugate base M" => "Målt pH; syre M; korresponderende base M",
        "pKa of conjugate acid; total base concentration M" => "pKa for korresponderende syre; total basekoncentration M",
        "Final strong-acid concentration M; final strong-base concentration M" => "Slutkoncentration af stærk syre M; slutkoncentration af stærk base M",
        "pKa; weak-acid concentration M" => "pKa; koncentration af svag syre M",
        "pKa is entered from your class data or a cited table. No structure-based pKa prediction is made." => "pKa indtastes fra undervisningsdata eller en kilde. Der beregnes ikke pKa ud fra strukturen.",
        "Local topological descriptors" => "Lokale topologiske deskriptorer",
        "Edit molecule data" => "Redigér molekyldata",
        "Copy as SMILES" => "Kopiér som SMILES",
        "Copy as SMARTS" => "Kopiér som SMARTS",
        "Copy drawing as image" => "Kopiér tegning som billede",
        "Show implicit hydrogens" => "Vis implicitte hydrogenatomer",
        "Implicit hydrogen labels" => "Visning af implicitte hydrogenatomer",
        "Shown" => "Vises",
        "Hidden" => "Skjult",
        "Interface language" => "Sprog i brugerfladen",
        "Appearance" => "Udseende",
        "Language and appearance are saved automatically." => "Sprog og udseende gemmes automatisk.",
        "KEY=VALUE saves or updates · -KEY removes · NAME sets molecule name" => "KEY=VALUE gemmer eller opdaterer · -KEY fjerner · NAME angiver molekylenavn",
        "Reaction analysis" => "Reaktionsanalyse",
        "Reaction coefficient: custom" => "Reaktionskoefficient: brugerdefineret",
        "Gasteiger partial charges" => "Gasteiger-partialladninger",
        "Set a whole-number coefficient from 1 to 99 for this component." => "Angiv en heltalskoefficient fra 1 til 99 for denne komponent.",
        "Auto-balance reaction" => "Afstem reaktion automatisk",
        "Auto-balance all reaction steps" => "Afstem alle reaktionstrin automatisk",
        "View previous reaction step" => "Vis forrige reaktionstrin",
        "View next reaction step" => "Vis næste reaktionstrin",
        "View all reaction steps" => "Vis alle reaktionstrin",
        "RDKit descriptors" => "RDKit-deskriptorer",
        "Enumerate tautomers" => "Enumerér tautomere former",
        "Enumerate stereoisomers" => "Enumerér stereoisomerer",
        "Fetch PubChem properties" => "Hent PubChem-data",
        "Search PubChem in browser" => "Søg i PubChem i browseren",
        "Assign R/S labels with RDKit" => "Tildel R/S-labels med RDKit",
        "Integration status" => "Integrationsstatus",
        "Quick start guide" => "Hurtig vejledning",
        "Shortcuts" => "Tastaturgenveje",
        "Feature status" => "Funktionsstatus",
        "Settings" => "Indstillinger",
        "Edit text above" => "Redigér tekst over pilen",
        "Edit text below" => "Redigér tekst under pilen",
        "Change arrow type" => "Skift piltype",
        "Document" => "Dokument",
        "Local graph properties" => "Lokale grafegenskaber",
        "Elements:" => "Grundstoffer:",
        "Atoms:" => "Atomer:",
        "Bonds:" => "Bindinger:",
        "Method: local valence/H count and embedded mass table." => "Metode: lokal valens-/H-beregning og indbygget massetabel.",
        "Method: graph counts and embedded mass table." => "Metode: graftælling og indbygget massetabel.",
        "Unusual valence, query chemistry and unsupported elements may be unresolved." => "Usædvanlig valens, søgekemi og ikke-understøttede grundstoffer kan være uafklarede.",
        "logP, TPSA, pKa, NMR and 3D require a verified chemistry engine." => "logP, TPSA, pKa, NMR og 3D kræver en verificeret kemimotor.",
        "Select and drag atoms or graphic objects (Esc)" => "Vælg, eller træk atomer og grafiske objekter (Esc)",
        "Delete one atom, bond or object (E)" => "Slet et atom, en binding eller et objekt (E)",
        "Draw a single bond (1)" => "Tegn en enkeltbinding (1)",
        "Draw or change a bond to a double bond (2)" => "Tegn eller ændr en binding til en dobbeltbinding (2)",
        "Draw a triple bond (3)" => "Tegn en tripelbinding (3)",
        "Draw an aromatic bond (4)" => "Tegn en aromatisk binding (4)",
        "Draw a filled wedge (5)" => "Tegn en fyldt kile (5)",
        "Draw a hashed wedge (6)" => "Tegn en stiplet kile (6)",
        "Draw an unspecified stereo bond (7)" => "Tegn en binding med uspecificeret stereokemi (7)",
        "Changed existing bond to a double bond" => "Ændrede den eksisterende binding til en dobbeltbinding",
        "Place or replace an atom (element symbol)" => "Placér eller udskift et atom (grundstofsymbol)",
        "Draw a zigzag carbon chain (drag between endpoints)" => "Tegn en zigzag-kæde af kulstof (træk mellem endepunkter)",
        "Place a six-membered ring (Shift+R)" => "Placér en seksleddet ring (Shift+R)",
        "Place a five-membered ring" => "Placér en femleddet ring",
        "Increase atom charge (+)" => "Forøg atomets ladning (+)",
        "Decrease atom charge (-)" => "Sænk atomets ladning (-)",
        "Drag to draw a reaction arrow" => "Træk for at tegne en reaktionspil",
        "Drag to draw a single-electron fishhook arrow" => "Træk for at tegne en fiskekrogspil for ét elektron",
        "Drag to draw a two-electron curved arrow" => "Træk for at tegne en buet pil for to elektroner",
        "Place a reaction plus sign" => "Placér et plustegn mellem reaktanter",
        "Place a text label" => "Placér en tekstetiket",
        "Drag to draw a line" => "Træk for at tegne en linje",
        "Drag to draw a rectangle" => "Træk for at tegne et rektangel",
        "Drag to draw an ellipse" => "Træk for at tegne en ellipse",
        "Click canvas to place template" => "Klik på lærredet for at placere skabelonen",
        "MOLECULA · QUICK START" => "MOLECULA · HURTIG VEJLEDNING",
        "01  DRAW" => "01  TEGN",
        "Choose a bond on the left, then drag or click two points. Choose Chain near the top of the tools to draw a zigzag carbon chain." => "Vælg en binding til venstre, og træk eller klik på to punkter. Vælg Kæde for at tegne en zigzag-kæde af kulstof.",
        "To change an existing bond, choose Double (=) and click its middle. Undo reverses the change." => "Hvis du vil ændre en binding, vælg Dobbeltbinding (=), og klik midt på bindingen. Fortryd ruller ændringen tilbage.",
        "02  CHANGE ATOMS" => "02  SKIFT ATOMER",
        "Choose C, N, O or another element below. Click an atom to replace it. Use the periodic table for all elements." => "Vælg C, N, O eller et andet grundstof nedenfor. Klik på et atom for at udskifte det. Brug periodesystemet til alle grundstoffer.",
        "With an atom selected, type @ followed by a mass number (for example @13); @0 clears the isotope. Type label:Ph to show an alias; label:- clears it." => "Når et atom er valgt, skriv @ efterfulgt af et massetal (fx @13); @0 fjerner isotopen. Skriv label:Ph for at vise et alias; label:- fjerner det.",
        "03  BUILD A REACTION" => "03  BYG EN REAKTION",
        "Place reactants left and products right. Insert → Arrow draws a reaction arrow. Use Calculations → Reaction analysis when you want formula, molar mass and balance. Use e⁻ tools for curved electron flow." => "Placér reaktanter til venstre og produkter til højre. Indsæt → Reaktionspil tegner en pil. Brug Beregninger → Reaktionsanalyse til formel, molmasse og balance. Brug e⁻-værktøjerne til buet elektronbevægelse.",
        "Automatic mapping: select every complete reactant and product component in one step; leave agents unselected, assign explicit roles, then choose Structure → Auto-map reaction. Requires optional AtomMapper + RDKit." => "Automatisk kortlægning: markér alle hele reaktant- og produktkomponenter i ét trin; lad agenter være umarkerede, angiv eksplicitte roller, og vælg derefter Struktur → Automatisk atomkortlægning. Kræver valgfri AtomMapper + RDKit.",
        "04  NAME AND CHECK" => "04  NAVNGIV OG KONTROLLÉR",
        "Structure → Name structure gives a name in the supported teaching subset. Molecular properties shows formula and mass." => "Struktur → Navngiv struktur giver et navn inden for det understøttede undervisningsområde. Molekyleegenskaber viser formel og masse.",
        "05  SAVE" => "05  GEM",
        "File → Save writes a .mols document with editable chemistry and graphics." => "Filer → Gem opretter et .mols-dokument med redigerbar kemi og grafik.",
        "UNDO Cmd/Ctrl+Z   ·   REDO Cmd/Ctrl+Y   ·   DELETE Delete   ·   CANCEL Esc" => "FORTRYD Cmd/Ctrl+Z   ·   GENTAG Cmd/Ctrl+Y   ·   SLET Delete   ·   ANNULLÉR Esc",
        "Sel" => "Vælg",
        "Del" => "Slet",
        "Chain" => "Kæde",
        "->" => "Pil",
        "e2" => "e-par",
        "e1" => "e−",
        "+R" => "Plus",
        "Txt" => "Tekst",
        "Wd" => "Kile",
        "Hs" => "Stiplet",
        "Wv" => "Bølge",
        "INSPECTOR" => "INSPEKTØR",
        "Properties" => "Egenskaber",
        "Canvas selection" => "Markering på lærredet",
        "Selected bond" => "Valgt binding",
        "Reaction arrow" => "Reaktionspil",
        "Electron-flow arrow" => "Elektronpil",
        "Graphic object" => "Grafisk objekt",
        "Nothing selected" => "Intet valgt",
        "Drag a box to select atoms" => "Træk en boks for at vælge atomer",
        "or click a bond/object" => "eller klik på binding/objekt",
        "Whole document" => "Hele dokumentet",
        "Analysis & information" => "Analyse og information",
        "Close" => "Luk",
        "Language preference saved" => "Sprogvalg gemt",
        "Could not save language preference" => "Sprogvalget kunne ikke gemmes",
        "Enter text" => "Indtast tekst",
        "Cancel" => "Annullér",
        "Apply" => "Anvend",
        "Periodic table - choose an element" => "Periodesystem – vælg et grundstof",
        "Choose a tool or drag on the canvas" => "Vælg et værktøj, eller træk på lærredet",
        "Start drawing chemistry" => "Begynd at tegne kemi",
        "Ready: draw a bond or choose a tool" => "Klar: tegn en binding, eller vælg et værktøj",
        "Analysis" => "Analyse",
        "ELEMENTS" => "GRUNDSTOFFER",
        "TOOLS" => "VÆRKTØJER",
        "English" | "Dansk" => source,
        "Language" => "Sprog",
        "Saved automatically" => "Gemmes automatisk",
        "New document" => "Nyt dokument",
        "Open MOL or MOLS" => "Åbn MOL eller MOLS",
        "Save molecule" => "Gem molekyle",
        "Export may lose information" => "Eksport kan miste information",
        "Save as .mols to preserve the full document." => "Gem som .mols for at bevare hele dokumentet.",
        "Use an absolute path or a path relative to the launch folder." => "Brug en absolut sti eller en sti relativ til programmets mappe.",
        "Draw a reaction, equilibrium, or retrosynthesis arrow first" => "Tegn først en reaktions-, ligevægts- eller retrosyntesepil",
        "RDKit and Open Babel are not installed." => "RDKit og Open Babel er ikke installeret.",
        "3D generation, IUPAC names, advanced descriptors, InChI and broad format conversion are unavailable." => "3D-generering, IUPAC-navne, avancerede deskriptorer, InChI og bred filkonvertering er utilgængelige.",
        "Local drawing, MOL and MOLS still work." => "Lokal tegning samt MOL og MOLS fungerer stadig.",
        "RDKit descriptor worker failed" => "RDKit-beregningen mislykkedes",
        "Calculating RDKit descriptors…" => "Beregner RDKit-deskriptorer…",
        "Draw a molecule first" => "Tegn først et molekyle",
        "Click canvas to place benzene" => "Klik på lærredet for at placere benzen",
        "Select an atom before assigning an R-group" => "Vælg et atom, før du tilføjer en R-gruppe",
        "Select an atom before assigning a map number" => "Vælg et atom, før du tilføjer et atomkortnummer",
        "Select an atom before assigning a query" => "Vælg et atom, før du tilføjer et søgekriterium",
        "Calculating with RDKit..." => "Beregner med RDKit...",
        "Arrow placed" => "Pil placeret",
        "Undid edit" => "Redigering fortrudt",
        "Redid edit" => "Redigering gentaget",
        "No implicit hydrogens to expand" => "Ingen implicitte hydrogenatomer at udvide",
        "Added explicit hydrogen atoms" => "Tilføjede eksplicitte hydrogenatomer",
        "No attached explicit hydrogens to remove" => "Ingen tilknyttede eksplicitte hydrogenatomer at fjerne",
        "Removed explicit hydrogen atoms" => "Fjernede eksplicitte hydrogenatomer",
        "Draw a structure first" => "Tegn først en struktur",
        "Draw at least two atoms before cleaning" => "Tegn mindst to atomer før oprydning",
        "Cleaned 2D coordinates using a spring layout" => "Ryddede op i 2D-koordinater med fjederlayout",
        "MOL V2000 cannot store more than 999 atoms or bonds" => "MOL V2000 kan højst gemme 999 atomer eller bindinger",
        "Unsupported extension. Use .mols, .mol, .sdf, .rxn, .smi, .smarts, .svg, .png or .pdf" => "Filtypen understøttes ikke. Brug .mols, .mol, .sdf, .rxn, .smi, .smarts, .svg, .png eller .pdf",
        "Calculation finished" => "Beregning færdig",
        "Discarded calculation for an earlier molecule" => "Kasserede beregning for et tidligere molekyle",
        "Click the arrow end point, or drag from the start" => "Klik på pilens slutpunkt, eller træk fra startpunktet",
        "GRAPH" => "GRAF",
        "Clear iso" => "Fjern isotop",
        "Query predicates" => "Søgekriterier",
        "Aliphatic" => "Alifatisk",
        "Unsaturated" => "Umættet",
        "Any atom" => "Vilkårligt atom",
        "Heteroatom" => "Heteroatom",
        "Clear query" => "Fjern søgning",
        "Radical" => "Radikal",
        "Implicit H" => "Implicit H",
        "Save selection as template" => "Gem markering som skabelon",
        "Built-in" => "Indbyggede",
        "My templates" => "Mine skabeloner",
        "Save selected fragment" => "Gem markeret fragment",
        "Previous" => "Forrige",
        "Next" => "Næste",
        "Confirm delete" => "Bekræft sletning",
        "Name custom template" => "Navngiv skabelon",
        "Rename custom template" => "Omdøb brugerdefineret skabelon",
        "Rename" => "Omdøb",
        "auto" => "automatisk",
        "none" => "ingen",
        "singlet" => "singlet",
        "doublet" => "doublet",
        "triplet" => "triplet",
        "Delete bond" => "Slet binding",
        "Order −" => "Orden −",
        "Order +" => "Orden +",
        "Make double (=)" => "Gør til dobbeltbinding (=)",
        "E/Z applies to double bonds" => "E/Z gælder dobbeltbindinger",
        "Export" => "Eksportér",
        "Type: reaction" => "Type: reaktion",
        "Type: equilibrium" => "Type: ligevægt",
        "Type: resonance" => "Type: resonans",
        "Type: retrosynthesis" => "Type: retrosyntese",
        "No basic valence warnings found. Advanced stereochemistry and structure checks require a chemistry engine." => "Ingen grundlæggende valensadvarsler. Avanceret stereokemi og strukturkontrol kræver en kemimotor.",
        "See STATUS.md in the source folder for implemented, partial and unavailable functions." => "Se STATUS.md i programmets mappe for funktioner, der er implementeret, delvise eller utilgængelige.",
        "Cmd/Ctrl+N: New | O: Open | S: Save | Z: Undo | Y: Redo" => "Cmd/Ctrl+N: Nyt | O: Åbn | S: Gem | Z: Fortryd | Y: Gentag",
        "Delete: delete selection | Esc: select tool | Shift+R: six-membered ring" => "Delete: slet markering | Esc: vælgerværktøj | Shift+R: seksleddet ring",
        "Drawing buffer: 1/2/3 orders, 4 aromatic, 5 wedge, 6 hash, 7 wavy" => "Tegneinput: 1/2/3 bindinger, 4 aromatisk, 5 kile, 6 stiplet, 7 bølget",
        "C, N, O, Cl... select element; + and - change selected charge" => "C, N, O, Cl... vælg grundstof; + og - ændrer ladningen",
        "R1...R32767 assigns R-group; M1... assigns a map; M0 clears; M+/M= assign unique maps" => "R1...R32767 angiver R-gruppe; M1... angiver atomkort",
        "Ctrl+E opens the periodic table." => "Ctrl+E åbner periodesystemet.",
        _ => return localized_dynamic(source),
    };
    Cow::Borrowed(exact)
}
fn quick_start_guide(language: Language) -> String {
    [
        "MOLECULA · QUICK START",
        "01  DRAW",
        "Choose a bond on the left, then drag or click two points. Choose Chain near the top of the tools to draw a zigzag carbon chain.",
        "To change an existing bond, choose Double (=) and click its middle. Undo reverses the change.",
        "02  CHANGE ATOMS",
        "Choose C, N, O or another element below. Click an atom to replace it. Use the periodic table for all elements.",
        "With an atom selected, type @ followed by a mass number (for example @13); @0 clears the isotope. Type label:Ph to show an alias; label:- clears it.",
        "03  BUILD A REACTION",
        "Place reactants left and products right. Insert → Arrow draws a reaction arrow. Use Calculations → Reaction analysis when you want formula, molar mass and balance. Use e⁻ tools for curved electron flow.",
        "Automatic mapping: select every complete reactant and product component in one step; leave agents unselected, assign explicit roles, then choose Structure → Auto-map reaction. Requires optional AtomMapper + RDKit.",
        "04  NAME AND CHECK",
        "Structure → Name structure gives a name in the supported teaching subset. Molecular properties shows formula and mass.",
        "05  SAVE",
        "File → Save writes a .mols document with editable chemistry and graphics.",
        "UNDO Cmd/Ctrl+Z   ·   REDO Cmd/Ctrl+Y   ·   DELETE Delete   ·   CANCEL Esc",
    ]
    .into_iter()
    .map(|line| localized(language, line).into_owned())
    .collect::<Vec<_>>()
    .join("\n\n")
}
fn localized_dynamic(source: &str) -> Cow<'_, str> {
    for (prefix, replacement) in [
        ("Charge:", "Ladning:"),
        ("Isotope:", "Isotop:"),
        ("Map:", "Atomkort:"),
        ("Query:", "Søgekriterium:"),
        ("Bond order:", "Bindingsorden:"),
        ("Style:", "Stil:"),
        ("Stereo:", "Stereokemi:"),
        ("Internal bonds:", "Interne bindinger:"),
        ("Formal charge:", "Formel ladning:"),
        ("Formula:", "Formel:"),
        ("Molar:", "Molmasse:"),
        ("Exact:", "Eksakt:"),
        ("Selected implicit H:", "Valgt implicit H:"),
        ("Warnings:", "Advarsler:"),
        ("Elements:", "Grundstoffer:"),
        ("Atoms:", "Atomer:"),
        ("Bonds:", "Bindinger:"),
        ("Reactants:", "Reaktanter:"),
        ("Products:", "Produkter:"),
        ("Agents:", "Reagenser:"),
        ("Balance unavailable:", "Balance utilgængelig:"),
        ("Method:", "Metode:"),
        ("Element: ", "Grundstof: "),
        ("Selected ", "Valgte "),
        ("Saved ", "Gemt "),
        ("Opened ", "Åbnede "),
        ("Drew a ", "Tegnede en "),
    ] {
        if let Some(rest) = source.strip_prefix(prefix) {
            return Cow::Owned(format!("{replacement}{rest}"));
        }
    }
    if let Some(count) = source.strip_suffix(" atoms") {
        return Cow::Owned(format!("{count} atomer"));
    }
    if source.contains(" atoms  ·  ") && source.ends_with(" bonds") {
        return Cow::Owned(
            source
                .replace(" atoms", " atomer")
                .replace(" bonds", " bindinger"),
        );
    }
    Cow::Borrowed(source)
}
fn ink() -> Color {
    UI_DARK_MODE.with(|dark| {
        if dark.get() {
            Color::from_rgba(232, 239, 234, 255)
        } else {
            Color::from_rgba(25, 49, 39, 255)
        }
    })
}
fn muted() -> Color {
    UI_DARK_MODE.with(|dark| {
        if dark.get() {
            Color::from_rgba(153, 171, 160, 255)
        } else {
            Color::from_rgba(103, 129, 115, 255)
        }
    })
}
fn navy() -> Color {
    UI_DARK_MODE.with(|dark| {
        if dark.get() {
            Color::from_rgba(15, 22, 18, 255)
        } else {
            Color::from_rgba(20, 49, 39, 255)
        }
    })
}
fn accent() -> Color {
    Color::from_rgba(36, 139, 94, 255)
}
fn surface() -> Color {
    UI_DARK_MODE.with(|dark| {
        if dark.get() {
            Color::from_rgba(22, 29, 25, 255)
        } else {
            Color::from_rgba(246, 250, 247, 255)
        }
    })
}
fn card_surface() -> Color {
    UI_DARK_MODE.with(|dark| {
        if dark.get() {
            Color::from_rgba(31, 40, 34, 255)
        } else {
            WHITE
        }
    })
}
fn panel_surface() -> Color {
    UI_DARK_MODE.with(|dark| {
        if dark.get() {
            Color::from_rgba(26, 34, 29, 255)
        } else {
            Color::from_rgba(244, 247, 252, 255)
        }
    })
}
fn round_rect(x: f32, y: f32, w: f32, h: f32, r: f32, c: Color) {
    let r = r.min(w / 2.).min(h / 2.);
    draw_rectangle(x + r, y, w - 2. * r, h, c);
    draw_rectangle(x, y + r, r, h - 2. * r, c);
    draw_rectangle(x + w - r, y + r, r, h - 2. * r, c);
    for (cx, cy) in [
        (x + r, y + r),
        (x + w - r, y + r),
        (x + r, y + h - r),
        (x + w - r, y + h - r),
    ] {
        draw_circle(cx, cy, r, c);
    }
}
fn round_border(x: f32, y: f32, w: f32, h: f32, r: f32, c: Color) {
    round_rect(x, y, w, h, r, c);
    round_rect(x + 1., y + 1., w - 2., h - 2., (r - 1.).max(0.), WHITE);
}
fn draw_shape(kind: ShapeKind, a: Vec2, b: Vec2, color: Color) {
    let x = a.x.min(b.x);
    let y = a.y.min(b.y);
    let w = (a.x - b.x).abs().max(1.);
    let h = (a.y - b.y).abs().max(1.);
    match kind {
        ShapeKind::Line => draw_line(a.x, a.y, b.x, b.y, 2., color),
        ShapeKind::Rectangle => draw_rectangle_lines(x, y, w, h, 2., color),
        ShapeKind::Ellipse => {
            for k in 0..48 {
                let t = std::f32::consts::TAU * k as f32 / 48.;
                let u = std::f32::consts::TAU * (k + 1) as f32 / 48.;
                draw_line(
                    x + w / 2. + w / 2. * t.cos(),
                    y + h / 2. + h / 2. * t.sin(),
                    x + w / 2. + w / 2. * u.cos(),
                    y + h / 2. + h / 2. * u.sin(),
                    2.,
                    color,
                );
            }
        }
    }
}
fn tool_button(x: f32, y: f32, w: f32, h: f32, t: Tool, active: bool) -> bool {
    let p = vec2(mouse_position().0, mouse_position().1);
    let hover = p.x >= x && p.x < x + w && p.y >= y && p.y < y + h;
    if active {
        round_rect(x, y, w, h, 11., accent());
    } else if hover {
        round_rect(x, y, w, h, 11., Color::from_rgba(37, 78, 61, 255));
    }
    let c = if active || hover {
        WHITE
    } else {
        Color::from_rgba(175, 205, 187, 255)
    };
    let cx = x + w / 2.;
    let cy = y + h / 2.;
    match t {
        Tool::Select => {
            draw_line(cx - 9., cy - 10., cx - 9., cy + 8., 2., c);
            draw_line(cx - 9., cy - 10., cx + 7., cy + 2., 2., c);
            draw_line(cx + 7., cy + 2., cx - 9., cy + 2., 2., c);
        }
        Tool::Erase => {
            draw_line(cx - 10., cy + 5., cx + 4., cy - 9., 4., c);
            draw_line(cx + 4., cy - 9., cx + 11., cy - 2., 4., c);
            draw_line(cx - 3., cy + 11., cx + 11., cy + 11., 1.5, c);
        }
        Tool::Bond(n) => {
            for k in 0..n {
                let o = (k as f32 - (n as f32 - 1.) / 2.) * 4.;
                draw_line(cx - 12., cy + 7. + o, cx + 12., cy - 7. + o, 2., c);
            }
        }
        Tool::StyledBond(BondStyle::Aromatic) => {
            draw_line(cx - 12., cy + 6., cx + 12., cy - 6., 2., c);
            draw_line(cx - 7., cy + 10., cx + 7., cy - 10., 1., c);
        }
        Tool::StyledBond(BondStyle::Wedge) => {
            draw_triangle(
                vec2(cx - 12., cy + 7.),
                vec2(cx + 11., cy - 10.),
                vec2(cx + 9., cy - 2.),
                c,
            );
        }
        Tool::StyledBond(BondStyle::Hash) => {
            for k in 0..5 {
                let xx = cx - 10. + k as f32 * 5.;
                draw_line(
                    xx,
                    cy + 7. - k as f32 * 3.,
                    xx + 3.,
                    cy + 9. - k as f32 * 2.,
                    1.5,
                    c,
                );
            }
        }
        Tool::StyledBond(BondStyle::Wavy) => {
            for k in 0..5 {
                let xx = cx - 12. + k as f32 * 5.;
                draw_line(
                    xx,
                    cy + if k % 2 == 0 { 3. } else { -3. },
                    xx + 5.,
                    cy + if k % 2 == 0 { -3. } else { 3. },
                    1.8,
                    c,
                );
            }
        }
        Tool::Atom => {
            draw_text("C", cx - 7., cy + 7., 21., c);
        }
        Tool::Chain => {
            draw_line(cx - 12., cy + 6., cx - 4., cy - 6., 2., c);
            draw_line(cx - 4., cy - 6., cx + 4., cy + 6., 2., c);
            draw_line(cx + 4., cy + 6., cx + 12., cy - 6., 2., c);
        }
        Tool::Ring(n) => {
            let count = n;
            for k in 0..count {
                let a = std::f32::consts::TAU * k as f32 / n as f32;
                let b = std::f32::consts::TAU * (k + 1) as f32 / n as f32;
                draw_line(
                    cx + 11. * a.cos(),
                    cy + 11. * a.sin(),
                    cx + 11. * b.cos(),
                    cy + 11. * b.sin(),
                    1.7,
                    c,
                );
            }
        }
        Tool::Arrow => {
            draw_line(cx - 12., cy, cx + 11., cy, 2., c);
            draw_line(cx + 11., cy, cx + 3., cy - 6., 2., c);
            draw_line(cx + 11., cy, cx + 3., cy + 6., 2., c);
        }
        Tool::Shape(ShapeKind::Line) => {
            draw_line(cx - 11., cy + 8., cx + 11., cy - 8., 2., c);
        }
        Tool::Shape(ShapeKind::Rectangle) => {
            draw_rectangle_lines(cx - 11., cy - 8., 22., 16., 2., c);
        }
        Tool::Shape(ShapeKind::Ellipse) => {
            for k in 0..24 {
                let a = std::f32::consts::TAU * k as f32 / 24.;
                let b = std::f32::consts::TAU * (k + 1) as f32 / 24.;
                draw_line(
                    cx + 11. * a.cos(),
                    cy + 8. * a.sin(),
                    cx + 11. * b.cos(),
                    cy + 8. * b.sin(),
                    1.7,
                    c,
                );
            }
        }
        Tool::Electron(n) => {
            let a = vec2(cx - 12., cy + 8.);
            let b = vec2(cx + 11., cy + 5.);
            let control = vec2(cx, cy - 13.);
            let mut prev = a;
            for k in 1..=12 {
                let t = k as f32 / 12.;
                let p = a * (1. - t).powi(2) + control * (2. * t * (1. - t)) + b * t * t;
                draw_line(prev.x, prev.y, p.x, p.y, 1.8, c);
                prev = p;
            }
            draw_line(b.x, b.y, b.x - 6., b.y - 6., 1.8, c);
            if n == 2 {
                draw_line(b.x, b.y, b.x - 7., b.y + 3., 1.8, c);
            }
        }
        Tool::Plus => {
            draw_line(cx - 9., cy, cx + 9., cy, 2., c);
            draw_line(cx, cy - 9., cx, cy + 9., 2., c);
        }
        Tool::Text => draw_text("T", cx - 6., cy + 7., 21., c),
        Tool::Charge(n) => {
            draw_text(if n > 0 { "+" } else { "−" }, cx - 7., cy + 7., 22., c);
        }
        _ => draw_circle(cx, cy, 4., c),
    }
    hover && is_mouse_button_pressed(MouseButton::Left) && UI_INPUT_ENABLED.with(Cell::get)
}
fn swap_stereo(s: &mut String) {
    let replacement = match s.as_str() {
        "R" => "S",
        "S" => "R",
        "r" => "s",
        "s" => "r",
        "CW" => "CCW",
        "CCW" => "CW",
        "E" => "Z",
        "Z" => "E",
        "@" => "@@",
        "@@" => "@",
        "cis" => "trans",
        "trans" => "cis",
        _ => return,
    };
    *s = replacement.into();
}
thread_local! {
    static UI_FONT: RefCell<Option<Font>> = const { RefCell::new(None) };
    static UI_LANGUAGE: Cell<Language> = const { Cell::new(Language::English) };
    static UI_INPUT_ENABLED: Cell<bool> = const { Cell::new(true) };
    static UI_DARK_MODE: Cell<bool> = const { Cell::new(false) };
}
fn draw_text(s: &str, x: f32, y: f32, size: f32, color: Color) {
    let text = UI_LANGUAGE.with(|language| localized(language.get(), s));
    UI_FONT.with(|f| {
        let borrowed = f.borrow();
        macroquad::text::draw_text_ex(
            text.as_ref(),
            x,
            y,
            TextParams {
                font: borrowed.as_ref(),
                font_size: size as u16,
                color,
                ..Default::default()
            },
        );
    });
}
fn measure_text(s: &str, font: Option<&Font>, size: u16, scale: f32) -> TextDimensions {
    UI_FONT.with(|f| {
        let borrowed = f.borrow();
        macroquad::text::measure_text(s, font.or(borrowed.as_ref()), size, scale)
    })
}

fn window_conf() -> Conf {
    // Install before Macroquad creates its WGL context so driver panics are
    // recorded too (the previous location ran only after context creation).
    #[cfg(target_os = "windows")]
    {
        install_startup_diagnostics();
        // The Windows installer ships Mesa's OpenGL runtime beside the app.
        // Force its CPU renderer so WGL works even when the system driver lacks
        // WGL_ARB_pixel_format (the crash reported on the affected PC).
        std::env::set_var("GALLIUM_DRIVER", "llvmpipe");
        std::env::set_var("LIBGL_ALWAYS_SOFTWARE", "1");
    }

    Conf {
        window_title: "molecula — Structure Studio".into(),
        window_width: 1400,
        window_height: 900,
        window_resizable: true,
        high_dpi: true,
        // Use the Molecula mark for the running window and macOS Dock as well
        // as the app bundle's Finder icon. Miniquad's default is its red/black
        // logo, which otherwise appears in the Dock while this app is open.
        icon: Some(macroquad::miniquad::conf::Icon {
            small: *include_bytes!("../assets/molecula-16.rgba"),
            medium: *include_bytes!("../assets/molecula-32.rgba"),
            big: *include_bytes!("../assets/molecula-64.rgba"),
        }),
        // Avoid requiring a 4x multisample WGL pixel format. Some older and
        // remote Windows graphics drivers fail during context creation when
        // that format is unavailable.
        sample_count: 1,
        ..Default::default()
    }
}

#[cfg(target_os = "windows")]
fn install_startup_diagnostics() {
    let directory = std::env::var_os("LOCALAPPDATA")
        .map(PathBuf::from)
        .unwrap_or_else(std::env::temp_dir)
        .join("Molecula");
    let _ = fs::create_dir_all(&directory);
    let log_path = directory.join("startup.log");
    let _ = fs::write(
        &log_path,
        "Molecula started; waiting for editor initialization.\n",
    );

    let previous_hook = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        if let Ok(mut file) = OpenOptions::new().create(true).append(true).open(&log_path) {
            let _ = writeln!(file, "Startup error: {info}");
            let _ = writeln!(
                file,
                "Backtrace:\n{}",
                std::backtrace::Backtrace::force_capture()
            );
        }
        previous_hook(info);
    }));
}

#[cfg(target_os = "windows")]
fn record_startup_note(note: &str) {
    let directory = std::env::var_os("LOCALAPPDATA")
        .map(PathBuf::from)
        .unwrap_or_else(std::env::temp_dir)
        .join("Molecula");
    if let Ok(mut file) = OpenOptions::new()
        .create(true)
        .append(true)
        .open(directory.join("startup.log"))
    {
        let _ = writeln!(file, "{note}");
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Atom {
    pos: Vec2,
    z: Option<f32>,
    coord3d: Option<Vec3>,
    element: String,
    charge: i8,
    isotope: Option<u16>,
    radical: u8,
    explicit_h: u8,
    implicit_h_override: Option<u8>,
    valence_override: Option<u8>,
    alias: Option<String>,
    atom_map: Option<u16>,
    aromatic: bool,
    stereo: Option<String>,
    r_group: Option<u16>,
    query: Option<String>,
    reaction_role: Option<String>,
    reaction_step: u16,
    coefficient: u8,
}
impl Atom {
    fn new(pos: Vec2, element: &str) -> Self {
        Self {
            pos,
            z: None,
            coord3d: None,
            element: element.into(),
            charge: 0,
            isotope: None,
            radical: 0,
            explicit_h: 0,
            implicit_h_override: None,
            valence_override: None,
            alias: None,
            atom_map: None,
            aromatic: false,
            stereo: None,
            r_group: None,
            query: None,
            reaction_role: None,
            reaction_step: 1,
            coefficient: 1,
        }
    }
}
fn cycle_valence_override(current: Option<u8>) -> Option<u8> {
    match current {
        None => Some(1),
        Some(0) => None,
        Some(1..=13) => current.map(|value| value + 1),
        Some(_) => Some(0),
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum BondStyle {
    Normal,
    Aromatic,
    Wedge,
    Hash,
    Wavy,
    Coordinate,
    Any,
}
impl BondStyle {
    fn code(self) -> &'static str {
        match self {
            Self::Normal => "normal",
            Self::Aromatic => "aromatic",
            Self::Wedge => "wedge",
            Self::Hash => "hash",
            Self::Wavy => "wavy",
            Self::Coordinate => "coordinate",
            Self::Any => "any",
        }
    }
    fn parse(s: &str) -> Option<Self> {
        Some(match s {
            "normal" => Self::Normal,
            "aromatic" => Self::Aromatic,
            "wedge" => Self::Wedge,
            "hash" => Self::Hash,
            "wavy" => Self::Wavy,
            "coordinate" => Self::Coordinate,
            "any" => Self::Any,
            _ => return None,
        })
    }
}
fn set_bond_order_value(bond: &mut Bond, order: u8) {
    let new_order = order.clamp(1, 4);
    if bond.order != new_order || new_order != 2 {
        bond.stereo = None;
    }
    bond.order = new_order;
    if bond.order == 4 {
        bond.style = BondStyle::Aromatic;
    } else {
        bond.style = BondStyle::Normal;
    }
    if let Some(query) = bond.query.as_deref() {
        let retained: Vec<&str> = query
            .split(';')
            .filter(|term| {
                !matches!(
                    *term,
                    "any" | "single-or-double" | "single-or-aromatic" | "double-or-aromatic"
                )
            })
            .collect();
        bond.query = (!retained.is_empty()).then(|| retained.join(";"));
    }
}
fn next_bond_style(style: BondStyle) -> BondStyle {
    match style {
        BondStyle::Normal | BondStyle::Aromatic => BondStyle::Wedge,
        BondStyle::Wedge => BondStyle::Hash,
        BondStyle::Hash => BondStyle::Wavy,
        BondStyle::Wavy => BondStyle::Coordinate,
        BondStyle::Coordinate => BondStyle::Any,
        BondStyle::Any => BondStyle::Normal,
    }
}
fn next_double_bond_stereo(stereo: Option<&str>) -> Option<String> {
    match stereo {
        None | Some("unspecified") => Some("E".into()),
        Some("E") => Some("Z".into()),
        Some("Z") => None,
        Some("cis") => Some("trans".into()),
        Some("trans") => Some("cis".into()),
        Some(_) => Some("E".into()),
    }
}
#[derive(Clone, Debug, PartialEq)]
struct Bond {
    a: usize,
    b: usize,
    order: u8,
    style: BondStyle,
    stereo: Option<String>,
    query: Option<String>,
}
#[derive(Clone, Debug, Default, PartialEq)]
struct Molecule {
    atoms: Vec<Atom>,
    bonds: Vec<Bond>,
    name: String,
    properties: Vec<(String, String)>,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ShapeKind {
    Line,
    Rectangle,
    Ellipse,
}
impl ShapeKind {
    fn code(&self) -> &'static str {
        match self {
            Self::Line => "line",
            Self::Rectangle => "rectangle",
            Self::Ellipse => "ellipse",
        }
    }
    fn parse(s: &str) -> Option<Self> {
        Some(match s {
            "line" => Self::Line,
            "rectangle" => Self::Rectangle,
            "ellipse" => Self::Ellipse,
            _ => return None,
        })
    }
}
#[derive(Clone, Debug, PartialEq)]
enum Graphic {
    Arrow {
        start: Vec2,
        end: Vec2,
        kind: String,
        above: String,
        below: String,
    },
    Text {
        pos: Vec2,
        text: String,
        size: f32,
        bold: bool,
        italic: bool,
    },
    Plus {
        pos: Vec2,
    },
    ElectronArrow {
        start: Vec2,
        control: Vec2,
        end: Vec2,
        electrons: u8,
    },
    Shape {
        kind: ShapeKind,
        start: Vec2,
        end: Vec2,
    },
}
#[derive(Clone)]
struct Snapshot {
    mol: Molecule,
    graphics: Vec<Graphic>,
}
#[derive(Clone)]
struct DocumentTab {
    title: String,
    path: Option<String>,
    snapshot: Snapshot,
    dirty: bool,
    undo: Vec<Snapshot>,
    redo: Vec<Snapshot>,
}
impl DocumentTab {
    fn blank() -> Self {
        Self {
            title: "Untitled".into(),
            path: None,
            snapshot: Snapshot {
                mol: Molecule::default(),
                graphics: Vec::new(),
            },
            dirty: false,
            undo: Vec::new(),
            redo: Vec::new(),
        }
    }
}
#[derive(Clone, Copy, PartialEq)]
enum Tool {
    Select,
    Erase,
    Bond(u8),
    StyledBond(BondStyle),
    Atom,
    Chain,
    Ring(usize),
    Charge(i8),
    Arrow,
    Electron(u8),
    Plus,
    Text,
    Template(&'static str),
    Shape(ShapeKind),
}
#[derive(Clone, Copy)]
enum TextTarget {
    New(Vec2),
    ArrowAbove(usize),
    ArrowBelow(usize),
    TemplateName,
    RenameUserTemplate(usize),
    MoleculeProperty,
    ReactionCoefficient,
    BufferPh,
    PkaFromBuffer,
    WeakAcidPh,
    WeakBasePh,
    StrongMixturePh,
    PolyproticSpeciation,
    BjerrumDiagram,
    PolyproticTitrationPh,
    PolyproticTitrationCurve,
}
fn acid_base_prompt(target: TextTarget) -> Option<&'static str> {
    match target {
        TextTarget::BufferPh => Some("Buffer pH (pKa; acid M; conjugate base M)"),
        TextTarget::PkaFromBuffer => Some("Measured pH; acid M; conjugate base M"),
        TextTarget::WeakAcidPh => Some("pKa; weak-acid concentration M"),
        TextTarget::WeakBasePh => Some("pKa of conjugate acid; total base concentration M"),
        TextTarget::StrongMixturePh => {
            Some("Final strong-acid concentration M; final strong-base concentration M")
        }
        TextTarget::PolyproticSpeciation => {
            Some("pKa values (in increasing order) separated by |; then pH")
        }
        TextTarget::BjerrumDiagram => Some("pKa values separated by |; pH minimum; maximum; step"),
        TextTarget::PolyproticTitrationPh => {
            Some("pKa values separated by |; acid volume mL; acid M; base M; added base volume mL")
        }
        TextTarget::PolyproticTitrationCurve => {
            Some("pKa values separated by |; acid volume mL; acid M; base M; max equivalents")
        }
        _ => None,
    }
}
fn acid_base_example(target: TextTarget) -> Option<&'static str> {
    match target {
        TextTarget::BufferPh => Some("Example: 4.76; 0.10; 0.20"),
        TextTarget::PkaFromBuffer => Some("Example: 5.06; 0.10; 0.20"),
        TextTarget::WeakAcidPh => Some("Example: 4.76; 0.10"),
        TextTarget::WeakBasePh => Some("Example: 9.25; 0.10"),
        TextTarget::StrongMixturePh => Some("Example: 0.10; 0.05"),
        TextTarget::PolyproticSpeciation => Some("Example: 2.15|7.20|12.35; 7.00"),
        TextTarget::BjerrumDiagram => Some("Example: 2.15|7.20; 0; 14; 0.5"),
        TextTarget::PolyproticTitrationPh => Some("Example: 4.76; 25; 0.10; 0.10; 12.5"),
        TextTarget::PolyproticTitrationCurve => Some("Example: 4.76; 25; 0.10; 0.10; 1.5"),
        _ => None,
    }
}
fn localize_acid_base_report(language: Language, report: String) -> String {
    if language != Language::Danish {
        return report;
    }
    [
        ("Buffer pH =", "Buffer-pH ="),
        ("Measured pH =", "Målt pH ="),
        ("Weak-acid pH =", "pH for svag syre ="),
        ("Weak-base pH =", "pH for svag base ="),
        ("Strong acid/base mixture pH =", "pH for blanding af stærk syre/base ="),
        ("Polyprotic acid speciation at pH", "Speciering af polyprot syre ved pH"),
        ("Bjerrum distribution table", "Bjerrumdiagram: fordelingsdata"),
        ("Titration pH =", "Titrerings-pH ="),
        ("Titration curve points", "Titrerkurve: datapunkter"),
        ("pKa values:", "pKa-værdier:"),
        ("pKa =", "pKa ="),
        ("Initial acid concentration", "Begyndelseskoncentration af syre"),
        ("Initial concentration =", "Begyndelseskoncentration ="),
        ("Initial base concentration =", "Begyndelseskoncentration af base ="),
        ("pKa of BH⁺ =", "pKa for BH⁺ ="),
        ("Initial acid:", "Startsyre:"),
        ("Added strong base:", "Tilsat stærk base:"),
        ("Final strong-acid concentration", "Slutkoncentration af stærk syre"),
        ("Final strong-base concentration", "Slutkoncentration af stærk base"),
        ("[HA] =", "[HA] ="),
        ("[A⁻] =", "[A⁻] ="),
        ("Method: Henderson–Hasselbalch, pH = pKa + log₁₀([A⁻]/[HA]). Aqueous ideal-solution calculation.", "Metode: Henderson–Hasselbalch, pH = pKa + log₁₀([A⁻]/[HA]). Ideel vandig opløsning."),
        ("Method: rearranged Henderson–Hasselbalch. This derives a pKa from supplied measurement data; it is not a structure-based prediction.", "Metode: Henderson–Hasselbalch omskrevet til pKa. pKa udledes af de indtastede måledata og forudsiges ikke ud fra strukturen."),
        ("Method: ideal aqueous monoprotic equilibrium solved by charge balance, including water autoionization; Kw = 1.0×10⁻¹⁴ at 25 °C.", "Metode: ideel vandig monoprotisk ligevægt løst ved ladningsbalance, inklusive vands autoprotolyse; Kw = 1,0×10⁻¹⁴ ved 25 °C."),
        ("Method: ideal aqueous conjugate-acid equilibrium solved by charge balance, including water autoionization; Kw = 1.0×10⁻¹⁴ at 25 °C.", "Metode: ideel vandig ligevægt for den korresponderende syre løst ved ladningsbalance, inklusive vands autoprotolyse; Kw = 1,0×10⁻¹⁴ ved 25 °C."),
        ("Method: neutralization and water ion balance, Kw = 1.0×10⁻¹⁴ at 25 °C. Enter concentrations after any volume change.", "Metode: neutralisation og ionbalance i vand, Kw = 1,0×10⁻¹⁴ ved 25 °C. Indtast koncentrationerne efter eventuel volumenændring."),
        ("Method: successive ideal aqueous dissociation equilibria; fractions sum to 100%.", "Metode: successive ideelle vandige dissociationsligevægte; andelene summerer til 100 %."),
        ("Method: polyprotic distribution plus electroneutrality after dilution; ideal aqueous conditions, 25 °C.", "Metode: polyprot fordeling og elektroneutralitet efter fortynding; ideelle vandige forhold ved 25 °C."),
        ("Method: equilibrium charge balance for a polyprotic acid titrated with a strong base; ideal aqueous conditions, 25 °C.", "Metode: ladningsbalance i ligevægt for en polyprot syre titreret med en stærk base; ideelle vandige forhold ved 25 °C."),
        ("Method:", "Metode:"),
        ("Aqueous ideal-solution calculation.", "Beregning for en ideel vandig opløsning."),
        ("at 25 °C.", "ved 25 °C."),
        ("pH\t", "pH\t"),
        ("α (%)", "α (%)"),
        ("V(base), mL", "V(base), mL"),
        ("Each α is the fraction of one successive protonation state. Ideal aqueous model; supplied pKa values are not predicted.", "Hver α er andelen af én protoneringstilstand. Ideel vandig model; indtastede pKa-værdier forudsiges ikke."),
        ("Table values can be plotted externally.", "Værdierne vises også som diagram og kan kopieres som CSV."),
    ]
    .into_iter()
    .fold(report, |text, (from, to)| text.replace(from, to))
}
fn valid_property_key(key: &str) -> bool {
    !key.is_empty()
        && key.len() <= 64
        && key
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || byte == b'_' || byte == b'-')
}
fn clipboard_candidates() -> Vec<(&'static str, &'static [&'static str])> {
    #[cfg(target_os = "macos")]
    {
        vec![("pbcopy", &[])]
    }
    #[cfg(target_os = "windows")]
    {
        vec![("clip.exe", &[])]
    }
    #[cfg(all(unix, not(target_os = "macos")))]
    {
        vec![
            ("wl-copy", &[]),
            ("xclip", &["-selection", "clipboard"]),
            ("xsel", &["--clipboard", "--input"]),
        ]
    }
}
fn copy_text_to_clipboard(text: &str) -> Result<(), String> {
    let candidates = clipboard_candidates();
    let mut failure = String::new();
    for (program, args) in candidates {
        let child = Command::new(program)
            .args(args)
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn();
        let Ok(mut child) = child else {
            failure = format!("{program} is unavailable");
            continue;
        };
        if let Some(mut input) = child.stdin.take() {
            if let Err(error) = input.write_all(text.as_bytes()) {
                drop(input);
                let _ = child.wait();
                failure = format!("Could not write to {program}: {error}");
                continue;
            }
        }
        match child.wait() {
            Ok(status) if status.success() => return Ok(()),
            Ok(_) => failure = format!("{program} could not access the system clipboard"),
            Err(error) => failure = format!("Could not finish {program}: {error}"),
        }
    }
    Err(if failure.is_empty() {
        "No supported system clipboard command is available".into()
    } else {
        format!("Could not copy to clipboard: {failure}")
    })
}
#[cfg(any(target_os = "macos", target_os = "windows"))]
fn process_failure(output: &std::process::Output) -> String {
    let stderr = String::from_utf8_lossy(&output.stderr).trim().to_string();
    if !stderr.is_empty() {
        return stderr;
    }
    let stdout = String::from_utf8_lossy(&output.stdout).trim().to_string();
    if !stdout.is_empty() {
        return stdout;
    }
    format!("process exited with {}", output.status)
}
#[cfg(any(target_os = "macos", target_os = "windows"))]
fn write_clipboard_png(png: &[u8]) -> Result<PathBuf, String> {
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_err(|error| error.to_string())?
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "molecula-clipboard-{}-{nonce}.png",
        std::process::id()
    ));
    fs::write(&path, png).map_err(|error| format!("Could not stage PNG for clipboard: {error}"))?;
    Ok(path)
}
fn copy_png_to_clipboard(png: &[u8]) -> Result<(), String> {
    #[cfg(target_os = "macos")]
    {
        let path = write_clipboard_png(png)?;
        let script = "on run argv\n set the clipboard to (read (POSIX file (item 1 of argv)) as «class PNGf»)\nend run";
        let result = Command::new("osascript")
            .arg("-e")
            .arg(script)
            .arg(&path)
            .output()
            .map_err(|error| format!("Could not start osascript: {error}"))
            .and_then(|output| {
                if output.status.success() {
                    Ok(())
                } else {
                    Err(process_failure(&output))
                }
            });
        let _ = fs::remove_file(path);
        return result;
    }
    #[cfg(target_os = "windows")]
    {
        let path = write_clipboard_png(png)?;
        let script = "Add-Type -AssemblyName System.Windows.Forms; Add-Type -AssemblyName System.Drawing; $ErrorActionPreference = 'Stop'; $image = [System.Drawing.Image]::FromFile($env:MOLECULA_CLIPBOARD_PNG); try { [System.Windows.Forms.Clipboard]::SetImage($image) } finally { $image.Dispose() }";
        let result = Command::new("powershell.exe")
            .args(["-NoProfile", "-STA", "-Command", script])
            .env("MOLECULA_CLIPBOARD_PNG", &path)
            .output()
            .map_err(|error| format!("Could not start PowerShell: {error}"))
            .and_then(|output| {
                if output.status.success() {
                    Ok(())
                } else {
                    Err(process_failure(&output))
                }
            });
        let _ = fs::remove_file(path);
        return result;
    }
    #[cfg(all(unix, not(target_os = "macos")))]
    {
        let mut failure = String::new();
        for (program, args) in [
            ("wl-copy", &(["--type", "image/png"][..])),
            (
                "xclip",
                &(["-selection", "clipboard", "-t", "image/png", "-i"][..]),
            ),
        ] {
            let child = Command::new(program)
                .args(args)
                .stdin(Stdio::piped())
                .stdout(Stdio::null())
                .stderr(Stdio::null())
                .spawn();
            let Ok(mut child) = child else {
                failure = format!("{program} is unavailable");
                continue;
            };
            if let Some(mut input) = child.stdin.take() {
                if let Err(error) = input.write_all(png) {
                    drop(input);
                    let _ = child.wait();
                    failure = format!("Could not write PNG to {program}: {error}");
                    continue;
                }
            }
            match child.wait() {
                Ok(status) if status.success() => return Ok(()),
                Ok(_) => failure = format!("{program} could not access the image clipboard"),
                Err(error) => failure = format!("Could not finish {program}: {error}"),
            }
        }
        return Err(if failure.is_empty() {
            "No supported image clipboard command is available".into()
        } else {
            format!("Could not copy image to clipboard: {failure}")
        });
    }
    #[allow(unreachable_code)]
    Err("Image clipboard is unavailable on this platform".into())
}
enum LessonContent {
    Naming { name: String, steps: Vec<String> },
    Reaction { feedback: Option<String> },
}
struct TeachingLesson {
    content: LessonContent,
    step: usize,
}
#[derive(Clone)]
struct ChemistryChart {
    title: String,
    x_label: String,
    y_label: String,
    series: Vec<String>,
    points: Vec<(f64, Vec<f64>)>,
}
fn chemistry_chart_csv(chart: &ChemistryChart) -> String {
    let mut csv = chart.x_label.clone();
    for label in &chart.series {
        csv.push(',');
        csv.push_str(label);
    }
    csv.push('\n');
    for (x, values) in &chart.points {
        csv.push_str(&format!("{x:.6}"));
        for value in values {
            csv.push_str(&format!(",{value:.8}"));
        }
        csv.push('\n');
    }
    csv
}
struct Editor {
    mol: Molecule,
    graphics: Vec<Graphic>,
    documents: Vec<DocumentTab>,
    active_document: usize,
    dirty: bool,
    recent_files: Vec<String>,
    recent_open: bool,
    pending_close: Option<usize>,
    close_after_save: Option<usize>,
    pending_recovery: Option<(Vec<DocumentTab>, usize)>,
    pending_quit: bool,
    quit_when_clean: bool,
    quit_now: bool,
    undo: Vec<Snapshot>,
    redo: Vec<Snapshot>,
    tool: Tool,
    element: String,
    selected: Option<usize>,
    selected_atoms: Vec<usize>,
    selected_bond: Option<usize>,
    selected_graphic: Option<usize>,
    selected_bonds: Vec<usize>,
    selected_graphics: Vec<usize>,
    zoom: f32,
    preferred_zoom: f32,
    zoom_preference_dirty: bool,
    pan: Vec2,
    panning: bool,
    naming_highlight: Vec<usize>,
    reaction_highlight_bonds: Vec<usize>,
    reaction_step_view: Option<u16>,
    view_3d_open: bool,
    rotation_3d: Vec2,
    zoom_3d: f32,
    show_3d_atom_symbols: bool,
    last_3d_pointer: Option<Vec2>,
    view_3d_style: u8,
    drag_start: Option<Vec2>,
    selection_box: Option<Vec2>,
    moved: bool,
    file_dialog: Option<bool>,
    path_input: String,
    status: String,
    quick: shortcuts::QuickInput,
    pending_export: Option<String>,
    periodic_open: bool,
    templates_open: bool,
    template_query: String,
    user_templates: Vec<UserTemplate>,
    template_tab_custom: bool,
    template_page: usize,
    template_rotation: u16,
    template_preview: Option<Molecule>,
    template_attachment_target: Option<usize>,
    pending_delete_template: Option<usize>,
    pending_custom_template: Option<Molecule>,
    text_target: Option<TextTarget>,
    text_buffer: String,
    bond_query: Option<String>,
    menu_open: Option<usize>,
    palette_open: bool,
    palette_query: String,
    palette_selection: usize,
    info_overlay: Option<String>,
    acid_base_open: bool,
    acid_base_chart: Option<ChemistryChart>,
    info_scroll: usize,
    lesson: Option<TeachingLesson>,
    lesson_open: bool,
    settings_open: bool,
    language: Language,
    dark_theme: bool,
    show_implicit_hydrogens: bool,
    arrow_kind: String,
    pending_arrow: Option<Vec2>,
    calc_rx: Option<(Molecule, Receiver<Result<String, String>>)>,
    auto_map_pending: Option<(Vec<usize>, Vec<usize>, String)>,
}
impl Editor {
    fn refuse_overlapping_worker(&mut self) -> bool {
        if self.calc_rx.is_some() {
            self.status =
                "Wait for the current background chemistry job to finish before starting another."
                    .into();
            true
        } else {
            false
        }
    }

    fn set_language(&mut self, language: Language) {
        self.language = language;
        UI_LANGUAGE.with(|current| current.set(language));
        self.status = match save_language(language) {
            Ok(()) => "Language preference saved".into(),
            Err(_) => "Could not save language preference".into(),
        };
    }
    fn set_dark_theme(&mut self, dark: bool) {
        self.dark_theme = dark;
        UI_DARK_MODE.with(|mode| mode.set(dark));
        self.status = match save_dark_theme(dark) {
            Ok(()) => "Appearance saved".into(),
            Err(_) => "Could not save appearance preference".into(),
        };
    }
    fn persist_zoom_preference(&mut self) {
        self.persist_zoom_preference_to(&settings_path());
    }
    fn persist_zoom_preference_to(&mut self, path: &Path) {
        if !self.zoom_preference_dirty {
            return;
        }
        self.zoom_preference_dirty = false;
        if save_zoom_preference_to(path, self.preferred_zoom).is_err() {
            self.status = "Could not save zoom preference".into();
        }
    }
    fn set_implicit_hydrogen_visibility(&mut self, show: bool) {
        self.show_implicit_hydrogens = show;
        self.status = match save_implicit_hydrogen_visibility(show) {
            Ok(()) => "Hydrogen display preference saved".into(),
            Err(_) => "Could not save hydrogen display preference".into(),
        };
    }
    fn toggle_3d_atom_symbols(&mut self) {
        self.show_3d_atom_symbols = !self.show_3d_atom_symbols;
        let visible = if self.show_3d_atom_symbols {
            "3D atom symbols shown"
        } else {
            "3D atom symbols hidden"
        };
        self.status = match save_3d_atom_symbols_visibility(self.show_3d_atom_symbols) {
            Ok(()) => visible.into(),
            Err(error) => format!("{visible}, but the preference could not be saved: {error}"),
        };
    }
    fn finish_arrow(&mut self, p: Vec2) {
        let drag = self.drag_start.unwrap_or(p);
        let start = if drag.distance(p) > 18. {
            Some(drag)
        } else {
            self.pending_arrow.filter(|a| a.distance(p) > 18.)
        };
        if let Some(start) = start {
            self.checkpoint();
            match self.tool {
                Tool::Arrow => self.graphics.push(Graphic::Arrow {
                    start,
                    end: p,
                    kind: self.arrow_kind.clone(),
                    above: String::new(),
                    below: String::new(),
                }),
                Tool::Electron(electrons) => self.graphics.push(Graphic::ElectronArrow {
                    start,
                    control: electron_control(start, p),
                    end: p,
                    electrons,
                }),
                _ => return,
            }
            self.pending_arrow = None;
            self.status = "Arrow placed".into();
        } else {
            self.pending_arrow = Some(p);
            self.status = "Click the arrow end point, or drag from the start".into();
        }
    }
    fn new() -> Self {
        let preferred_zoom = load_zoom_preference();
        let (user_templates, template_load_error) =
            match load_user_templates_from(&user_templates_path()) {
                Ok(templates) => (templates, None),
                Err(error) => (Vec::new(), Some(error)),
            };
        let mut editor = Self {
            mol: Molecule::default(),
            graphics: vec![],
            documents: vec![DocumentTab::blank()],
            active_document: 0,
            dirty: false,
            recent_files: load_recent_files(),
            recent_open: false,
            pending_close: None,
            close_after_save: None,
            pending_recovery: None,
            pending_quit: false,
            quit_when_clean: false,
            quit_now: false,
            undo: vec![],
            redo: vec![],
            tool: Tool::Bond(1),
            element: "C".into(),
            selected: None,
            selected_atoms: Vec::new(),
            selected_bond: None,
            selected_graphic: None,
            selected_bonds: Vec::new(),
            selected_graphics: Vec::new(),
            zoom: preferred_zoom,
            preferred_zoom,
            zoom_preference_dirty: false,
            pan: vec2(0., 0.),
            panning: false,
            naming_highlight: Vec::new(),
            reaction_highlight_bonds: Vec::new(),
            reaction_step_view: None,
            view_3d_open: false,
            rotation_3d: Vec2::ZERO,
            zoom_3d: 1.,
            show_3d_atom_symbols: load_3d_atom_symbols_visibility(),
            last_3d_pointer: None,
            view_3d_style: 1,
            drag_start: None,
            selection_box: None,
            moved: false,
            file_dialog: None,
            path_input: format!(
                "{}/Documents/molecule.mols",
                std::env::var("HOME").unwrap_or_else(|_| ".".into())
            ),
            status: "Ready: draw a bond or choose a tool".into(),
            quick: shortcuts::QuickInput::default(),
            pending_export: None,
            periodic_open: false,
            templates_open: false,
            template_query: String::new(),
            user_templates,
            template_tab_custom: false,
            template_page: 0,
            template_rotation: 0,
            template_preview: None,
            template_attachment_target: None,
            pending_delete_template: None,
            pending_custom_template: None,
            text_target: None,
            text_buffer: String::new(),
            bond_query: None,
            menu_open: None,
            palette_open: false,
            palette_query: String::new(),
            palette_selection: 0,
            info_overlay: None,
            acid_base_open: false,
            acid_base_chart: None,
            info_scroll: 0,
            lesson: None,
            lesson_open: false,
            settings_open: false,
            language: load_language(),
            dark_theme: load_dark_theme(),
            show_implicit_hydrogens: load_implicit_hydrogen_visibility(),
            arrow_kind: "reaction".into(),
            pending_arrow: None,
            calc_rx: None,
            auto_map_pending: None,
        };
        if let Some(error) = template_load_error {
            editor.status = format!("Could not load custom templates: {error}");
        }
        editor.restore_autosave();
        editor
    }
    fn view_origin(w: f32, h: f32) -> Vec2 {
        vec2((LEFT + w - RIGHT) / 2., (TOP + h - BOTTOM) / 2.)
    }
    fn to_world(&self, p: Vec2, w: f32, h: f32) -> Vec2 {
        let origin = Self::view_origin(w, h);
        origin + (p - origin - self.pan) / self.zoom
    }
    fn to_screen(&self, p: Vec2, w: f32, h: f32) -> Vec2 {
        let origin = Self::view_origin(w, h);
        origin + (p - origin) * self.zoom + self.pan
    }
    fn adjust_zoom(&mut self, delta: f32, anchor: Vec2, w: f32, h: f32) {
        let before = self.to_world(anchor, w, h);
        self.zoom = (self.zoom * delta).clamp(0.2, 4.0);
        self.preferred_zoom = self.zoom;
        self.zoom_preference_dirty = true;
        let origin = Self::view_origin(w, h);
        self.pan = anchor - origin - (before - origin) * self.zoom;
        self.pan.x = self.pan.x.clamp(-w * 2., w * 2.);
        self.pan.y = self.pan.y.clamp(-h * 2., h * 2.);
    }
    fn checkpoint(&mut self) {
        if self.lesson.is_some() && !self.lesson_open {
            self.lesson = None;
            self.naming_highlight.clear();
        }
        self.reaction_highlight_bonds.clear();
        self.undo.push(Snapshot {
            mol: self.mol.clone(),
            graphics: self.graphics.clone(),
        });
        if self.undo.len() > 100 {
            self.undo.remove(0);
        }
        self.redo.clear();
        self.dirty = true;
    }
    fn undo(&mut self) {
        if self.lesson.is_some() && !self.lesson_open {
            self.lesson = None;
            self.naming_highlight.clear();
        }
        if let Some(m) = self.undo.pop() {
            self.redo.push(Snapshot {
                mol: self.mol.clone(),
                graphics: self.graphics.clone(),
            });
            self.mol = m.mol;
            self.graphics = m.graphics;
            self.dirty = true;
            self.selected = None;
            self.selected_atoms.clear();
            self.selected_bond = None;
            self.selected_graphic = None;
            self.status = "Undid edit".into();
        }
    }
    fn redo(&mut self) {
        if self.lesson.is_some() && !self.lesson_open {
            self.lesson = None;
            self.naming_highlight.clear();
        }
        if let Some(m) = self.redo.pop() {
            self.undo.push(Snapshot {
                mol: self.mol.clone(),
                graphics: self.graphics.clone(),
            });
            self.mol = m.mol;
            self.graphics = m.graphics;
            self.dirty = true;
            self.selected = None;
            self.selected_atoms.clear();
            self.selected_bond = None;
            self.selected_graphic = None;
            self.status = "Redid edit".into();
        }
    }
    fn hit_atom(&self, p: Vec2) -> Option<usize> {
        self.mol
            .atoms
            .iter()
            .enumerate()
            .rev()
            .find(|(_, a)| a.pos.distance(p) < 15.)
            .map(|(i, _)| i)
    }
    fn hit_bond(&self, p: Vec2) -> Option<usize> {
        self.mol
            .bonds
            .iter()
            .enumerate()
            .find(|(_, b)| {
                distance_segment(p, self.mol.atoms[b.a].pos, self.mol.atoms[b.b].pos) < 9.
            })
            .map(|(i, _)| i)
    }
    fn atoms_in_rect(&self, a: Vec2, b: Vec2) -> Vec<usize> {
        let (left, right) = (a.x.min(b.x), a.x.max(b.x));
        let (top, bottom) = (a.y.min(b.y), a.y.max(b.y));
        self.mol
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(i, atom)| {
                (atom.pos.x >= left
                    && atom.pos.x <= right
                    && atom.pos.y >= top
                    && atom.pos.y <= bottom)
                    .then_some(i)
            })
            .collect()
    }
    fn select_at(&mut self, p: Vec2, additive: bool) {
        let atom = self.hit_atom(p);
        let bond = atom.is_none().then(|| self.hit_bond(p)).flatten();
        let graphic = (atom.is_none() && bond.is_none())
            .then(|| self.hit_graphic(p))
            .flatten();
        if !additive {
            self.selected_atoms.clear();
            self.selected_bonds.clear();
            self.selected_graphics.clear();
        }
        if let Some(i) = atom {
            if !self.selected_atoms.contains(&i) {
                self.selected_atoms.push(i);
            }
        } else if let Some(i) = bond {
            if !self.selected_bonds.contains(&i) {
                self.selected_bonds.push(i);
            }
        } else if let Some(i) = graphic {
            if !self.selected_graphics.contains(&i) {
                self.selected_graphics.push(i);
            }
        }
        self.selected = (self.selected_atoms.len() == 1).then(|| self.selected_atoms[0]);
        self.selected_bond = self.selected_bonds.first().copied();
        self.selected_graphic = self.selected_graphics.first().copied();
        self.selection_box = if atom.is_none() && bond.is_none() && graphic.is_none() {
            Some(p)
        } else {
            None
        };
    }
    fn select_box(&mut self, a: Vec2, b: Vec2) {
        self.selected_atoms = self.atoms_in_rect(a, b);
        let inside = |p: Vec2| {
            p.x >= a.x.min(b.x) && p.x <= a.x.max(b.x) && p.y >= a.y.min(b.y) && p.y <= a.y.max(b.y)
        };
        self.selected_bonds = self
            .mol
            .bonds
            .iter()
            .enumerate()
            .filter_map(|(i, bond)| {
                let p = (self.mol.atoms[bond.a].pos + self.mol.atoms[bond.b].pos) / 2.;
                inside(p).then_some(i)
            })
            .collect();
        self.selected_graphics = self
            .graphics
            .iter()
            .enumerate()
            .filter_map(|(i, graphic)| {
                let selected = match graphic {
                    Graphic::Arrow { start, end, .. } => {
                        inside(*start) || inside(*end) || inside((*start + *end) / 2.)
                    }
                    Graphic::ElectronArrow {
                        start,
                        control,
                        end,
                        ..
                    } => inside(*start) || inside(*control) || inside(*end),
                    Graphic::Text { pos, .. } | Graphic::Plus { pos } => inside(*pos),
                    Graphic::Shape { start, end, .. } => {
                        inside(*start) || inside(*end) || inside((*start + *end) / 2.)
                    }
                };
                selected.then_some(i)
            })
            .collect();
        self.selected = (self.selected_atoms.len() == 1).then(|| self.selected_atoms[0]);
        self.selected_bond = self.selected_bonds.first().copied();
        self.selected_graphic = self.selected_graphics.first().copied();
        self.status = if self.selected_atoms.is_empty()
            && self.selected_bonds.is_empty()
            && self.selected_graphics.is_empty()
        {
            "Nothing selected".into()
        } else {
            format!(
                "Selected {} atom(s), {} bond(s), {} object(s)",
                self.selected_atoms.len(),
                self.selected_bonds.len(),
                self.selected_graphics.len()
            )
        };
    }
    fn has_selection(&self) -> bool {
        !self.selected_atoms.is_empty()
            || !self.selected_bonds.is_empty()
            || !self.selected_graphics.is_empty()
    }
    fn move_selection(&mut self, delta: Vec2) {
        let mut moved_atoms = self.selected_atoms.clone();
        for &bi in &self.selected_bonds {
            if let Some(b) = self.mol.bonds.get(bi) {
                moved_atoms.push(b.a);
                moved_atoms.push(b.b);
            }
        }
        moved_atoms.sort_unstable();
        moved_atoms.dedup();
        for &i in &moved_atoms {
            if let Some(atom) = self.mol.atoms.get_mut(i) {
                atom.pos += delta;
            }
        }
        for &i in &self.selected_graphics {
            if let Some(g) = self.graphics.get_mut(i) {
                match g {
                    Graphic::Arrow { start, end, .. } => {
                        *start += delta;
                        *end += delta;
                    }
                    Graphic::ElectronArrow {
                        start,
                        control,
                        end,
                        ..
                    } => {
                        *start += delta;
                        *control += delta;
                        *end += delta;
                    }
                    Graphic::Text { pos, .. } | Graphic::Plus { pos } => *pos += delta,
                    Graphic::Shape { start, end, .. } => {
                        *start += delta;
                        *end += delta;
                    }
                }
            }
        }
    }
    fn delete_selection(&mut self) {
        if !self.has_selection() {
            return;
        }
        self.checkpoint();
        let atom_ids: std::collections::BTreeSet<usize> =
            self.selected_atoms.iter().copied().collect();
        let bond_ids: std::collections::BTreeSet<usize> =
            self.selected_bonds.iter().copied().collect();
        let graphic_ids: std::collections::BTreeSet<usize> =
            self.selected_graphics.iter().copied().collect();
        let mut remap = vec![None; self.mol.atoms.len()];
        let mut atoms = Vec::new();
        for (i, atom) in self.mol.atoms.drain(..).enumerate() {
            if !atom_ids.contains(&i) {
                remap[i] = Some(atoms.len());
                atoms.push(atom);
            }
        }
        self.mol.atoms = atoms;
        let mut new_bonds = Vec::new();
        for (i, mut bond) in self.mol.bonds.drain(..).enumerate() {
            if bond_ids.contains(&i) || atom_ids.contains(&bond.a) || atom_ids.contains(&bond.b) {
                continue;
            }
            let (Some(a), Some(b)) = (remap[bond.a], remap[bond.b]) else {
                continue;
            };
            bond.a = a;
            bond.b = b;
            new_bonds.push(bond);
        }
        self.mol.bonds = new_bonds;
        self.graphics = self
            .graphics
            .drain(..)
            .enumerate()
            .filter_map(|(i, g)| (!graphic_ids.contains(&i)).then_some(g))
            .collect();
        self.clear_selection();
        self.status = "Deleted selection".into();
    }
    fn clear_selection(&mut self) {
        self.selected = None;
        self.selected_atoms.clear();
        self.selected_bond = None;
        self.selected_bonds.clear();
        self.selected_graphic = None;
        self.selected_graphics.clear();
    }
    fn select_all(&mut self) {
        self.selected_atoms = (0..self.mol.atoms.len()).collect();
        self.selected_bonds = (0..self.mol.bonds.len()).collect();
        self.selected_graphics = (0..self.graphics.len()).collect();
        self.selected = None;
        self.selected_bond = self.selected_bonds.first().copied();
        self.selected_graphic = self.selected_graphics.first().copied();
        self.status = format!(
            "Selected {} atoms, {} bonds and {} objects",
            self.selected_atoms.len(),
            self.selected_bonds.len(),
            self.selected_graphics.len()
        );
    }
    fn set_reaction_component(
        &mut self,
        role: Option<&str>,
        step: Option<u16>,
        coefficient: Option<u8>,
    ) {
        let mut seeds = self.selected_atoms.clone();
        for &i in &self.selected_bonds {
            if let Some(b) = self.mol.bonds.get(i) {
                seeds.push(b.a);
                seeds.push(b.b);
            }
        }
        seeds.sort_unstable();
        seeds.dedup();
        let Some(root) = seeds.first().copied() else {
            self.status = "Select an atom or bond in the reaction component".into();
            return;
        };
        let mut seen = vec![false; self.mol.atoms.len()];
        let mut stack = vec![root];
        seen[root] = true;
        while let Some(i) = stack.pop() {
            for b in &self.mol.bonds {
                let j = if b.a == i {
                    Some(b.b)
                } else if b.b == i {
                    Some(b.a)
                } else {
                    None
                };
                if let Some(j) = j {
                    if !seen[j] {
                        seen[j] = true;
                        stack.push(j);
                    }
                }
            }
        }
        self.checkpoint();
        for (i, atom) in self.mol.atoms.iter_mut().enumerate() {
            if seen[i] {
                if let Some(role) = role {
                    atom.reaction_role = if role.is_empty() {
                        None
                    } else {
                        Some(role.into())
                    };
                }
                if let Some(step) = step {
                    atom.reaction_step = step.max(1);
                }
                if let Some(coefficient) = coefficient {
                    atom.coefficient = coefficient.clamp(1, 99);
                }
            }
        }
        self.status = "Updated reaction component metadata".into();
    }
    fn assign_unique_atom_maps(&mut self) {
        let mut targets = self.selected_atoms.clone();
        if let Some(atom) = self.selected {
            targets.push(atom);
        }
        targets.retain(|&index| index < self.mol.atoms.len());
        targets.sort_unstable();
        targets.dedup();
        if targets.is_empty() {
            targets.extend(0..self.mol.atoms.len());
        }
        if targets.is_empty() {
            self.status = "Draw or select atoms before assigning map numbers".into();
            return;
        }
        let target_set: BTreeSet<usize> = targets.iter().copied().collect();
        let used: BTreeSet<u16> = self
            .mol
            .atoms
            .iter()
            .enumerate()
            .filter(|(index, _)| !target_set.contains(index))
            .filter_map(|(_, atom)| atom.atom_map)
            .collect();
        let assignments: Vec<u16> = (1..=32_767)
            .filter(|number| !used.contains(number))
            .take(targets.len())
            .collect();
        if assignments.len() != targets.len() {
            self.status = "No free atom-map numbers remain in the supported range 1–32767".into();
            return;
        }
        if targets
            .iter()
            .zip(&assignments)
            .all(|(&index, &number)| self.mol.atoms[index].atom_map == Some(number))
        {
            self.status = "Selected atoms already have these unique map numbers".into();
            return;
        }
        self.checkpoint();
        for (index, number) in targets.into_iter().zip(assignments) {
            self.mol.atoms[index].atom_map = Some(number);
        }
        self.status = "Assigned unique atom-map numbers; undo restores previous maps".into();
    }
    fn duplicate_selection(&mut self) {
        if !self.has_selection() {
            self.status = "Select atoms, bonds or objects to duplicate".into();
            return;
        }
        self.checkpoint();
        let mut selected_atoms = self.selected_atoms.clone();
        for &bi in &self.selected_bonds {
            if let Some(b) = self.mol.bonds.get(bi) {
                selected_atoms.push(b.a);
                selected_atoms.push(b.b);
            }
        }
        selected_atoms.sort_unstable();
        selected_atoms.dedup();
        let atom_set: std::collections::BTreeSet<usize> = selected_atoms.iter().copied().collect();
        let offset = vec2(BOND * 1.1, BOND * 0.8);
        let mut remap = std::collections::BTreeMap::new();
        let old_len = self.mol.atoms.len();
        for i in selected_atoms {
            let mut atom = self.mol.atoms[i].clone();
            atom.pos += offset;
            atom.atom_map = None;
            remap.insert(i, self.mol.atoms.len());
            self.mol.atoms.push(atom);
        }
        let selected_bonds: std::collections::BTreeSet<usize> =
            self.selected_bonds.iter().copied().collect();
        let old_bond_len = self.mol.bonds.len();
        let bonds: Vec<Bond> = self
            .mol
            .bonds
            .iter()
            .enumerate()
            .filter_map(|(i, b)| {
                let chosen = (atom_set.contains(&b.a) && atom_set.contains(&b.b))
                    || selected_bonds.contains(&i);
                if !chosen {
                    return None;
                }
                let (&a, &c) = (remap.get(&b.a)?, remap.get(&b.b)?);
                let mut copy = b.clone();
                copy.a = a;
                copy.b = c;
                Some(copy)
            })
            .collect();
        self.mol.bonds.extend(bonds);
        let new_bond_len = self.mol.bonds.len();
        let selected_graphics = self.selected_graphics.clone();
        for i in selected_graphics {
            if let Some(mut g) = self.graphics.get(i).cloned() {
                match &mut g {
                    Graphic::Arrow { start, end, .. } => {
                        *start += offset;
                        *end += offset;
                    }
                    Graphic::ElectronArrow {
                        start,
                        control,
                        end,
                        ..
                    } => {
                        *start += offset;
                        *control += offset;
                        *end += offset;
                    }
                    Graphic::Text { pos, .. } | Graphic::Plus { pos } => *pos += offset,
                    Graphic::Shape { start, end, .. } => {
                        *start += offset;
                        *end += offset;
                    }
                }
                self.graphics.push(g);
            }
        }
        self.selected_atoms = (old_len..self.mol.atoms.len()).collect();
        self.selected_bonds = (old_bond_len..new_bond_len).collect();
        self.selected_graphics.clear();
        self.selected = None;
        self.selected_bond = self.selected_bonds.first().copied();
        self.selected_graphic = None;
        self.status = "Duplicated selection".into();
    }
    fn rotate_selection(&mut self, angle: f32) {
        if self.selected_atoms.is_empty()
            && self.selected_bonds.is_empty()
            && self.selected_graphics.is_empty()
        {
            self.status = "Select atoms or objects to rotate".into();
            return;
        }
        self.checkpoint();
        let mut atoms_to_rotate = self.selected_atoms.clone();
        for &bi in &self.selected_bonds {
            if let Some(b) = self.mol.bonds.get(bi) {
                atoms_to_rotate.push(b.a);
                atoms_to_rotate.push(b.b);
            }
        }
        atoms_to_rotate.sort_unstable();
        atoms_to_rotate.dedup();
        let mut points: Vec<Vec2> = atoms_to_rotate
            .iter()
            .filter_map(|&i| self.mol.atoms.get(i).map(|a| a.pos))
            .collect();
        for &i in &self.selected_graphics {
            if let Some(g) = self.graphics.get(i) {
                points.push(match g {
                    Graphic::Arrow { start, end, .. } | Graphic::Shape { start, end, .. } => {
                        (*start + *end) / 2.
                    }
                    Graphic::ElectronArrow { control, .. } => *control,
                    Graphic::Text { pos, .. } | Graphic::Plus { pos } => *pos,
                });
            }
        }
        let center = points.iter().copied().sum::<Vec2>() / points.len() as f32;
        let (s, c) = angle.sin_cos();
        let rot = |p: Vec2| {
            let q = p - center;
            center + vec2(q.x * c - q.y * s, q.x * s + q.y * c)
        };
        for &i in &atoms_to_rotate {
            if let Some(a) = self.mol.atoms.get_mut(i) {
                a.pos = rot(a.pos);
            }
        }
        for &i in &self.selected_graphics {
            if let Some(g) = self.graphics.get_mut(i) {
                match g {
                    Graphic::Arrow { start, end, .. } | Graphic::Shape { start, end, .. } => {
                        *start = rot(*start);
                        *end = rot(*end);
                    }
                    Graphic::ElectronArrow {
                        start,
                        control,
                        end,
                        ..
                    } => {
                        *start = rot(*start);
                        *control = rot(*control);
                        *end = rot(*end);
                    }
                    Graphic::Text { pos, .. } | Graphic::Plus { pos } => *pos = rot(*pos),
                }
            }
        }
        self.status = "Rotated selected atoms".into();
    }
    fn display_snapshot(&self, w: f32, h: f32) -> Snapshot {
        let mut snapshot = Snapshot {
            mol: self.mol.clone(),
            graphics: self.graphics.clone(),
        };
        for atom in &mut snapshot.mol.atoms {
            atom.pos = self.to_screen(atom.pos, w, h);
        }
        for graphic in &mut snapshot.graphics {
            match graphic {
                Graphic::Arrow { start, end, .. } => {
                    *start = self.to_screen(*start, w, h);
                    *end = self.to_screen(*end, w, h);
                }
                Graphic::ElectronArrow {
                    start,
                    control,
                    end,
                    ..
                } => {
                    *start = self.to_screen(*start, w, h);
                    *control = self.to_screen(*control, w, h);
                    *end = self.to_screen(*end, w, h);
                }
                Graphic::Text { pos, .. } | Graphic::Plus { pos } => {
                    *pos = self.to_screen(*pos, w, h)
                }
                Graphic::Shape { start, end, .. } => {
                    *start = self.to_screen(*start, w, h);
                    *end = self.to_screen(*end, w, h);
                }
            }
        }
        snapshot
    }
    fn sync_active_document(&mut self) {
        if let Some(doc) = self.documents.get_mut(self.active_document) {
            doc.snapshot = Snapshot {
                mol: self.mol.clone(),
                graphics: self.graphics.clone(),
            };
            doc.dirty = self.dirty;
            doc.undo = self.undo.clone();
            doc.redo = self.redo.clone();
        }
    }
    fn switch_document(&mut self, index: usize) {
        if index >= self.documents.len() || index == self.active_document {
            return;
        }
        self.sync_active_document();
        self.lesson = None;
        self.lesson_open = false;
        self.naming_highlight.clear();
        self.active_document = index;
        let (mol, graphics, dirty, title, undo, redo) = {
            let doc = &self.documents[index];
            (
                doc.snapshot.mol.clone(),
                doc.snapshot.graphics.clone(),
                doc.dirty,
                doc.title.clone(),
                doc.undo.clone(),
                doc.redo.clone(),
            )
        };
        self.mol = mol;
        self.graphics = graphics;
        self.dirty = dirty;
        self.undo = undo;
        self.redo = redo;
        self.clear_selection();
        self.zoom = self.preferred_zoom;
        self.pan = Vec2::ZERO;
        self.status = format!("Opened tab: {title}");
    }
    fn request_close_document(&mut self, index: usize) {
        if index >= self.documents.len() {
            return;
        }
        if self.documents[index].dirty {
            self.pending_close = Some(index);
        } else {
            self.close_document(index);
        }
    }
    fn close_document(&mut self, index: usize) {
        if index >= self.documents.len() {
            return;
        }
        self.sync_active_document();
        let previous = self.active_document;
        self.documents.remove(index);
        if self.documents.is_empty() {
            self.documents.push(DocumentTab::blank());
            self.active_document = 0;
        } else if index < previous {
            self.active_document = previous - 1;
        } else if index == previous {
            self.active_document = index.min(self.documents.len() - 1);
        } else {
            self.active_document = previous.min(self.documents.len() - 1);
        }
        let doc = &self.documents[self.active_document];
        self.mol = doc.snapshot.mol.clone();
        self.graphics = doc.snapshot.graphics.clone();
        self.dirty = doc.dirty;
        self.undo = doc.undo.clone();
        self.redo = doc.redo.clone();
        self.clear_selection();
        self.zoom = self.preferred_zoom;
        self.pan = Vec2::ZERO;
        self.pending_close = None;
        self.status = "Closed document tab".into();
    }
    fn close_pending_with_save(&mut self) {
        let Some(index) = self.pending_close else {
            return;
        };
        let Some(path) = self.documents.get(index).and_then(|doc| doc.path.clone()) else {
            self.switch_document(index);
            self.close_after_save = Some(index);
            self.pending_close = None;
            self.file_dialog = Some(true);
            return;
        };
        self.switch_document(index);
        self.pending_close = None;
        self.close_after_save = Some(index);
        self.save(&path);
        if self.close_after_save == Some(index) {
            self.close_after_save = None;
        }
    }
    fn save_before_quit(&mut self) {
        self.pending_quit = false;
        self.quit_when_clean = true;
        self.save_all();
        if let Some(index) = self
            .documents
            .iter()
            .position(|doc| doc.dirty && doc.path.is_none())
        {
            self.switch_document(index);
            self.file_dialog = Some(true);
        }
        if !self.documents.iter().any(|doc| doc.dirty) {
            self.quit_now = true;
        }
    }
    fn analyze_reaction(&mut self) {
        match reaction::analyze_all(&self.mol, &self.graphics) {
            Ok(mut reports) => {
                if let Some(step) = self.reaction_step_view {
                    reports.retain(|report| report.step == step);
                    if reports.is_empty() {
                        self.reaction_highlight_bonds.clear();
                        self.info_overlay = Some(format!(
                            "Reaction step {step} does not exist. Draw another reaction arrow or show all steps."
                        ));
                        return;
                    }
                }
                let pairs: std::collections::BTreeSet<(u16, u16)> = reports
                    .iter()
                    .flat_map(|r| r.changed_map_pairs.iter().copied())
                    .collect();
                self.reaction_highlight_bonds = self
                    .mol
                    .bonds
                    .iter()
                    .enumerate()
                    .filter_map(|(i, b)| {
                        let a = self.mol.atoms.get(b.a)?.atom_map?;
                        let c = self.mol.atoms.get(b.b)?.atom_map?;
                        pairs.contains(&(a.min(c), a.max(c))).then_some(i)
                    })
                    .collect();
                self.info_overlay = Some(reaction::report_all_text(&reports));
            }
            Err(e) => {
                self.reaction_highlight_bonds.clear();
                self.info_overlay = Some(e);
            }
        }
    }
    fn auto_balance_reaction(&mut self) {
        let step = self.reaction_step_view.unwrap_or(1);
        match reaction::propose_balancing_coefficients(&self.mol, &self.graphics, step) {
            Ok(proposals) => {
                let changed = proposals.iter().any(|proposal| {
                    proposal.atom_indices.iter().any(|&index| {
                        self.mol
                            .atoms
                            .get(index)
                            .is_some_and(|atom| atom.coefficient != proposal.coefficient)
                    })
                });
                if !changed {
                    self.status = format!(
                        "Reaction step {step} already has the smallest balanced coefficients"
                    );
                    return;
                }
                self.checkpoint();
                let mut lines = vec![format!(
                    "Balanced reaction step {step} by elements and formal charge:"
                )];
                for proposal in proposals {
                    for &index in &proposal.atom_indices {
                        if let Some(atom) = self.mol.atoms.get_mut(index) {
                            atom.coefficient = proposal.coefficient;
                        }
                    }
                    lines.push(format!(
                        "{}: coefficient {}",
                        proposal.formula, proposal.coefficient
                    ));
                }
                self.info_scroll = 0;
                self.info_overlay = Some(lines.join("\n"));
                self.status =
                    format!("Reaction step {step} balanced; undo restores previous coefficients");
            }
            Err(error) => {
                self.info_scroll = 0;
                self.info_overlay = Some(format!("Automatic balancing unavailable\n{error}"));
            }
        }
    }
    fn auto_balance_all_reactions(&mut self) {
        let step_count = self
            .graphics
            .iter()
            .filter(|graphic| matches!(graphic, Graphic::Arrow { kind, .. } if matches!(kind.as_str(), "reaction" | "equilibrium" | "retrosynthesis")))
            .count();
        if step_count == 0 {
            self.info_overlay =
                Some("Draw at least one reaction arrow before balancing all steps.".into());
            return;
        }

        // Solve every step before editing anything. If one step is ambiguous or
        // incomplete, preserve all current coefficients and explain the step.
        let mut by_step = Vec::with_capacity(step_count);
        for step in 1..=step_count {
            match reaction::propose_balancing_coefficients(&self.mol, &self.graphics, step as u16) {
                Ok(proposals) => by_step.push((step, proposals)),
                Err(error) => {
                    self.info_overlay = Some(format!(
                        "No coefficients were changed. Reaction step {step} could not be uniquely balanced.\n\n{error}"
                    ));
                    return;
                }
            }
        }

        let has_changes = by_step.iter().any(|(_, proposals)| {
            proposals.iter().any(|proposal| {
                proposal.atom_indices.iter().any(|&index| {
                    self.mol
                        .atoms
                        .get(index)
                        .is_some_and(|atom| atom.coefficient != proposal.coefficient)
                })
            })
        });
        if !has_changes {
            self.status =
                "All reaction steps already have the smallest balanced coefficients".into();
            return;
        }

        self.checkpoint();
        let mut lines = vec![format!(
            "Balanced all {step_count} reaction steps by elements and formal charge:"
        )];
        for (step, proposals) in by_step {
            lines.push(format!("\nStep {step}"));
            for proposal in proposals {
                for &index in &proposal.atom_indices {
                    if let Some(atom) = self.mol.atoms.get_mut(index) {
                        atom.coefficient = proposal.coefficient;
                    }
                }
                lines.push(format!(
                    "{}: coefficient {}",
                    proposal.formula, proposal.coefficient
                ));
            }
        }
        self.info_scroll = 0;
        self.info_overlay = Some(lines.join("\n"));
        self.status = format!("Balanced all {step_count} reaction steps");
    }
    fn start_naming_lesson(&mut self) {
        match naming::name(&self.mol) {
            Ok(result) => {
                self.naming_highlight = result.highlighted_atoms;
                self.lesson = Some(TeachingLesson { content: LessonContent::Naming { name: result.name, steps: result.steps }, step: 0 });
                self.lesson_open = true;
            }
            Err(error) => self.info_overlay = Some(format!("Naming lesson unavailable\n{error}\n\nDraw a supported acyclic carbon structure with the functional groups this teaching subset recognizes.")),
        }
    }
    fn start_reaction_lesson(&mut self) {
        self.lesson = Some(TeachingLesson {
            content: LessonContent::Reaction { feedback: None },
            step: 0,
        });
        self.lesson_open = true;
    }
    fn start_rdkit_graph_operation(&mut self, operation: &'static str) {
        if self.refuse_overlapping_worker() {
            return;
        }
        if self.mol.atoms.is_empty() {
            self.info_overlay = Some("Draw a structure first".into());
            return;
        }
        if self.mol.atoms.len() > 999 || self.mol.bonds.len() > 999 {
            self.info_overlay=Some("The optional RDKit editor bridge currently accepts MOL V2000 structures up to 999 atoms and bonds.".into());
            return;
        }
        if self.mol.atoms.iter().any(|a| {
            a.alias.is_some()
                || a.query.is_some()
                || a.r_group.is_some()
                || a.stereo.is_some()
                || a.atom_map.is_some()
                || a.radical > 0
        }) || self.mol.bonds.iter().any(|b| {
            b.query.is_some()
                || b.stereo.is_some()
                || matches!(
                    b.style,
                    BondStyle::Wedge
                        | BondStyle::Hash
                        | BondStyle::Wavy
                        | BondStyle::Coordinate
                        | BondStyle::Any
                )
        }) {
            self.info_overlay=Some("This optional conversion only accepts ordinary atom and bond graphs. Advanced query, map, alias, radical or stereochemical data would not be preserved, so the structure was left unchanged.".into());
            return;
        }
        let block = self.molfile();
        let captured = self.mol.clone();
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || {
            let result = match operation {
                "aromatize" => rdkit::aromatize(&block),
                _ => rdkit::kekulize(&block),
            };
            let _ = tx.send(result);
        });
        self.calc_rx = Some((captured, rx));
        self.status = format!("Running optional RDKit {operation} operation…");
    }
    fn start_stereo_report(&mut self) {
        self.start_stereo_check(false);
    }
    fn start_stereo_assignment(&mut self) {
        self.start_stereo_check(true);
    }
    fn start_copy_smiles(&mut self) {
        if self.refuse_overlapping_worker() {
            return;
        }
        if self.mol.atoms.is_empty() {
            self.info_overlay = Some("Draw or open a molecule before copying SMILES".into());
            return;
        }
        if self.mol.atoms.iter().any(|atom| {
            atom.alias.is_some()
                || atom.r_group.is_some()
                || atom.query.is_some()
                || !chemistry::valid_element(&atom.element)
        }) || self.mol.bonds.iter().any(|bond| {
            bond.query.is_some()
                || matches!(
                    bond.style,
                    BondStyle::Coordinate | BondStyle::Any | BondStyle::Wavy
                )
        }) {
            self.info_overlay = Some(
                "This structure contains aliases, query features or nonstandard bonds that cannot be represented safely in SMILES. The clipboard was left unchanged.".into(),
            );
            return;
        }
        let block = self.molfile();
        let captured = self.mol.clone();
        let smiles_molecule = captured.clone();
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || {
            let result = smiles::export(&smiles_molecule)
                .or_else(|_| rdkit::to_smiles_with_stereo(&block, &smiles_molecule))
                .and_then(|smiles| {
                    copy_text_to_clipboard(&smiles)?;
                    Ok(format!("COPIED_SMILES\n{smiles}"))
                });
            let _ = tx.send(result);
        });
        self.calc_rx = Some((captured, rx));
        self.status =
            "Creating SMILES locally when supported (optional RDKit fallback otherwise)…".into();
    }
    fn start_copy_smarts(&mut self) {
        match smarts::export(&self.mol) {
            Ok(pattern) => match copy_text_to_clipboard(&pattern) {
                Ok(()) => {
                    self.info_overlay = Some(format!(
                        "SMARTS copied (supported local subset)\n\n{pattern}\n\nMolecule metadata and drawing objects are not part of SMARTS."
                    ));
                }
                Err(error) => self.info_overlay = Some(error),
            },
            Err(error) => {
                self.info_overlay = Some(format!(
                    "SMARTS export unavailable for this structure\n{error}\n\nThe clipboard was left unchanged."
                ));
            }
        }
    }
    fn start_stereo_check(&mut self, apply_labels: bool) {
        if self.refuse_overlapping_worker() {
            return;
        }
        if self.mol.atoms.is_empty() {
            self.info_overlay = Some("Draw a structure first".into());
            return;
        }
        if self.mol.atoms.len() > 999
            || self.mol.bonds.len() > 999
            || self.mol.atoms.iter().any(|a| {
                a.query.is_some()
                    || a.r_group.is_some()
                    || a.alias.is_some()
                    || !chemistry::valid_element(&a.element)
            })
            || self.mol.bonds.iter().any(|b| b.query.is_some())
        {
            self.info_overlay=Some("The optional RDKit stereochemistry check accepts ordinary MOL V2000 structures only. No report was calculated for this structure.".into());
            return;
        }
        let block = self.molfile();
        let captured = self.mol.clone();
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || {
            let result = rdkit::stereo_report(&block).map(|report| {
                if apply_labels {
                    format!("ASSIGN_STEREO\n{report}")
                } else {
                    report
                }
            });
            let _ = tx.send(result);
        });
        self.calc_rx = Some((captured, rx));
        self.info_overlay = Some(if apply_labels {
            "Assigning CIP R/S labels with optional RDKit…".into()
        } else {
            "Checking stereochemistry with optional RDKit…".into()
        });
    }
    fn apply_stereo_labels(&mut self, report: &str) {
        let mut labels = BTreeMap::new();
        for line in report.lines() {
            let Some(rest) = line.strip_prefix("ATOM\t") else {
                continue;
            };
            let mut fields = rest.split('\t');
            let Ok(index) = fields.next().unwrap_or_default().parse::<usize>() else {
                self.info_overlay =
                    Some("RDKit returned an invalid atom index; no labels were applied.".into());
                self.status = "Rejected invalid CIP-label output; structure unchanged".into();
                return;
            };
            let label = fields.next().unwrap_or_default();
            if index == 0
                || index > self.mol.atoms.len()
                || !matches!(label, "R" | "S")
                || fields.next().is_some()
                || labels.insert(index - 1, label).is_some()
            {
                self.info_overlay = Some(
                    "RDKit returned invalid or duplicate CIP labels; no labels were applied."
                        .into(),
                );
                self.status = "Rejected invalid CIP-label output; structure unchanged".into();
                return;
            }
        }
        let changes = self.mol.atoms.iter().enumerate().any(|(index, atom)| {
            let current = atom
                .stereo
                .as_deref()
                .filter(|value| matches!(*value, "R" | "S"));
            current != labels.get(&index).copied()
        });
        if changes {
            self.checkpoint();
            for (index, atom) in self.mol.atoms.iter_mut().enumerate() {
                if matches!(atom.stereo.as_deref(), Some("R" | "S")) {
                    atom.stereo = None;
                }
                if let Some(label) = labels.get(&index) {
                    atom.stereo = Some((*label).to_string());
                }
            }
        }
        let count = labels.len();
        self.info_overlay = Some(format!(
            "{}\n\nApplied {count} CIP label(s) from RDKit. Drawing wedge/hash marks define the stereochemistry; labels are undoable annotations.",
            report.trim()
        ));
        self.info_scroll = 0;
        self.status = format!("Applied {count} calculated CIP label(s)");
    }
    fn start_rdkit_chemical_report(&mut self, operation: &'static str) {
        if self.refuse_overlapping_worker() {
            return;
        }
        if self.mol.atoms.is_empty() {
            self.info_overlay = Some("Draw a structure first".into());
            return;
        }
        if self.mol.atoms.len() > 999
            || self.mol.bonds.len() > 999
            || self.mol.atoms.iter().any(|atom| {
                atom.alias.is_some()
                    || atom.query.is_some()
                    || atom.r_group.is_some()
                    || !chemistry::valid_element(&atom.element)
            })
            || self.mol.bonds.iter().any(|bond| bond.query.is_some())
        {
            self.info_overlay = Some("This optional RDKit analysis accepts ordinary atom and bond graphs up to the MOL V2000 limits. No candidates were calculated.".into());
            return;
        }
        let block = self.molfile();
        let captured = self.mol.clone();
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || {
            let _ = tx.send(rdkit::import_external(operation, &block));
        });
        self.calc_rx = Some((captured, rx));
        self.info_overlay = Some(format!("Running optional RDKit {operation} analysis…"));
        self.info_scroll = 0;
    }
    fn start_rdkit_partial_charges(&mut self) {
        if self.refuse_overlapping_worker() {
            return;
        }
        if self.mol.atoms.is_empty() {
            self.info_overlay = Some("Draw a structure before calculating partial charges".into());
            return;
        }
        if self.mol.atoms.len() > 999
            || self.mol.bonds.len() > 999
            || self.mol.atoms.iter().any(|atom| {
                atom.alias.is_some()
                    || atom.query.is_some()
                    || atom.r_group.is_some()
                    || !chemistry::valid_element(&atom.element)
            })
            || self.mol.bonds.iter().any(|bond| bond.query.is_some())
        {
            self.info_overlay = Some("The optional RDKit Gasteiger calculation accepts ordinary MOL V2000 atom/bond graphs only. Query atoms, aliases and R-groups are refused; the drawing was not changed.".into());
            return;
        }
        let block = self.molfile();
        let captured = self.mol.clone();
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || {
            let _ = tx.send(rdkit::partial_charge_report(&block));
        });
        self.calc_rx = Some((captured, rx));
        self.info_overlay =
            Some("Calculating Gasteiger-Marsili partial charges with optional RDKit…".into());
        self.info_scroll = 0;
        self.status = "Calculating modelled partial charges…".into();
    }
    fn start_pubchem_lookup(&mut self) {
        if self.refuse_overlapping_worker() {
            return;
        }
        if self.mol.atoms.is_empty() {
            self.info_overlay =
                Some("Draw or open a structure before looking it up in PubChem".into());
            return;
        }
        let smiles = match smiles::export(&self.mol) {
            Ok(smiles) => smiles,
            Err(error) => {
                self.info_overlay = Some(format!(
                    "PubChem lookup cannot safely identify this drawing as a complete structure. No data was sent.\n\n{error}"
                ));
                return;
            }
        };
        let captured = self.mol.clone();
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || {
            let result =
                pubchem::lookup(&smiles).map_err(|message| format!("PUBCHEM_LOOKUP\n{message}"));
            let _ = tx.send(result);
        });
        self.calc_rx = Some((captured, rx));
        self.info_overlay = Some(
            "Sending the current structure to PubChem through the official NCBI PUG REST API. The drawing will not be changed.".into(),
        );
        self.info_scroll = 0;
        self.status = "Fetching matching PubChem database records…".into();
    }
    fn apply_rdkit_molblock(&mut self, metadata: &str, block: &str) {
        let operation = metadata
            .lines()
            .find_map(|line| line.strip_prefix("OPERATION\t"));
        let generated = match parse_mol_at(block, vec2(0., 0.)) {
            Ok(m) => m,
            Err(error) => {
                self.info_overlay =
                    Some(format!("RDKit returned an unreadable structure: {error}"));
                return;
            }
        };
        if generated.atoms.len() != self.mol.atoms.len()
            || generated.bonds.len() != self.mol.bonds.len()
        {
            self.info_overlay=Some("RDKit returned a graph with a different atom or bond count. The source document was left unchanged.".into());
            return;
        }
        if let Some(operation) = operation {
            self.checkpoint();
            for (bond, converted) in self.mol.bonds.iter_mut().zip(&generated.bonds) {
                bond.order = converted.order;
                bond.style = converted.style;
            }
            let aromatic_atoms: BTreeSet<usize> = self
                .mol
                .bonds
                .iter()
                .filter(|b| b.order == 4)
                .flat_map(|b| [b.a, b.b])
                .collect();
            for (i, atom) in self.mol.atoms.iter_mut().enumerate() {
                atom.aromatic = aromatic_atoms.contains(&i);
            }
            let engine = metadata
                .lines()
                .find(|l| l.starts_with("ENGINE\t"))
                .unwrap_or("ENGINE\tRDKit")
                .trim_start_matches("ENGINE\t");
            self.info_overlay=Some(format!("{} completed with {engine}. Atom and bond counts were preserved; this operation can be undone.",if operation=="aromatize"{"Aromatization"}else{"Kekulé conversion"}));
        } else {
            self.checkpoint();
            for (atom, coords) in self.mol.atoms.iter_mut().zip(generated.atoms) {
                atom.coord3d = coords.coord3d;
            }
            let engine = metadata
                .lines()
                .find(|line| line.starts_with("ENGINE\t"))
                .unwrap_or("ENGINE\tRDKit");
            let method = metadata
                .lines()
                .find(|line| line.starts_with("METHOD\t"))
                .unwrap_or("METHOD\tunknown force field");
            self.info_overlay=Some(format!("3D coordinates generated\n{}\n{}\nCoordinates stored separately from the 2D drawing.",engine.replace("ENGINE\t",""),method.replace("METHOD\t","")));
        }
    }
    fn check_reaction_lesson(&mut self) {
        let result = match reaction::analyze_all(&self.mol, &self.graphics) {
            Ok(reports) if reports.iter().all(|r|r.balanced==Some(true)) => format!("Correct: all analyzed reaction steps balance by formula and formal charge.\n\n{}",reaction::report_all_text(&reports)),
            Ok(reports) if reports.iter().any(|r|r.balanced==Some(false)) => format!("Not balanced yet. Compare the reactant/product element and charge totals below, then adjust component roles or coefficients.\n\n{}",reaction::report_all_text(&reports)),
            Ok(reports) => format!("Not evaluated: at least one step has incomplete or unsupported chemistry, so no pass/fail result is assigned.\n\n{}",reaction::report_all_text(&reports)),
            Err(error) => format!("Not evaluated: {error}. Add a reaction arrow and at least one reactant and one product, then try again."),
        };
        if let Some(TeachingLesson {
            content: LessonContent::Reaction { feedback },
            ..
        }) = &mut self.lesson
        {
            *feedback = Some(result);
        }
    }
    fn discard_and_quit(&mut self) {
        let _ = fs::remove_file(autosave_path());
        self.quit_now = true;
    }
    fn new_document(&mut self) {
        self.sync_active_document();
        let doc = DocumentTab::blank();
        self.documents.push(doc);
        self.active_document = self.documents.len() - 1;
        self.mol = Molecule::default();
        self.graphics.clear();
        self.dirty = false;
        self.undo.clear();
        self.redo.clear();
        self.clear_selection();
        self.zoom = self.preferred_zoom;
        self.pan = Vec2::ZERO;
        self.status = "New document tab".into();
    }
    fn open_snapshot_tab(&mut self, snapshot: Snapshot, path: Option<String>) {
        self.sync_active_document();
        if self.documents.len() == 1
            && !self.dirty
            && self.mol.atoms.is_empty()
            && self.graphics.is_empty()
            && self.documents[0].path.is_none()
        {
            self.documents.pop();
        }
        let title = path
            .as_deref()
            .and_then(|p| Path::new(p).file_name())
            .and_then(|n| n.to_str())
            .unwrap_or("Imported structure")
            .to_string();
        self.documents.push(DocumentTab {
            title,
            path,
            snapshot: snapshot.clone(),
            dirty: false,
            undo: Vec::new(),
            redo: Vec::new(),
        });
        self.active_document = self.documents.len() - 1;
        self.mol = snapshot.mol;
        self.graphics = snapshot.graphics;
        self.dirty = false;
        self.undo.clear();
        self.redo.clear();
        self.clear_selection();
        self.zoom = self.preferred_zoom;
        self.pan = Vec2::ZERO;
    }
    fn autosave(&mut self) {
        self.sync_active_document();
        if !self.documents.iter().any(|d| d.dirty) {
            let _ = fs::remove_file(autosave_path());
            return;
        }
        let path = autosave_path();
        let contents = native::write_workspace(&self.documents, self.active_document);
        let result = path
            .parent()
            .map(fs::create_dir_all)
            .unwrap_or(Ok(()))
            .and_then(|_| fs::write(&path, contents));
        if result.is_err() {
            self.status = "Autosave failed; save the document manually".into();
        }
    }
    fn restore_autosave(&mut self) {
        let Ok(text) = fs::read_to_string(autosave_path()) else {
            return;
        };
        if self.stage_recovery(&text).is_err() {
            return;
        }
        self.status = "An autosaved workspace is available to recover".into();
    }
    fn stage_recovery(&mut self, text: &str) -> Result<(), String> {
        let (docs, active) = native::read_workspace(text)?;
        if docs.is_empty() {
            return Err("Autosave workspace is empty".into());
        }
        self.pending_recovery = Some((docs, active));
        Ok(())
    }
    fn apply_recovery(&mut self, recover: bool) {
        let Some((mut docs, active)) = self.pending_recovery.take() else {
            return;
        };
        if recover {
            for doc in &mut docs {
                doc.dirty = true;
            }
            self.documents = docs;
            self.active_document = active.min(self.documents.len().saturating_sub(1));
            let doc = &self.documents[self.active_document];
            self.mol = doc.snapshot.mol.clone();
            self.graphics = doc.snapshot.graphics.clone();
            self.undo.clear();
            self.redo.clear();
            self.dirty = true;
            self.status = format!(
                "Recovered autosaved workspace · {} tab(s)",
                self.documents.len()
            );
        } else {
            let _ = fs::remove_file(autosave_path());
            self.status = "Autosaved recovery discarded; original files were not changed".into();
        }
    }
    fn record_recent(&mut self, path: &str) {
        self.recent_files.retain(|p| p != path);
        self.recent_files.insert(0, path.to_string());
        self.recent_files.truncate(10);
        let _ = save_recent_files(&self.recent_files);
    }
    fn hit_graphic(&self, p: Vec2) -> Option<usize> {
        self.graphics
            .iter()
            .enumerate()
            .rev()
            .find(|(_, g)| match g {
                Graphic::Arrow { start, end, .. } => distance_segment(p, *start, *end) < 10.,
                Graphic::ElectronArrow {
                    start,
                    control,
                    end,
                    ..
                } => (0..16).any(|k| {
                    let t = k as f32 / 16.;
                    let a =
                        *start * (1. - t).powi(2) + *control * (2. * t * (1. - t)) + *end * t * t;
                    let u = (k + 1) as f32 / 16.;
                    let b =
                        *start * (1. - u).powi(2) + *control * (2. * u * (1. - u)) + *end * u * u;
                    distance_segment(p, a, b) < 10.
                }),
                Graphic::Plus { pos } => p.distance(*pos) < 18.,
                Graphic::Shape {
                    kind: ShapeKind::Line,
                    start,
                    end,
                } => distance_segment(p, *start, *end) < 9.,
                Graphic::Shape {
                    kind: ShapeKind::Rectangle,
                    start,
                    end,
                } => {
                    let (x0, x1) = (start.x.min(end.x), start.x.max(end.x));
                    let (y0, y1) = (start.y.min(end.y), start.y.max(end.y));
                    p.x >= x0 - 7.
                        && p.x <= x1 + 7.
                        && p.y >= y0 - 7.
                        && p.y <= y1 + 7.
                        && ((p.x - x0).abs() < 9.
                            || (p.x - x1).abs() < 9.
                            || (p.y - y0).abs() < 9.
                            || (p.y - y1).abs() < 9.)
                }
                Graphic::Shape {
                    kind: ShapeKind::Ellipse,
                    start,
                    end,
                } => {
                    let rx = (start.x - end.x).abs() / 2.;
                    let ry = (start.y - end.y).abs() / 2.;
                    let c = (*start + *end) / 2.;
                    rx > 1.
                        && ry > 1.
                        && (((p.x - c.x) / rx).powi(2) + ((p.y - c.y) / ry).powi(2) - 1.).abs()
                            < 0.22
                }
                Graphic::Text {
                    pos, text, size, ..
                } => {
                    p.x >= pos.x
                        && p.x <= pos.x + text.len() as f32 * size * 0.6
                        && p.y >= pos.y - size
                        && p.y <= pos.y
                }
            })
            .map(|(i, _)| i)
    }
    fn add_atom(&mut self, p: Vec2, element: &str) -> usize {
        self.mol.atoms.push(Atom::new(p, element));
        self.mol.atoms.len() - 1
    }
    fn draw_chain(&mut self, start: Vec2, finish: Vec2) {
        let a_existing = self.hit_atom(start);
        let origin = a_existing.map(|i| self.mol.atoms[i].pos).unwrap_or(start);
        let delta = finish - origin;
        let length = delta.length();
        let direction = if length > 1. {
            delta / length
        } else {
            vec2(1., 0.)
        };
        let count = if length < 18. {
            3
        } else {
            ((length / BOND).round() as usize).clamp(2, 10)
        };
        let end = if length < 18. {
            origin + direction * (BOND * count as f32)
        } else {
            finish
        };
        let normal = vec2(-direction.y, direction.x);
        let mut points = Vec::with_capacity(count + 1);
        for k in 0..=count {
            let t = k as f32 / count as f32;
            let zig = if k == 0 || k == count {
                0.
            } else if k % 2 == 1 {
                -1.
            } else {
                1.
            };
            points.push(origin + (end - origin) * t + normal * zig * (BOND * 0.22));
        }
        points[0] = origin;
        points[count] = end;
        let first = a_existing.unwrap_or_else(|| self.add_atom(origin, "C"));
        let last_existing = self.hit_atom(end).filter(|i| *i != first);
        let mut prev = first;
        for point in points.iter().take(count).skip(1) {
            let next = self.add_atom(*point, "C");
            self.add_bond(prev, next, 1);
            prev = next;
        }
        let last = last_existing.unwrap_or_else(|| self.add_atom(end, "C"));
        self.add_bond(prev, last, 1);
        self.status = format!("Drew a {}-bond carbon chain", count);
    }
    fn add_bond(&mut self, a: usize, b: usize, order: u8) {
        if a == b {
            return;
        }
        if let Some(x) = self
            .mol
            .bonds
            .iter_mut()
            .find(|x| (x.a == a && x.b == b) || (x.a == b && x.b == a))
        {
            set_bond_order_value(x, order);
        } else {
            self.mol.bonds.push(Bond {
                a,
                b,
                order,
                style: BondStyle::Normal,
                stereo: None,
                query: None,
            });
        }
    }
    fn add_styled_bond(&mut self, a: usize, b: usize, style: BondStyle) {
        self.add_bond(a, b, 1);
        if let Some(bond) = self
            .mol
            .bonds
            .iter_mut()
            .find(|x| (x.a == a && x.b == b) || (x.a == b && x.b == a))
        {
            bond.style = style;
            bond.order = if style == BondStyle::Aromatic { 4 } else { 1 };
        }
    }
    fn set_bond_order(&mut self, index: usize, order: u8) {
        let order = order.clamp(1, 4);
        if let Some(bond) = self.mol.bonds.get(index) {
            let has_order_query = bond.query.as_deref().is_some_and(|query| {
                query.split(';').any(|term| {
                    matches!(
                        term,
                        "any" | "single-or-double" | "single-or-aromatic" | "double-or-aromatic"
                    )
                })
            });
            let expected_style = if order == 4 {
                BondStyle::Aromatic
            } else {
                BondStyle::Normal
            };
            if bond.order == order && bond.style == expected_style && !has_order_query {
                return;
            }
        } else {
            return;
        }
        self.checkpoint();
        set_bond_order_value(&mut self.mol.bonds[index], order);
        self.refresh_aromatic_atoms();
    }
    fn set_bond_query(&mut self, index: usize, query: &str) {
        let Some(bond) = self.mol.bonds.get(index) else {
            return;
        };
        let mut terms: Vec<&str> = bond
            .query
            .as_deref()
            .unwrap_or("")
            .split(';')
            .filter(|term| matches!(*term, "topology-ring" | "topology-chain"))
            .collect();
        terms.extend(query.split(';').filter(|term| !term.is_empty()));
        terms.sort_unstable();
        terms.dedup();
        let combined = terms.join(";");
        if bond.order == 1
            && bond.style == BondStyle::Any
            && bond.stereo.is_none()
            && bond.query.as_deref() == Some(combined.as_str())
        {
            return;
        }
        self.checkpoint();
        let bond = &mut self.mol.bonds[index];
        bond.order = 1;
        bond.style = BondStyle::Any;
        bond.stereo = None;
        bond.query = Some(combined);
        self.refresh_aromatic_atoms();
    }
    fn set_bond_style(&mut self, index: usize, style: BondStyle) {
        let Some(bond) = self.mol.bonds.get(index) else {
            return;
        };
        let order = if style == BondStyle::Aromatic { 4 } else { 1 };
        if bond.order == order && bond.style == style {
            return;
        }
        self.checkpoint();
        let bond = &mut self.mol.bonds[index];
        set_bond_order_value(bond, order);
        bond.style = style;
        self.refresh_aromatic_atoms();
    }
    fn refresh_aromatic_atoms(&mut self) {
        let aromatic_bonds: BTreeSet<usize> = self
            .mol
            .bonds
            .iter()
            .filter(|bond| bond.order == 4)
            .flat_map(|bond| [bond.a, bond.b])
            .collect();
        for (atom_index, atom) in self.mol.atoms.iter_mut().enumerate() {
            atom.aromatic = aromatic_bonds.contains(&atom_index);
        }
    }
    fn apply_bond_tool_to_hit_bond(&mut self, point: Vec2) -> bool {
        // Atom hits keep the usual click-and-drag workflow for drawing from an
        // atom. Clicking the middle of an existing bond edits that bond in one
        // step, matching the common chemistry-editor interaction.
        if self.hit_atom(point).is_some() {
            return false;
        }
        let Some(index) = self.hit_bond(point) else {
            return false;
        };
        match self.tool {
            Tool::Bond(order) => self.set_bond_order(index, order),
            Tool::StyledBond(BondStyle::Any) if self.bond_query.is_some() => {
                let query = self.bond_query.clone().unwrap();
                self.set_bond_query(index, &query);
            }
            Tool::StyledBond(style) => self.set_bond_style(index, style),
            _ => return false,
        }
        self.selected = None;
        self.selected_atoms.clear();
        self.selected_bond = Some(index);
        self.selected_bonds = vec![index];
        self.selected_graphic = None;
        self.selected_graphics.clear();
        let label = match self.mol.bonds[index].style {
            BondStyle::Aromatic => "aromatic",
            BondStyle::Wedge => "wedge",
            BondStyle::Hash => "hashed wedge",
            BondStyle::Wavy => "wavy",
            BondStyle::Coordinate => "coordinate",
            BondStyle::Any => "any-order",
            BondStyle::Normal => match self.mol.bonds[index].order {
                1 => "single",
                2 => "double",
                3 => "triple",
                _ => "bond",
            },
        };
        self.status = if self.tool == Tool::Bond(2) {
            localized(self.language, "Changed existing bond to a double bond").into()
        } else {
            format!("Changed existing bond to {label}")
        };
        true
    }
    fn cycle_bond_style(&mut self, index: usize) {
        let Some(bond) = self.mol.bonds.get(index) else {
            return;
        };
        let style = next_bond_style(bond.style);
        self.checkpoint();
        let bond = &mut self.mol.bonds[index];
        if bond.order == 4 {
            bond.order = 1;
        }
        if matches!(style, BondStyle::Wedge | BondStyle::Hash | BondStyle::Wavy) && bond.order != 1
        {
            bond.order = 1;
            bond.stereo = None;
        }
        bond.style = style;
        let aromatic_bonds: BTreeSet<usize> = self
            .mol
            .bonds
            .iter()
            .filter(|bond| bond.order == 4)
            .flat_map(|bond| [bond.a, bond.b])
            .collect();
        for (atom_index, atom) in self.mol.atoms.iter_mut().enumerate() {
            atom.aromatic = aromatic_bonds.contains(&atom_index);
        }
    }
    fn cycle_bond_topology(&mut self, index: usize) {
        let Some(bond) = self.mol.bonds.get(index) else {
            return;
        };
        let current = bond
            .query
            .as_deref()
            .unwrap_or("")
            .split(';')
            .find(|term| matches!(*term, "topology-ring" | "topology-chain"));
        let next = match current {
            Some("topology-ring") => Some("topology-chain"),
            Some("topology-chain") => None,
            _ => Some("topology-ring"),
        };
        let mut terms: Vec<&str> = bond
            .query
            .as_deref()
            .unwrap_or("")
            .split(';')
            .filter(|term| !term.is_empty() && !matches!(*term, "topology-ring" | "topology-chain"))
            .collect();
        if let Some(next) = next {
            terms.push(next);
        }
        let query = (!terms.is_empty()).then(|| terms.join(";"));
        if bond.query == query {
            return;
        }
        self.checkpoint();
        self.mol.bonds[index].query = query;
    }
    fn clear_atom_query(&mut self, index: usize) {
        if self
            .mol
            .atoms
            .get(index)
            .is_some_and(|atom| atom.query.is_some())
        {
            self.checkpoint();
            self.mol.atoms[index].query = None;
        }
    }
    fn cycle_atom_radical(&mut self, index: usize) {
        if let Some(atom) = self.mol.atoms.get(index) {
            let next = (atom.radical + 1) % 4;
            if next != atom.radical {
                self.checkpoint();
                self.mol.atoms[index].radical = next;
            }
        }
    }
    fn cycle_implicit_h_override(&mut self, index: usize) {
        if let Some(atom) = self.mol.atoms.get(index) {
            let next = match atom.implicit_h_override {
                None => Some(0),
                Some(count @ 0..=4) => Some(count + 1),
                Some(_) => None,
            };
            if atom.implicit_h_override != next {
                self.checkpoint();
                self.mol.atoms[index].implicit_h_override = next;
            }
        }
    }
    fn cycle_double_bond_stereo(&mut self, index: usize) {
        let Some(bond) = self.mol.bonds.get(index) else {
            return;
        };
        if bond.order != 2 {
            return;
        }
        self.checkpoint();
        let bond = &mut self.mol.bonds[index];
        bond.stereo = next_double_bond_stereo(bond.stereo.as_deref());
    }
    fn selected_fragment_for_template(&self) -> Result<Molecule, String> {
        let mut selected = self.selected_atoms.iter().copied().collect::<BTreeSet<_>>();
        if selected.is_empty() {
            if let Some(index) = self.selected {
                selected.insert(index);
            }
        }
        selected.retain(|index| *index < self.mol.atoms.len());
        if selected.is_empty() {
            return Err("Select a connected atom fragment before saving a template".into());
        }
        let center = selected
            .iter()
            .map(|index| self.mol.atoms[*index].pos)
            .sum::<Vec2>()
            / selected.len() as f32;
        let mut remap = vec![None; self.mol.atoms.len()];
        let mut fragment = Molecule::default();
        for source in &selected {
            let mut atom = self.mol.atoms[*source].clone();
            atom.pos -= center;
            atom.atom_map = None;
            atom.reaction_role = None;
            atom.reaction_step = 1;
            atom.coefficient = 1;
            remap[*source] = Some(fragment.atoms.len());
            fragment.atoms.push(atom);
        }
        for bond in &self.mol.bonds {
            if let (Some(a), Some(b)) = (remap[bond.a], remap[bond.b]) {
                let mut bond = bond.clone();
                bond.a = a;
                bond.b = b;
                fragment.bonds.push(bond);
            }
        }
        if molecule_components(&fragment).len() != 1 {
            return Err("A saved template must be one connected fragment".into());
        }
        Ok(fragment)
    }
    fn save_selected_user_template(&mut self, name: &str) -> Result<(), String> {
        let fragment = self.selected_fragment_for_template()?;
        let template = save_user_template_to(&user_templates_path(), name, &fragment)?;
        self.user_templates.push(template);
        self.user_templates
            .sort_by_key(|template| template.name.to_lowercase());
        self.status = format!("Saved custom template: {}", name.trim());
        Ok(())
    }
    fn delete_user_template(&mut self, index: usize) -> Result<(), String> {
        let template = self
            .user_templates
            .get(index)
            .ok_or("Select a custom template to delete")?;
        fs::remove_file(&template.path).map_err(|error| error.to_string())?;
        let name = template.name.clone();
        self.user_templates.remove(index);
        self.template_page = self
            .template_page
            .min(self.user_templates.len().saturating_sub(1) / 6);
        self.status = format!("Deleted custom template: {name}");
        Ok(())
    }
    fn rename_user_template(&mut self, index: usize, name: &str) -> Result<(), String> {
        self.rename_user_template_in(index, name, &user_templates_path())
    }
    fn rename_user_template_in(
        &mut self,
        index: usize,
        name: &str,
        directory: &Path,
    ) -> Result<(), String> {
        let original = self
            .user_templates
            .get(index)
            .cloned()
            .ok_or("Select a custom template to rename")?;
        if original.name == name.trim() {
            return Ok(());
        }
        let renamed = save_user_template_to(directory, name, &original.mol)?;
        if let Err(error) = fs::remove_file(&original.path) {
            let cleanup = fs::remove_file(&renamed.path);
            return Err(match cleanup {
                Ok(()) => error.to_string(),
                Err(cleanup_error) => format!(
                    "Could not remove the old template ({error}); the new copy is at {} ({cleanup_error})",
                    renamed.path.display()
                ),
            });
        }
        self.user_templates.remove(index);
        self.user_templates.push(renamed);
        self.user_templates
            .sort_by_key(|template| template.name.to_lowercase());
        self.status = format!("Renamed custom template to {}", name.trim());
        Ok(())
    }
    fn insert_custom_template(&mut self, mut template: Molecule, center: Vec2) {
        if template.atoms.is_empty() {
            self.status = "The selected custom template is empty".into();
            return;
        }
        self.checkpoint();
        let origin = template
            .atoms
            .iter()
            .fold(vec2(0., 0.), |sum, atom| sum + atom.pos)
            / template.atoms.len() as f32;
        let rotation = self.template_rotation;
        for atom in &mut template.atoms {
            atom.pos += center - origin;
            atom.pos = rotate_template_point(atom.pos, center, rotation);
            atom.atom_map = None;
            atom.reaction_role = None;
            atom.reaction_step = 1;
            atom.coefficient = 1;
        }
        let offset = self.mol.atoms.len();
        self.mol.atoms.extend(template.atoms);
        self.mol
            .bonds
            .extend(template.bonds.into_iter().map(|mut bond| {
                bond.a += offset;
                bond.b += offset;
                bond
            }));
        self.selected_atoms = (offset..self.mol.atoms.len()).collect();
        self.selected = self.selected_atoms.first().copied();
        self.template_rotation = 0;
        self.status = if rotation == 0 {
            format!("Placed custom template: {}", template.name)
        } else {
            format!(
                "Placed custom template rotated {rotation}°: {}",
                template.name
            )
        };
    }
    fn activate_builtin_template(&mut self, kind: &'static str) {
        let mut scratch = Self::new();
        scratch.insert_template(kind, vec2(0., 0.));
        self.template_preview = Some(scratch.mol);
        self.pending_custom_template = None;
        self.template_attachment_target = None;
        self.template_rotation = 0;
        self.quick.pending.clear();
        self.tool = Tool::Template(kind);
        self.templates_open = false;
        self.status = format!("Click canvas to place {kind}; press R to rotate");
    }
    fn activate_custom_template(&mut self, template: Molecule) {
        let name = template.name.clone();
        let mut preview = template.clone();
        if !preview.atoms.is_empty() {
            let origin = preview
                .atoms
                .iter()
                .fold(vec2(0., 0.), |sum, atom| sum + atom.pos)
                / preview.atoms.len() as f32;
            for atom in &mut preview.atoms {
                atom.pos -= origin;
            }
        }
        self.template_preview = Some(preview);
        self.pending_custom_template = Some(template);
        self.template_attachment_target = None;
        self.template_rotation = 0;
        self.quick.pending.clear();
        self.tool = Tool::Template("User");
        self.templates_open = false;
        self.status = format!("Click canvas to place {name}; press R to rotate");
    }
    fn template_anchor(&self) -> Option<(usize, Vec2)> {
        let preview = self.template_preview.as_ref()?;
        let (index, atom) =
            preview.atoms.iter().enumerate().min_by(|(_, a), (_, b)| {
                a.pos.length_squared().total_cmp(&b.pos.length_squared())
            })?;
        Some((index, atom.pos))
    }
    fn template_placement_center(&self, cursor: Vec2) -> Vec2 {
        let Some(target) = self.template_attachment_target else {
            return cursor;
        };
        let Some(target_atom) = self.mol.atoms.get(target) else {
            return cursor;
        };
        let Some((_, anchor)) = self.template_anchor() else {
            return cursor;
        };
        let direction = cursor - target_atom.pos;
        let direction = if direction.length() < 1. {
            vec2(1., 0.)
        } else {
            direction.normalize()
        };
        target_atom.pos + direction * BOND
            - rotate_template_point(anchor, vec2(0., 0.), self.template_rotation)
    }
    fn connect_template_atom(&mut self, target: Option<usize>, inserted: usize, anchor: usize) {
        let Some(target) = target.filter(|index| *index < inserted) else {
            return;
        };
        let anchor = inserted + anchor;
        if anchor >= self.mol.atoms.len() {
            return;
        }
        self.mol.bonds.push(Bond {
            a: target,
            b: anchor,
            order: 1,
            style: BondStyle::Normal,
            stereo: None,
            query: None,
        });
        self.selected_atoms = vec![target, anchor];
        self.selected = Some(anchor);
        self.status = "Template attached to the selected atom with a single bond".into();
    }
    fn toggle_template_attachment(&mut self) {
        if !matches!(self.tool, Tool::Template(_)) {
            return;
        }
        if self.template_attachment_target.is_some() {
            self.template_attachment_target = None;
            self.status = "Template will be placed as a separate fragment".into();
            return;
        }
        if self.selected_atoms.len() != 1 {
            self.status = "Select exactly one atom, then press A to attach the template".into();
            return;
        }
        let target = self.selected_atoms[0];
        if target >= self.mol.atoms.len() {
            self.status = "Selected attachment atom is no longer available".into();
            return;
        }
        self.template_attachment_target = Some(target);
        self.status =
            "Template will attach to the selected atom; move pointer to set direction".into();
    }
    fn cancel_template_placement(&mut self) {
        if matches!(self.tool, Tool::Template(_)) {
            self.tool = Tool::Select;
            self.template_preview = None;
            self.pending_custom_template = None;
            self.template_rotation = 0;
            self.template_attachment_target = None;
            self.quick.pending.clear();
            self.status = "Template placement cancelled".into();
        }
    }
    fn remove_atom(&mut self, i: usize) {
        self.mol.atoms.remove(i);
        self.mol.bonds.retain(|b| b.a != i && b.b != i);
        for b in &mut self.mol.bonds {
            if b.a > i {
                b.a -= 1;
            }
            if b.b > i {
                b.b -= 1;
            }
        }
    }
    fn ring(&mut self, center: Vec2, n: usize) {
        self.ring_with_snap(center, n, true);
    }
    fn ring_with_snap(&mut self, center: Vec2, n: usize, snap_to_existing: bool) {
        let radius = BOND / (2. * (std::f32::consts::PI / n as f32).sin());
        let mut ids = Vec::new();
        for k in 0..n {
            let t = -std::f32::consts::FRAC_PI_2 + 2. * std::f32::consts::PI * k as f32 / n as f32;
            let p = center + vec2(t.cos(), t.sin()) * radius;
            let existing = snap_to_existing.then(|| self.hit_atom(p)).flatten();
            let id = existing.unwrap_or_else(|| self.add_atom(p, "C"));
            ids.push(id);
        }
        for k in 0..n {
            self.add_bond(ids[k], ids[(k + 1) % n], 1);
        }
    }
    fn fused_ring_template(
        &mut self,
        center: Vec2,
        first_size: usize,
        second_size: usize,
        aromatic: bool,
        heteroatom: Option<(&str, usize, u8)>,
    ) {
        let apothem = |size: usize| BOND / (2. * (std::f32::consts::PI / size as f32).tan());
        let first_apothem = apothem(first_size);
        let second_apothem = apothem(second_size);
        let mut first = Vec::with_capacity(first_size);
        for k in 0..first_size {
            let angle = std::f32::consts::FRAC_PI_2 - std::f32::consts::PI / first_size as f32
                + 2. * std::f32::consts::PI * k as f32 / first_size as f32;
            let pos = center + vec2(angle.cos() * BOND, angle.sin() * BOND - first_apothem);
            first.push(self.add_atom(pos, "C"));
        }
        let mut second = Vec::with_capacity(second_size);
        for k in 0..second_size {
            if k == 0 {
                second.push(first[1]);
                continue;
            }
            if k == 1 {
                second.push(first[0]);
                continue;
            }
            let angle = 3. * std::f32::consts::FRAC_PI_2
                - std::f32::consts::PI / second_size as f32
                + 2. * std::f32::consts::PI * k as f32 / second_size as f32;
            let pos = center + vec2(angle.cos() * BOND, angle.sin() * BOND + second_apothem);
            let element = heteroatom
                .filter(|(_, ring_position, _)| *ring_position == k)
                .map_or("C", |(element, _, _)| element);
            let id = self.add_atom(pos, element);
            if let Some((_, ring_position, hydrogen)) = heteroatom {
                if ring_position == k {
                    self.mol.atoms[id].explicit_h = hydrogen;
                }
            }
            second.push(id);
        }
        for ring in [&first[..], &second[..]] {
            for k in 0..ring.len() {
                let a = ring[k];
                let b = ring[(k + 1) % ring.len()];
                if self
                    .mol
                    .bonds
                    .iter()
                    .any(|bond| (bond.a == a && bond.b == b) || (bond.a == b && bond.b == a))
                {
                    continue;
                }
                self.add_bond(a, b, if aromatic { 4 } else { 1 });
                if aromatic {
                    if let Some(bond) = self.mol.bonds.last_mut() {
                        bond.style = BondStyle::Aromatic;
                    }
                }
            }
        }
        if aromatic {
            for atom in first.iter().chain(second.iter()) {
                self.mol.atoms[*atom].aromatic = true;
            }
        }
    }
    fn norbornane_template(&mut self, center: Vec2) {
        let points = [
            center + vec2(-BOND, 0.),
            center + vec2(BOND, 0.),
            center + vec2(-BOND * 0.5, -BOND * 0.8660254),
            center + vec2(BOND * 0.5, -BOND * 0.8660254),
            center + vec2(-BOND * 0.5, BOND * 0.8660254),
            center + vec2(BOND * 0.5, BOND * 0.8660254),
            center,
        ];
        let ids: Vec<usize> = points.iter().map(|pos| self.add_atom(*pos, "C")).collect();
        for (a, b) in [
            (0, 2),
            (2, 3),
            (3, 1),
            (0, 4),
            (4, 5),
            (5, 1),
            (0, 6),
            (6, 1),
        ] {
            self.add_bond(ids[a], ids[b], 1);
        }
    }
    fn aromatic_tricyclic_template(&mut self, center: Vec2, angular: bool) {
        let root_three = 1.7320508 * BOND;
        let translation = if angular {
            vec2(-0.75 * BOND, 0.4330127 * BOND)
        } else {
            vec2(0., 0.)
        };
        let ring_centers = if angular {
            [
                center + translation + vec2(0., -root_three),
                center + translation,
                center + translation + vec2(1.5 * BOND, 0.5 * root_three),
            ]
        } else {
            [
                center + vec2(0., -root_three),
                center,
                center + vec2(0., root_three),
            ]
        };
        let mut coordinate_atoms = BTreeMap::<(i32, i32), usize>::new();
        let mut rings = Vec::with_capacity(3);
        for ring_center in ring_centers {
            let mut ring = Vec::with_capacity(6);
            for k in 0..6 {
                let angle = k as f32 * std::f32::consts::PI / 3.;
                let pos = ring_center + vec2(angle.cos(), angle.sin()) * BOND;
                let key = (
                    (pos.x * 1000.).round() as i32,
                    (pos.y * 1000.).round() as i32,
                );
                let atom = *coordinate_atoms
                    .entry(key)
                    .or_insert_with(|| self.add_atom(pos, "C"));
                ring.push(atom);
            }
            rings.push(ring);
        }
        let mut added = BTreeSet::new();
        for ring in rings {
            for k in 0..6 {
                let (a, b) = (ring[k], ring[(k + 1) % 6]);
                let edge = (a.min(b), a.max(b));
                if added.insert(edge) {
                    self.add_bond(a, b, 4);
                    if let Some(bond) = self.mol.bonds.last_mut() {
                        bond.style = BondStyle::Aromatic;
                    }
                }
            }
        }
        for atom in coordinate_atoms.values() {
            self.mol.atoms[*atom].aromatic = true;
        }
    }
    fn add_explicit_hydrogens(&mut self) {
        let analysis = chemistry::analyze(&self.mol);
        let original_len = self.mol.atoms.len();
        let mut planned = Vec::new();
        for i in 0..original_len {
            let atom = &self.mol.atoms[i];
            if atom.element == "H" {
                continue;
            }
            let Some(implicit) = analysis.implicit_h.get(i).copied().flatten() else {
                continue;
            };
            let total = implicit.saturating_add(atom.explicit_h).min(8);
            if total == 0 {
                continue;
            }
            let origin = atom.pos;
            let neighbors: Vec<Vec2> = self
                .mol
                .bonds
                .iter()
                .filter_map(|b| {
                    let other = if b.a == i {
                        Some(b.b)
                    } else if b.b == i {
                        Some(b.a)
                    } else {
                        None
                    }?;
                    let d = self.mol.atoms[other].pos - origin;
                    (d.length() > 0.1).then(|| d.normalize())
                })
                .collect();
            let sum = neighbors.iter().copied().sum::<Vec2>();
            let base = if sum.length() > 0.1 {
                -sum.normalize()
            } else {
                vec2(1., 0.)
            };
            for k in 0..total {
                let angle = base.y.atan2(base.x) + std::f32::consts::TAU * k as f32 / total as f32;
                planned.push((i, origin + vec2(angle.cos(), angle.sin()) * (BOND * 0.72)));
            }
        }
        if planned.is_empty() {
            self.status = "No implicit hydrogens to expand".into();
            return;
        }
        self.checkpoint();
        for (parent, pos) in planned {
            let h = self.add_atom(pos, "H");
            self.add_bond(parent, h, 1);
        }
        for atom in &mut self.mol.atoms[..original_len] {
            if atom.element != "H" {
                atom.explicit_h = 0;
                atom.implicit_h_override = None;
            }
        }
        self.status = "Added explicit hydrogen atoms".into();
    }
    fn remove_explicit_hydrogens(&mut self) {
        let ids: Vec<usize> = self
            .mol
            .atoms
            .iter()
            .enumerate()
            .filter_map(|(i, a)| {
                if a.element != "H" {
                    return None;
                }
                self.mol
                    .bonds
                    .iter()
                    .any(|b| {
                        (b.a == i && self.mol.atoms[b.b].element != "H")
                            || (b.b == i && self.mol.atoms[b.a].element != "H")
                    })
                    .then_some(i)
            })
            .collect();
        if ids.is_empty() {
            self.status = "No attached explicit hydrogens to remove".into();
            return;
        }
        self.checkpoint();
        for i in ids.into_iter().rev() {
            self.remove_atom(i);
        }
        self.status = "Removed explicit hydrogen atoms".into();
    }
    fn transform_document(&mut self, angle: f32, flip_x: bool, flip_y: bool) {
        if self.mol.atoms.is_empty() {
            self.status = "Draw a structure first".into();
            return;
        }
        self.checkpoint();
        let center =
            self.mol.atoms.iter().map(|a| a.pos).sum::<Vec2>() / self.mol.atoms.len() as f32;
        let transform = |p: Vec2| {
            let mut q = p - center;
            if flip_x {
                q.x = -q.x;
            }
            if flip_y {
                q.y = -q.y;
            }
            let (s, c) = angle.sin_cos();
            center + vec2(q.x * c - q.y * s, q.x * s + q.y * c)
        };
        for atom in &mut self.mol.atoms {
            atom.pos = transform(atom.pos);
        }
        for graphic in &mut self.graphics {
            match graphic {
                Graphic::Arrow { start, end, .. } => {
                    *start = transform(*start);
                    *end = transform(*end);
                }
                Graphic::Text { pos, .. } | Graphic::Plus { pos } => *pos = transform(*pos),
                Graphic::Shape { start, end, .. } => {
                    *start = transform(*start);
                    *end = transform(*end);
                }
                Graphic::ElectronArrow {
                    start,
                    control,
                    end,
                    ..
                } => {
                    *start = transform(*start);
                    *control = transform(*control);
                    *end = transform(*end);
                }
            }
        }
        if flip_x ^ flip_y {
            for atom in &mut self.mol.atoms {
                if let Some(s) = &mut atom.stereo {
                    swap_stereo(s);
                }
            }
            for bond in &mut self.mol.bonds {
                if let Some(s) = &mut bond.stereo {
                    swap_stereo(s);
                }
            }
        }
        self.status = if flip_x || flip_y {
            "Mirrored document; stereo labels inverted".into()
        } else {
            "Rotated document in 2D".into()
        };
    }
    fn clean_2d(&mut self) {
        let n = self.mol.atoms.len();
        if n < 2 {
            self.status = "Draw at least two atoms before cleaning".into();
            return;
        }
        if n > 400 {
            self.status = "Clean 2D is limited to 400 atoms for responsiveness".into();
            return;
        }
        self.checkpoint();
        let mut seen = vec![false; n];
        for root in 0..n {
            if seen[root] {
                continue;
            }
            let mut stack = vec![root];
            let mut group = Vec::new();
            seen[root] = true;
            while let Some(i) = stack.pop() {
                group.push(i);
                for b in &self.mol.bonds {
                    let j = if b.a == i {
                        Some(b.b)
                    } else if b.b == i {
                        Some(b.a)
                    } else {
                        None
                    };
                    if let Some(j) = j {
                        if !seen[j] {
                            seen[j] = true;
                            stack.push(j);
                        }
                    }
                }
            }
            if group.len() < 2 {
                continue;
            }
            let map: std::collections::BTreeMap<usize, usize> =
                group.iter().enumerate().map(|(j, &i)| (i, j)).collect();
            let edges: Vec<(usize, usize)> = self
                .mol
                .bonds
                .iter()
                .filter_map(|b| Some((*map.get(&b.a)?, *map.get(&b.b)?)))
                .collect();
            let mut p: Vec<Vec2> = group.iter().map(|&i| self.mol.atoms[i].pos).collect();
            let k = BOND;
            for step in 0..100 {
                let mut d = vec![vec2(0., 0.); p.len()];
                for i in 0..p.len() {
                    for j in i + 1..p.len() {
                        let mut v = p[i] - p[j];
                        let mut len = v.length();
                        if len < 0.01 {
                            v = vec2((i + 1) as f32, (j + 1) as f32);
                            len = v.length();
                        }
                        let f = k * k / len;
                        let force = v / len * f;
                        d[i] += force;
                        d[j] -= force;
                    }
                }
                for &(i, j) in &edges {
                    let v = p[i] - p[j];
                    let len = v.length().max(0.01);
                    let force = v / len * (len * len / k);
                    d[i] -= force;
                    d[j] += force;
                }
                let temp = 8. * (1. - step as f32 / 100.);
                for i in 0..p.len() {
                    let len = d[i].length();
                    if len > 0.001 {
                        p[i] += d[i] / len * len.min(temp);
                    }
                }
            }
            let old =
                group.iter().map(|&i| self.mol.atoms[i].pos).sum::<Vec2>() / group.len() as f32;
            let new = p.iter().copied().sum::<Vec2>() / p.len() as f32;
            for (j, &i) in group.iter().enumerate() {
                self.mol.atoms[i].pos = p[j] + old - new;
            }
        }
        self.status = "Cleaned 2D coordinates using a spring layout".into();
    }
    fn add_aromatic_template_ring(
        &mut self,
        center: Vec2,
        size: usize,
        heteroatoms: &[(usize, &str, bool)],
    ) -> Vec<usize> {
        let radius = BOND / (2.0 * (std::f32::consts::PI / size as f32).sin());
        let mut atoms = Vec::with_capacity(size);
        for index in 0..size {
            let angle =
                -std::f32::consts::FRAC_PI_2 + std::f32::consts::TAU * index as f32 / size as f32;
            let element = heteroatoms
                .iter()
                .find(|(position, _, _)| *position == index)
                .map(|(_, element, _)| *element)
                .unwrap_or("C");
            let atom = self.add_atom(center + vec2(angle.cos(), angle.sin()) * radius, element);
            self.mol.atoms[atom].aromatic = true;
            if heteroatoms
                .iter()
                .any(|(position, _, explicit_h)| *position == index && *explicit_h)
            {
                self.mol.atoms[atom].explicit_h = 1;
            }
            atoms.push(atom);
        }
        for index in 0..size {
            self.add_bond(atoms[index], atoms[(index + 1) % size], 4);
            if let Some(bond) = self.mol.bonds.iter_mut().find(|bond| {
                (bond.a == atoms[index] && bond.b == atoms[(index + 1) % size])
                    || (bond.b == atoms[index] && bond.a == atoms[(index + 1) % size])
            }) {
                bond.style = BondStyle::Aromatic;
            }
        }
        atoms
    }
    fn insert_amino_acid(&mut self, kind: &str, center: Vec2) {
        let at = |x: f32, y: f32| center + vec2(BOND * x, BOND * y);
        let n = self.add_atom(at(-1.0, 0.0), "N");
        let ca = self.add_atom(center, "C");
        let carbonyl = self.add_atom(at(1.0, 0.0), "C");
        let carbonyl_o = self.add_atom(at(1.5, -0.866), "O");
        let hydroxyl_o = self.add_atom(at(1.5, 0.866), "O");
        self.add_bond(n, ca, 1);
        self.add_bond(ca, carbonyl, 1);
        self.add_bond(carbonyl, carbonyl_o, 2);
        self.add_bond(carbonyl, hydroxyl_o, 1);

        let add = |editor: &mut Self, pos: Vec2, element: &str| editor.add_atom(pos, element);
        let attach = |editor: &mut Self, a: usize, b: usize, order| editor.add_bond(a, b, order);
        match kind {
            "Glycine" => {}
            "Alanine" => {
                let beta = add(self, at(0.0, 1.0), "C");
                attach(self, ca, beta, 1);
            }
            "Valine" => {
                let beta = add(self, at(0.0, 1.0), "C");
                let methyl_a = add(self, at(-0.866, 1.5), "C");
                let methyl_b = add(self, at(0.866, 1.5), "C");
                attach(self, ca, beta, 1);
                attach(self, beta, methyl_a, 1);
                attach(self, beta, methyl_b, 1);
            }
            "Leucine" => {
                let beta = add(self, at(0.0, 1.0), "C");
                let gamma = add(self, at(0.866, 1.5), "C");
                let delta_a = add(self, at(1.732, 1.0), "C");
                let delta_b = add(self, at(0.866, 2.5), "C");
                attach(self, ca, beta, 1);
                attach(self, beta, gamma, 1);
                attach(self, gamma, delta_a, 1);
                attach(self, gamma, delta_b, 1);
            }
            "Isoleucine" => {
                let beta = add(self, at(0.0, 1.0), "C");
                let gamma = add(self, at(0.866, 1.5), "C");
                let methyl_a = add(self, at(-0.866, 1.5), "C");
                let methyl_b = add(self, at(1.732, 1.0), "C");
                attach(self, ca, beta, 1);
                attach(self, beta, gamma, 1);
                attach(self, beta, methyl_a, 1);
                attach(self, gamma, methyl_b, 1);
            }
            "Serine" | "Cysteine" => {
                let beta = add(self, at(0.0, 1.0), "C");
                let end = add(
                    self,
                    at(0.866, 1.5),
                    if kind == "Serine" { "O" } else { "S" },
                );
                attach(self, ca, beta, 1);
                attach(self, beta, end, 1);
            }
            "Threonine" => {
                let beta = add(self, at(0.0, 1.0), "C");
                let hydroxyl = add(self, at(-0.866, 1.5), "O");
                let methyl = add(self, at(0.866, 1.5), "C");
                attach(self, ca, beta, 1);
                attach(self, beta, hydroxyl, 1);
                attach(self, beta, methyl, 1);
            }
            "Methionine" => {
                let beta = add(self, at(0.0, 1.0), "C");
                let gamma = add(self, at(0.866, 1.5), "C");
                let sulfur = add(self, at(1.732, 1.0), "S");
                let methyl = add(self, at(2.598, 1.5), "C");
                attach(self, ca, beta, 1);
                attach(self, beta, gamma, 1);
                attach(self, gamma, sulfur, 1);
                attach(self, sulfur, methyl, 1);
            }
            "Aspartic acid" | "Asparagine" | "Glutamic acid" | "Glutamine" => {
                let beta = add(self, at(0.0, 1.0), "C");
                attach(self, ca, beta, 1);
                let mut previous = beta;
                if matches!(kind, "Glutamic acid" | "Glutamine") {
                    let gamma = add(self, at(0.866, 1.5), "C");
                    attach(self, beta, gamma, 1);
                    previous = gamma;
                }
                let side_carbonyl = add(self, at(1.732, 1.0), "C");
                let side_oxo = add(self, at(2.232, 0.134), "O");
                let side_hetero = add(
                    self,
                    at(2.232, 1.866),
                    if matches!(kind, "Asparagine" | "Glutamine") {
                        "N"
                    } else {
                        "O"
                    },
                );
                attach(self, previous, side_carbonyl, 1);
                attach(self, side_carbonyl, side_oxo, 2);
                attach(self, side_carbonyl, side_hetero, 1);
            }
            "Lysine" => {
                let b = add(self, at(0.0, 1.0), "C");
                let g = add(self, at(0.866, 1.5), "C");
                let d = add(self, at(1.732, 1.0), "C");
                let e = add(self, at(2.598, 1.5), "C");
                let z = add(self, at(3.464, 1.0), "N");
                attach(self, ca, b, 1);
                attach(self, b, g, 1);
                attach(self, g, d, 1);
                attach(self, d, e, 1);
                attach(self, e, z, 1);
            }
            "Arginine" => {
                let b = add(self, at(0.0, 1.0), "C");
                let g = add(self, at(0.866, 1.5), "C");
                let d = add(self, at(1.732, 1.0), "C");
                let e = add(self, at(2.598, 1.5), "N");
                let z = add(self, at(3.464, 1.0), "C");
                let n_double = add(self, at(4.33, 0.5), "N");
                let n_single_a = add(self, at(4.33, 1.5), "N");
                attach(self, ca, b, 1);
                attach(self, b, g, 1);
                attach(self, g, d, 1);
                attach(self, d, e, 1);
                attach(self, e, z, 1);
                attach(self, z, n_double, 2);
                attach(self, z, n_single_a, 1);
            }
            "Phenylalanine" | "Tyrosine" => {
                let beta = add(self, at(0.0, 1.0), "C");
                attach(self, ca, beta, 1);
                let ring = self.add_aromatic_template_ring(at(0.0, 3.0), 6, &[]);
                attach(self, beta, ring[0], 1);
                if kind == "Tyrosine" {
                    let oxygen = add(self, at(1.0, 4.0), "O");
                    attach(self, ring[3], oxygen, 1);
                }
            }
            "Histidine" => {
                let beta = add(self, at(0.0, 1.0), "C");
                attach(self, ca, beta, 1);
                let ring = self.add_aromatic_template_ring(
                    at(0.0, 3.0),
                    5,
                    &[(1, "N", false), (3, "N", true)],
                );
                attach(self, beta, ring[0], 1);
            }
            "Tryptophan" => {
                let beta = add(self, at(0.0, 1.0), "C");
                attach(self, ca, beta, 1);
                let benzene = self.add_aromatic_template_ring(at(0.0, 3.0), 6, &[]);
                attach(self, beta, benzene[0], 1);
                let b = self.mol.atoms[benzene[3]].pos;
                let positions = [
                    b + vec2(-0.743, -0.669) * BOND,
                    b + vec2(-0.336, -1.583) * BOND,
                    b + vec2(0.659, -1.478) * BOND,
                ];
                let c = add(self, positions[0], "C");
                self.mol.atoms[c].aromatic = true;
                let nh = add(self, positions[1], "N");
                self.mol.atoms[nh].aromatic = true;
                self.mol.atoms[nh].explicit_h = 1;
                let d = add(self, positions[2], "C");
                self.mol.atoms[d].aromatic = true;
                for (u, v) in [(benzene[2], c), (c, nh), (nh, d), (d, benzene[3])] {
                    attach(self, u, v, 4);
                    if let Some(bond) =
                        self.mol.bonds.iter_mut().find(|bond| {
                            (bond.a == u && bond.b == v) || (bond.a == v && bond.b == u)
                        })
                    {
                        bond.style = BondStyle::Aromatic;
                    }
                }
            }
            "Proline" => {
                let beta = add(self, at(0.0, 1.0), "C");
                let gamma = add(self, at(-0.866, 1.5), "C");
                let delta = add(self, at(-1.589, 0.809), "C");
                attach(self, ca, beta, 1);
                attach(self, beta, gamma, 1);
                attach(self, gamma, delta, 1);
                attach(self, delta, n, 1);
            }
            _ => {}
        }
    }
    fn suggest_endpoint(&self, i: usize) -> Vec2 {
        let origin = self.mol.atoms[i].pos;
        let neighbors: Vec<Vec2> = self
            .mol
            .bonds
            .iter()
            .filter_map(|b| {
                if b.a == i {
                    Some(self.mol.atoms[b.b].pos - origin)
                } else if b.b == i {
                    Some(self.mol.atoms[b.a].pos - origin)
                } else {
                    None
                }
            })
            .collect();
        if neighbors.is_empty() {
            return origin + vec2(BOND, 0.);
        }
        let sum = neighbors.iter().copied().sum::<Vec2>();
        let dir = if sum.length() < 0.1 {
            vec2(0., -1.)
        } else {
            -sum.normalize()
        };
        let angle = if neighbors.len() == 1 {
            std::f32::consts::FRAC_PI_3
        } else {
            0.
        };
        let out = vec2(
            dir.x * angle.cos() - dir.y * angle.sin(),
            dir.x * angle.sin() + dir.y * angle.cos(),
        );
        origin + out * BOND
    }
    fn insert_template(&mut self, kind: &str, center: Vec2) {
        if kind == "User" {
            if let Some(template) = self.pending_custom_template.clone() {
                let inserted = self.mol.atoms.len();
                let attachment_target = self.template_attachment_target;
                let anchor = self.template_anchor().map(|(index, _)| index).unwrap_or(0);
                let center = self.template_placement_center(center);
                self.insert_custom_template(template, center);
                if attachment_target.is_some() {
                    self.connect_template_atom(attachment_target, inserted, anchor);
                    self.template_attachment_target = None;
                }
            } else {
                self.status = "Choose a custom template in the template library first".into();
            }
            return;
        }
        let center = self.template_placement_center(center);
        self.checkpoint();
        let first_atom = self.mol.atoms.len();
        let rotation = self.template_rotation;
        let attachment_target = self.template_attachment_target;
        let anchor = self.template_anchor().map(|(index, _)| index).unwrap_or(0);
        match kind {
            "Benzene" | "Cyclohexane" | "Pyridine" | "Pyrimidine" | "Imidazole" | "Pyrrole"
            | "Furan" | "Thiophene" => {
                let first_b = self.mol.bonds.len();
                let first_a = self.mol.atoms.len();
                let size = if matches!(kind, "Imidazole" | "Pyrrole" | "Furan" | "Thiophene") {
                    5
                } else {
                    6
                };
                self.ring_with_snap(center, size, rotation == 0 && attachment_target.is_none());
                if kind != "Cyclohexane" {
                    for b in self.mol.bonds.iter_mut().skip(first_b) {
                        b.style = BondStyle::Aromatic;
                        b.order = 4;
                    }
                    for a in self.mol.atoms.iter_mut().skip(first_a) {
                        a.aromatic = true;
                    }
                }
                if kind == "Pyridine" && self.mol.atoms.len() > first_a {
                    self.mol.atoms[first_a].element = "N".into();
                }
                if kind == "Pyrimidine" && self.mol.atoms.len() >= first_a + 4 {
                    self.mol.atoms[first_a].element = "N".into();
                    self.mol.atoms[first_a + 3].element = "N".into();
                }
                if kind == "Imidazole" && self.mol.atoms.len() >= first_a + 3 {
                    self.mol.atoms[first_a].element = "N".into();
                    self.mol.atoms[first_a + 2].element = "N".into();
                    // Imidazole has one pyrrolic [NH] and one pyridine-like N.
                    // The charge-neutral pyrrolic nitrogen contributes one H.
                    self.mol.atoms[first_a].explicit_h = 1;
                }
                if kind == "Pyrrole" && self.mol.atoms.len() >= first_a + 5 {
                    self.mol.atoms[first_a].element = "N".into();
                    self.mol.atoms[first_a].explicit_h = 1;
                }
                if kind == "Furan" && self.mol.atoms.len() >= first_a + 5 {
                    self.mol.atoms[first_a].element = "O".into();
                }
                if kind == "Thiophene" && self.mol.atoms.len() >= first_a + 5 {
                    self.mol.atoms[first_a].element = "S".into();
                }
            }
            "Cyclopropane" | "Cyclobutane" | "Cyclopentane" | "Cycloheptane" | "Cyclooctane" => {
                let n = match kind {
                    "Cyclopropane" => 3,
                    "Cyclobutane" => 4,
                    "Cyclopentane" => 5,
                    "Cycloheptane" => 7,
                    _ => 8,
                };
                self.ring_with_snap(center, n, rotation == 0 && attachment_target.is_none());
            }
            "Naphthalene" | "Decalin" => {
                self.fused_ring_template(center, 6, 6, kind == "Naphthalene", None)
            }
            "Indole" => self.fused_ring_template(center, 6, 5, true, Some(("N", 2, 1))),
            "Quinoline" => self.fused_ring_template(center, 6, 6, true, Some(("N", 2, 0))),
            "Anthracene" => self.aromatic_tricyclic_template(center, false),
            "Phenanthrene" => self.aromatic_tricyclic_template(center, true),
            "Norbornane" => self.norbornane_template(center),
            "Ethanol" => {
                let a = self.add_atom(center - vec2(BOND, 0.), "C");
                let b = self.add_atom(center, "C");
                let c = self.add_atom(center + vec2(BOND, 0.), "O");
                self.add_bond(a, b, 1);
                self.add_bond(b, c, 1);
            }
            "Carboxyl" => {
                let a = self.add_atom(center, "C");
                let b = self.add_atom(center + vec2(BOND * 0.8, -BOND * 0.6), "O");
                let c = self.add_atom(center + vec2(BOND * 0.8, BOND * 0.6), "O");
                self.add_bond(a, b, 2);
                self.add_bond(a, c, 1);
            }
            "Hydroxyl" | "Amino" | "Carbonyl" | "Aldehyde" | "Nitrile" => {
                let (terminal, order) = match kind {
                    "Hydroxyl" => ("O", 1),
                    "Amino" => ("N", 1),
                    "Carbonyl" => ("O", 2),
                    "Aldehyde" => ("O", 2),
                    _ => ("N", 3),
                };
                let center_atom = self.add_atom(center, "C");
                let end = self.add_atom(center + vec2(BOND, 0.), terminal);
                self.add_bond(center_atom, end, order);
                if kind == "Aldehyde" {
                    let stem = self.add_atom(center - vec2(BOND, 0.), "C");
                    self.add_bond(stem, center_atom, 1);
                }
            }
            "Ester" | "Amide" => {
                let carbon = self.add_atom(center, "C");
                let carbonyl = self.add_atom(center + vec2(BOND * 0.75, -BOND * 0.7), "O");
                let hetero = self.add_atom(
                    center + vec2(BOND * 0.8, BOND * 0.65),
                    if kind == "Ester" { "O" } else { "N" },
                );
                let methyl = self.add_atom(center - vec2(BOND, 0.), "C");
                self.add_bond(carbon, carbonyl, 2);
                self.add_bond(carbon, hetero, 1);
                self.add_bond(carbon, methyl, 1);
            }
            "Glycine" | "Alanine" | "Valine" | "Leucine" | "Isoleucine" | "Serine"
            | "Threonine" | "Cysteine" | "Methionine" | "Aspartic acid" | "Asparagine"
            | "Glutamic acid" | "Glutamine" | "Lysine" | "Arginine" | "Histidine"
            | "Phenylalanine" | "Tyrosine" | "Tryptophan" | "Proline" => {
                self.insert_amino_acid(kind, center);
            }
            _ => {}
        }
        for atom in self.mol.atoms.iter_mut().skip(first_atom) {
            atom.pos = rotate_template_point(atom.pos, center, rotation);
        }
        if attachment_target.is_some() {
            self.connect_template_atom(attachment_target, first_atom, anchor);
            self.template_attachment_target = None;
        }
        self.template_rotation = 0;
    }
    fn save(&mut self, path: &str) {
        let losses = self.export_losses(path);
        if !losses.is_empty() {
            self.pending_export = Some(path.to_string());
            self.status = format!("Export will lose: {}", losses.join(", "));
            return;
        }
        self.save_confirmed(path);
    }
    fn save_confirmed(&mut self, path: &str) {
        let lower = path.to_lowercase();
        let export_only = [
            ".svg",
            ".png",
            ".pdf",
            ".smi",
            ".smiles",
            ".smarts",
            ".cml",
            ".pdb",
            ".inchi",
            ".inchikey",
        ]
        .iter()
        .any(|ext| lower.ends_with(ext));
        let result = if lower.ends_with(".svg") {
            self.write_svg(path, screen_width(), screen_height())
        } else if lower.ends_with(".png") {
            render_png(&self.mol, &self.graphics, 1920, 1080)
                .and_then(|bytes| fs::write(path, bytes).map_err(|e| e.to_string()))
                .map_err(std::io::Error::other)
        } else if lower.ends_with(".pdf") {
            fs::write(path, render_pdf(&self.mol, &self.graphics, 842., 595.))
        } else if lower.ends_with(".mols") {
            fs::write(
                path,
                native::write(&Snapshot {
                    mol: self.mol.clone(),
                    graphics: self.graphics.clone(),
                }),
            )
        } else if lower.ends_with(".mol") {
            fs::write(path, self.molfile())
        } else if lower.ends_with(".sdf") {
            let mut text = molfile_for(&self.mol);
            for (key, value) in &self.mol.properties {
                text += &format!(">  <{}>\n{}\n\n", key.replace(['<', '>', '\n'], ""), value);
            }
            text += "$$$$\n";
            fs::write(path, text)
        } else if lower.ends_with(".rxn") {
            match rxn_for(&self.mol) {
                Ok(text) => fs::write(path, text),
                Err(error) => {
                    self.status = error;
                    return;
                }
            }
        } else if lower.ends_with(".smi") || lower.ends_with(".smiles") {
            let block = self.molfile();
            match smiles::export(&self.mol).or_else(|local_error| {
                rdkit::to_smiles_with_stereo(&block, &self.mol).map_err(|rdkit_error| {
                    format!("Local SMILES exporter: {local_error}; optional RDKit fallback: {rdkit_error}")
                })
            }) {
                Ok(smiles) => fs::write(path, format!("{}\n", smiles.trim())),
                Err(error) => {
                    self.status = format!("SMILES export unavailable: {error}");
                    return;
                }
            }
        } else if lower.ends_with(".smarts") {
            match smarts::export(&self.mol) {
                Ok(pattern) => fs::write(path, format!("{pattern}\n")),
                Err(error) => {
                    self.status = format!("SMARTS export unavailable: {error}");
                    return;
                }
            }
        } else if lower.ends_with(".cml") {
            fs::write(path, cml::write(&self.mol))
        } else if [".pdb", ".inchi", ".inchikey"]
            .iter()
            .any(|ext| lower.ends_with(ext))
        {
            let operation = if lower.ends_with(".pdb") {
                "to-pdb"
            } else if lower.ends_with(".inchikey") {
                "to-inchikey"
            } else {
                "to-inchi"
            };
            match rdkit::export_external(operation, &self.molfile()) {
                Ok(data) => fs::write(path, data),
                Err(error) => {
                    self.status = format!(
                        "{} export requires optional RDKit: {error}",
                        operation.trim_start_matches("to-")
                    );
                    return;
                }
            }
        } else {
            self.status =
                "Unsupported extension. Use .mols, .mol, .sdf, .rxn, .smi, .smarts, .cml, .pdb, .inchi, .inchikey, .svg, .png or .pdf"
                    .into();
            return;
        };
        self.status = match result {
            Ok(_) => {
                if !export_only {
                    let title = Path::new(path)
                        .file_name()
                        .and_then(|s| s.to_str())
                        .unwrap_or(path)
                        .to_string();
                    if let Some(doc) = self.documents.get_mut(self.active_document) {
                        doc.title = title;
                        doc.path = Some(path.into());
                        doc.dirty = false;
                    }
                    self.dirty = false;
                }
                self.record_recent(path);
                format!("{} {path}", if export_only { "Exported" } else { "Saved" })
            }
            Err(e) => format!("Save failed: {e}"),
        };
        if self.status.starts_with("Saved ") {
            if let Some(index) = self.close_after_save.take() {
                self.close_document(index);
            }
        }
        if !self.documents.iter().any(|doc| doc.dirty) {
            let _ = fs::remove_file(autosave_path());
        }
    }
    fn save_all(&mut self) {
        self.sync_active_document();
        for doc in &self.documents {
            let Some(path) = doc.path.as_deref() else {
                continue;
            };
            if !path.to_lowercase().ends_with(".mol") {
                continue;
            }
            let losses = Self::export_losses_for(&doc.snapshot.mol, &doc.snapshot.graphics, path);
            if !losses.is_empty() {
                self.pending_export = Some(path.to_string());
                self.status = format!(
                    "Save all needs export confirmation for {}: {}",
                    path,
                    losses.join(", ")
                );
                return;
            }
        }
        let mut saved = 0usize;
        for (i, doc) in self.documents.iter_mut().enumerate() {
            let Some(path) = doc.path.clone() else {
                continue;
            };
            let lower = path.to_lowercase();
            let result = if lower.ends_with(".mols") {
                fs::write(&path, native::write(&doc.snapshot))
            } else if lower.ends_with(".mol") {
                fs::write(&path, molfile_for(&doc.snapshot.mol))
            } else {
                continue;
            };
            match result {
                Ok(()) => {
                    doc.dirty = false;
                    saved += 1;
                }
                Err(err) => {
                    self.status = format!("Save all stopped at {path}: {err}");
                    self.dirty = self
                        .documents
                        .get(self.active_document)
                        .is_some_and(|d| d.dirty);
                    return;
                }
            }
            if i == self.active_document {
                self.dirty = false;
            }
        }
        self.status = if saved == 0 {
            "No named documents to save".into()
        } else {
            format!("Saved {saved} document(s)")
        };
    }
    fn load(&mut self, path: &str) {
        match fs::read_to_string(path) {
            Ok(s) => {
                let lower = path.to_lowercase();
                if lower.ends_with(".smarts") {
                    match smarts::parse(&s) {
                        Ok(mol) => {
                            self.open_snapshot_tab(
                                Snapshot {
                                    mol,
                                    graphics: vec![],
                                },
                                Some(path.into()),
                            );
                            self.record_recent(path);
                            self.status = format!("Opened SMARTS pattern {path}");
                        }
                        Err(error) => self.status = format!("SMARTS import failed: {error}"),
                    }
                    return;
                }
                if lower.ends_with(".smi") || lower.ends_with(".smiles") {
                    let smiles = s
                        .lines()
                        .find(|line| !line.trim().is_empty() && !line.trim_start().starts_with('#'))
                        .unwrap_or("")
                        .split_whitespace()
                        .next()
                        .unwrap_or("");
                    let imported = smiles::parse(smiles).or_else(|local_error| {
                        rdkit::from_smiles_with_stereo(smiles)
                            .and_then(|(molblock, stereo)| {
                                let mut mol = parse_mol(&molblock)?;
                                for (index, value) in stereo {
                                    let bond = mol.bonds.get_mut(index).ok_or(
                                        "RDKit stereo annotation references a missing bond",
                                    )?;
                                    if bond.order != 2 {
                                        return Err("RDKit assigned alkene stereo to a non-double bond".into());
                                    }
                                    bond.stereo = Some(value);
                                }
                                Ok(mol)
                            })
                            .map_err(|rdkit_error| {
                                format!("Local SMILES parser: {local_error}; optional RDKit fallback: {rdkit_error}")
                            })
                    });
                    match imported {
                        Ok(mol) => {
                            self.open_snapshot_tab(
                                Snapshot {
                                    mol,
                                    graphics: vec![],
                                },
                                Some(path.into()),
                            );
                            self.record_recent(path);
                            self.status = format!("Opened SMILES {path}");
                        }
                        Err(error) => self.status = format!("SMILES import unavailable: {error}"),
                    }
                    return;
                }
                let import_op = if lower.ends_with(".mol2") {
                    Some("from-mol2")
                } else if lower.ends_with(".pdb") {
                    Some("from-pdb")
                } else if lower.ends_with(".inchi") {
                    Some("from-inchi")
                } else {
                    None
                };
                if let Some(operation) = import_op {
                    if self.refuse_overlapping_worker() {
                        return;
                    }
                    let captured = self.mol.clone();
                    let filename = path.to_string();
                    let (tx, rx) = mpsc::channel();
                    std::thread::spawn(move || {
                        let result = rdkit::import_external(operation, &s)
                            .map(|report| format!("IMPORT\t{filename}\n{report}"));
                        let _ = tx.send(result);
                    });
                    self.calc_rx = Some((captured, rx));
                    self.status = format!(
                        "Importing {} using optional RDKit…",
                        operation.trim_start_matches("from-")
                    );
                    return;
                }
                if lower.ends_with(".sdf") {
                    match parse_sdf(&s) {
                        Ok(records) if !records.is_empty() => {
                            for mol in records {
                                self.open_snapshot_tab(
                                    Snapshot {
                                        mol,
                                        graphics: vec![],
                                    },
                                    Some(path.into()),
                                );
                            }
                            self.record_recent(path);
                            self.status = format!(
                                "Opened {} SDF record(s)",
                                s.matches("$$$$").count().max(1)
                            );
                        }
                        Ok(_) => self.status = "SDF contains no molecule records".into(),
                        Err(error) => self.status = format!("SDF import failed: {error}"),
                    }
                    return;
                }
                if lower.ends_with(".rxn") {
                    match parse_rxn(
                        &s,
                        vec2(screen_width() / 2., (TOP + screen_height() - BOTTOM) / 2.),
                    ) {
                        Ok(mol) => {
                            self.open_snapshot_tab(
                                Snapshot {
                                    mol,
                                    graphics: vec![],
                                },
                                Some(path.into()),
                            );
                            self.record_recent(path);
                            self.status = format!("Opened reaction {path}");
                        }
                        Err(error) => self.status = format!("RXN import failed: {error}"),
                    }
                    return;
                }
                if lower.ends_with(".cml") {
                    match cml::parse(
                        &s,
                        vec2(screen_width() / 2., (TOP + screen_height() - BOTTOM) / 2.),
                    ) {
                        Ok(mol) => {
                            self.open_snapshot_tab(
                                Snapshot {
                                    mol,
                                    graphics: vec![],
                                },
                                Some(path.into()),
                            );
                            self.record_recent(path);
                            self.status = format!("Opened CML molecule {path}");
                        }
                        Err(error) => self.status = format!("CML import failed: {error}"),
                    }
                    return;
                }
                let parsed = if lower.ends_with(".mols") {
                    native::read(&s)
                } else {
                    parse_mol(&s).map(|mol| Snapshot {
                        mol,
                        graphics: vec![],
                    })
                };
                match parsed {
                    Ok(snapshot) => {
                        self.open_snapshot_tab(snapshot, Some(path.into()));
                        self.record_recent(path);
                        self.status = format!("Opened {path}");
                    }
                    Err(e) => self.status = e,
                }
            }
            Err(e) => self.status = format!("Open failed: {e}"),
        }
    }
    fn export_losses(&self, path: &str) -> Vec<&'static str> {
        Self::export_losses_for(&self.mol, &self.graphics, path)
    }
    fn export_losses_for(mol: &Molecule, graphics: &[Graphic], path: &str) -> Vec<&'static str> {
        if path.to_lowercase().ends_with(".mols") {
            return vec![];
        }
        if path.to_lowercase().ends_with(".svg") {
            return vec!["editable chemical graph and objects"];
        }
        if path.to_lowercase().ends_with(".png") || path.to_lowercase().ends_with(".pdf") {
            return vec!["editable chemical graph and molecule metadata"];
        }
        if path.to_lowercase().ends_with(".rxn") {
            return vec!["reaction arrows, conditions, atom maps, coefficients, multi-step data and graphical annotations"];
        }
        if path.to_lowercase().ends_with(".smarts") {
            return vec!["document metadata, reaction annotation and drawing objects"];
        }
        if path.to_lowercase().ends_with(".cml") {
            return cml::export_losses(mol, !graphics.is_empty());
        }
        if [".smi", ".smiles"]
            .iter()
            .any(|extension| path.to_lowercase().ends_with(extension))
        {
            let mut losses = Vec::new();
            if !graphics.is_empty() {
                losses.push("graphics, reaction arrows, conditions and text");
            }
            if !mol.name.is_empty() || !mol.properties.is_empty() {
                losses.push("molecule name and data fields");
            }
            if mol.atoms.iter().any(|atom| {
                atom.pos.length_squared() > 0.01 || atom.z.is_some() || atom.coord3d.is_some()
            }) {
                losses.push("2D/3D coordinates and drawing layout");
            }
            if mol.atoms.iter().any(|atom| {
                atom.alias.is_some()
                    || atom.r_group.is_some()
                    || atom.query.is_some()
                    || atom.radical != 0
                    || atom.valence_override.is_some()
                    || atom.stereo.is_some()
            }) {
                losses.push("query, alias, R-group, radical, valence or stereo annotations");
            }
            if mol.bonds.iter().any(|bond| {
                bond.query.is_some()
                    || bond.stereo.is_some()
                    || !matches!(bond.style, BondStyle::Normal | BondStyle::Aromatic)
            }) {
                losses.push("bond query, stereo or display annotations");
            }
            if mol.atoms.iter().any(|atom| {
                atom.reaction_role.is_some() || atom.reaction_step != 1 || atom.coefficient != 1
            }) {
                losses.push("reaction roles, steps and stoichiometric coefficients");
            }
            return losses;
        }
        let is_mol_exchange = [".mol", ".sdf"]
            .iter()
            .any(|extension| path.to_lowercase().ends_with(extension));
        let mol_v3000 = is_mol_exchange && mol_uses_v3000(mol);
        let mol_v2000 = is_mol_exchange && !mol_v3000;
        let aromatic_atoms: BTreeSet<usize> = mol
            .bonds
            .iter()
            .filter(|bond| bond.order == 4)
            .flat_map(|bond| [bond.a, bond.b])
            .collect();
        let mut v = Vec::new();
        if !graphics.is_empty() {
            v.push("graphics, reaction arrows and text");
        }
        if mol.atoms.iter().enumerate().any(|(index, atom)| {
            atom.alias.is_some()
                || atom.query.as_deref().is_some_and(|query| {
                    !mol_v3000 || v3000_query_attributes(Some(query)).1 || v3000_atom_type(atom).1
                })
                || atom.r_group.is_some()
                || atom.z.is_some()
                || (atom.valence_override.is_some() && !is_mol_exchange)
                || (atom.valence_override.is_some_and(|v| v > 14) && mol_v2000)
                || (atom.implicit_h_override.is_some()
                    && (!is_mol_exchange
                        || (mol_v2000 && atom.implicit_h_override.is_some_and(|h| h > 4))
                        || (!mol_v2000 && atom.implicit_h_override.is_some_and(|h| h > 5))))
                || atom.explicit_h > 0
                || (atom.aromatic && !aromatic_atoms.contains(&index))
                || (atom.stereo.is_some()
                    && (!is_mol_exchange
                        || !matches!(atom.stereo.as_deref(), Some("@" | "@@" | "either"))))
        }) {
            v.push("advanced atom properties");
        }
        if mol.bonds.iter().any(|bond| {
            bond.stereo.is_some()
                || (bond.query.is_some() && (!mol_v3000 || v3000_bond_attributes(bond).2))
                || (is_mol_exchange && bond.style == BondStyle::Coordinate)
                || (is_mol_exchange && bond.style == BondStyle::Any && !mol_v3000)
                || (!is_mol_exchange && bond.style != BondStyle::Normal)
        }) {
            v.push("bond styles or query data");
        }
        if mol.atoms.iter().any(|a| a.atom_map.is_some()) {
            v.push("atom maps");
        }
        if mol
            .atoms
            .iter()
            .any(|a| a.reaction_role.is_some() || a.reaction_step != 1 || a.coefficient != 1)
        {
            v.push("reaction roles, steps or stoichiometric coefficients");
        }
        v
    }
    fn request_save(&mut self, path: String) {
        self.save(&path);
    }
    fn commit_text(&mut self) {
        if let Some(target) = self.text_target.take() {
            match target {
                TextTarget::TemplateName => {
                    match self.save_selected_user_template(&self.text_buffer.clone()) {
                        Ok(()) => self.template_tab_custom = true,
                        Err(error) => self.status = error,
                    }
                    self.text_buffer.clear();
                    return;
                }
                TextTarget::RenameUserTemplate(index) => {
                    if let Err(error) = self.rename_user_template(index, &self.text_buffer.clone())
                    {
                        self.status = error;
                    }
                    self.text_buffer.clear();
                    return;
                }
                TextTarget::MoleculeProperty => {
                    self.apply_molecule_property_edit(&self.text_buffer.clone());
                    self.text_buffer.clear();
                    return;
                }
                TextTarget::ReactionCoefficient => {
                    let value = self.text_buffer.trim().parse::<u8>();
                    match value {
                        Ok(value) if (1..=99).contains(&value) => {
                            self.set_reaction_component(None, None, Some(value));
                        }
                        _ => {
                            self.status =
                                "Enter a whole-number reaction coefficient from 1 to 99".into();
                        }
                    }
                    self.text_buffer.clear();
                    return;
                }
                TextTarget::BufferPh
                | TextTarget::PkaFromBuffer
                | TextTarget::WeakAcidPh
                | TextTarget::WeakBasePh
                | TextTarget::StrongMixturePh
                | TextTarget::PolyproticSpeciation
                | TextTarget::BjerrumDiagram
                | TextTarget::PolyproticTitrationPh
                | TextTarget::PolyproticTitrationCurve => {
                    self.calculate_acid_base(target, &self.text_buffer.clone());
                    self.text_buffer.clear();
                    return;
                }
                _ => {}
            }
            self.checkpoint();
            match target {
                TextTarget::New(pos) => self.graphics.push(Graphic::Text {
                    pos,
                    text: self.text_buffer.clone(),
                    size: 20.,
                    bold: false,
                    italic: false,
                }),
                TextTarget::ArrowAbove(i) => {
                    if let Some(Graphic::Arrow { above, .. }) = self.graphics.get_mut(i) {
                        *above = self.text_buffer.clone();
                    }
                }
                TextTarget::ArrowBelow(i) => {
                    if let Some(Graphic::Arrow { below, .. }) = self.graphics.get_mut(i) {
                        *below = self.text_buffer.clone();
                    }
                }
                TextTarget::TemplateName => unreachable!(),
                TextTarget::RenameUserTemplate(_) => unreachable!(),
                TextTarget::MoleculeProperty => unreachable!(),
                TextTarget::ReactionCoefficient => unreachable!(),
                TextTarget::BufferPh
                | TextTarget::PkaFromBuffer
                | TextTarget::WeakAcidPh
                | TextTarget::WeakBasePh
                | TextTarget::StrongMixturePh
                | TextTarget::PolyproticSpeciation
                | TextTarget::BjerrumDiagram
                | TextTarget::PolyproticTitrationPh
                | TextTarget::PolyproticTitrationCurve => unreachable!(),
            }
            self.text_buffer.clear();
        }
    }
    fn start_acid_base_input(&mut self, target: TextTarget) {
        self.acid_base_open = false;
        self.text_buffer.clear();
        self.text_target = Some(target);
    }
    fn calculate_acid_base(&mut self, target: TextTarget, input: &str) {
        self.acid_base_chart = None;
        let parse_fields = |count: usize| -> Result<Vec<f64>, String> {
            let fields: Vec<&str> = input.split(';').map(str::trim).collect();
            if fields.len() != count {
                return Err(format!("Enter exactly {count} semicolon-separated values"));
            }
            fields
                .iter()
                .map(|field| {
                    field
                        .replace(',', ".")
                        .parse::<f64>()
                        .map_err(|_| format!("Invalid numeric value: {field}"))
                })
                .collect()
        };
        let parse_pkas = |text: &str| -> Result<Vec<f64>, String> {
            text.split('|')
                .map(|value| {
                    value
                        .trim()
                        .replace(',', ".")
                        .parse::<f64>()
                        .map_err(|_| format!("Invalid pKa value: {value}"))
                })
                .collect()
        };
        let result = (|| -> Result<String, String> {
            match target {
                TextTarget::BufferPh => {
                    let v = parse_fields(3)?;
                    let ph = acid_base::buffer_ph(v[0], v[1], v[2])?;
                    Ok(format!(
                        "Buffer pH = {ph:.4}\n\npKa = {:.4}\n[HA] = {:.6} mol/L\n[A⁻] = {:.6} mol/L\nMethod: Henderson–Hasselbalch, pH = pKa + log₁₀([A⁻]/[HA]). Aqueous ideal-solution calculation.",
                        v[0], v[1], v[2]
                    ))
                }
                TextTarget::PkaFromBuffer => {
                    let v = parse_fields(3)?;
                    let pka = acid_base::pka_from_buffer(v[0], v[1], v[2])?;
                    Ok(format!(
                        "pKa = {pka:.4}\n\nMeasured pH = {:.4}\n[HA] = {:.6} mol/L\n[A⁻] = {:.6} mol/L\nMethod: rearranged Henderson–Hasselbalch. This derives a pKa from supplied measurement data; it is not a structure-based prediction.",
                        v[0], v[1], v[2]
                    ))
                }
                TextTarget::WeakAcidPh => {
                    let v = parse_fields(2)?;
                    let ph = acid_base::weak_acid_ph(v[0], v[1])?;
                    Ok(format!(
                        "Weak-acid pH = {ph:.4}\n\npKa = {:.4}\nInitial concentration = {:.6} mol/L\nMethod: ideal aqueous monoprotic equilibrium solved by charge balance, including water autoionization; Kw = 1.0×10⁻¹⁴ at 25 °C.",
                        v[0], v[1]
                    ))
                }
                TextTarget::WeakBasePh => {
                    let v = parse_fields(2)?;
                    let ph = acid_base::weak_base_ph(v[0], v[1])?;
                    Ok(format!(
                        "Weak-base pH = {ph:.4}\n\npKa of BH⁺ = {:.4}\nInitial base concentration = {:.6} mol/L\nMethod: ideal aqueous conjugate-acid equilibrium solved by charge balance, including water autoionization; Kw = 1.0×10⁻¹⁴ at 25 °C.",
                        v[0], v[1]
                    ))
                }
                TextTarget::StrongMixturePh => {
                    let v = parse_fields(2)?;
                    let ph = acid_base::strong_mixture_ph(v[0], v[1])?;
                    Ok(format!(
                        "Strong acid/base mixture pH = {ph:.4}\n\nFinal strong-acid concentration = {:.6} mol/L\nFinal strong-base concentration = {:.6} mol/L\nMethod: neutralization and water ion balance, Kw = 1.0×10⁻¹⁴ at 25 °C. Enter concentrations after any volume change.",
                        v[0], v[1]
                    ))
                }
                TextTarget::PolyproticSpeciation => {
                    let (pka_text, ph_text) = input
                        .split_once(';')
                        .ok_or("Separate the pKa list and pH with a semicolon")?;
                    let pkas = parse_pkas(pka_text)?;
                    let ph: f64 = ph_text
                        .trim()
                        .replace(',', ".")
                        .parse()
                        .map_err(|_| "Invalid pH value")?;
                    let fractions = acid_base::polyprotic_fractions(ph, &pkas)?;
                    let n = pkas.len();
                    let rows = fractions
                        .iter()
                        .enumerate()
                        .map(|(index, fraction)| {
                            let lost = index;
                            let formula = if lost == 0 {
                                if n == 1 {
                                    "HA".to_string()
                                } else {
                                    format!("H{n}A")
                                }
                            } else if lost == n {
                                format!("A{}−", charge_suffix(lost))
                            } else {
                                let hydrogen_count = n - lost;
                                let hydrogen = if hydrogen_count == 1 {
                                    "H".to_string()
                                } else {
                                    format!("H{hydrogen_count}")
                                };
                                format!("{hydrogen}A{}−", charge_suffix(lost))
                            };
                            format!("{formula}: {:.3}%", fraction * 100.0)
                        })
                        .collect::<Vec<_>>()
                        .join("\n");
                    Ok(format!(
                        "Polyprotic acid speciation at pH {ph:.3}\n\npKa values: {}\n\n{rows}\n\nMethod: successive ideal aqueous dissociation equilibria; fractions sum to 100%.",
                        pkas.iter().map(|pka| format!("{pka:.3}")).collect::<Vec<_>>().join(", ")
                    ))
                }
                TextTarget::BjerrumDiagram => {
                    let fields: Vec<&str> = input.split(';').map(str::trim).collect();
                    if fields.len() != 4 {
                        return Err("Enter pKa-list; pH-min; pH-max; step".into());
                    }
                    let pkas = parse_pkas(fields[0])?;
                    let numbers = fields[1..]
                        .iter()
                        .map(|field| {
                            field
                                .replace(',', ".")
                                .parse::<f64>()
                                .map_err(|_| format!("Invalid numeric value: {field}"))
                        })
                        .collect::<Result<Vec<_>, _>>()?;
                    let points =
                        acid_base::bjerrum_points(&pkas, numbers[0], numbers[1], numbers[2])?;
                    self.acid_base_chart = Some(ChemistryChart {
                        title: localized(self.language, "Bjerrum diagram · species distribution")
                            .into_owned(),
                        x_label: "pH".into(),
                        y_label: localized(self.language, "Fraction (%)").into_owned(),
                        series: (0..=pkas.len()).map(|index| format!("α{index}")).collect(),
                        points: points
                            .iter()
                            .map(|(ph, fractions)| {
                                (
                                    *ph,
                                    fractions.iter().map(|fraction| fraction * 100.0).collect(),
                                )
                            })
                            .collect(),
                    });
                    let mut output = format!(
                        "Bjerrum distribution table\npKa: {}\n\n",
                        pkas.iter()
                            .map(|pka| format!("{pka:.3}"))
                            .collect::<Vec<_>>()
                            .join(" | ")
                    );
                    output.push_str("pH");
                    for index in 0..=pkas.len() {
                        output.push_str(&format!("\tα{index} (%)"));
                    }
                    for (ph, fractions) in points {
                        output.push_str(&format!("\n{ph:.2}"));
                        for fraction in fractions {
                            output.push_str(&format!("\t{:.2}", fraction * 100.0));
                        }
                    }
                    output.push_str("\n\nEach α is the fraction of one successive protonation state. Ideal aqueous model; supplied pKa values are not predicted.");
                    Ok(output)
                }
                TextTarget::PolyproticTitrationPh => {
                    let fields: Vec<&str> = input.split(';').map(str::trim).collect();
                    if fields.len() != 5 {
                        return Err(
                            "Enter pKa-list; acid-volume-mL; acid-M; base-M; added-base-volume-mL"
                                .into(),
                        );
                    }
                    let pkas = parse_pkas(fields[0])?;
                    let values = fields[1..]
                        .iter()
                        .map(|field| {
                            field
                                .replace(',', ".")
                                .parse::<f64>()
                                .map_err(|_| format!("Invalid numeric value: {field}"))
                        })
                        .collect::<Result<Vec<_>, _>>()?;
                    let ph = acid_base::polyprotic_titration_ph(
                        &pkas,
                        values[0] / 1000.0,
                        values[1],
                        values[2],
                        values[3] / 1000.0,
                    )?;
                    Ok(format!(
                        "Titration pH = {ph:.4}\n\npKa values: {}\nInitial acid: {:.3} mL × {:.6} mol/L\nAdded strong base: {:.3} mL × {:.6} mol/L\nMethod: polyprotic distribution plus electroneutrality after dilution; ideal aqueous conditions, 25 °C.",
                        pkas.iter().map(|pka| format!("{pka:.3}")).collect::<Vec<_>>().join(", "),
                        values[0], values[1], values[3], values[2]
                    ))
                }
                TextTarget::PolyproticTitrationCurve => {
                    let fields: Vec<&str> = input.split(';').map(str::trim).collect();
                    if fields.len() != 5 {
                        return Err(
                            "Enter pKa-list; acid-volume-mL; acid-M; base-M; max-equivalents"
                                .into(),
                        );
                    }
                    let pkas = parse_pkas(fields[0])?;
                    let values = fields[1..]
                        .iter()
                        .map(|field| {
                            field
                                .replace(',', ".")
                                .parse::<f64>()
                                .map_err(|_| format!("Invalid numeric value: {field}"))
                        })
                        .collect::<Result<Vec<_>, _>>()?;
                    let curve = acid_base::polyprotic_titration_curve(
                        &pkas,
                        values[0] / 1000.0,
                        values[1],
                        values[2],
                        values[3],
                        41,
                    )?;
                    self.acid_base_chart = Some(ChemistryChart {
                        title: localized(self.language, "Polyprotic acid · titration curve")
                            .into_owned(),
                        x_label: localized(self.language, "Added strong base (mL)").into_owned(),
                        y_label: "pH".into(),
                        series: vec!["pH".into()],
                        points: curve
                            .iter()
                            .map(|(volume, ph)| (volume * 1000.0, vec![*ph]))
                            .collect(),
                    });
                    let mut output = format!(
                        "Titration curve points\npKa values: {}\n\nV(base), mL\tpH",
                        pkas.iter()
                            .map(|pka| format!("{pka:.3}"))
                            .collect::<Vec<_>>()
                            .join(", ")
                    );
                    for (volume, ph) in curve {
                        output.push_str(&format!("\n{:.4}\t{ph:.4}", volume * 1000.0));
                    }
                    output.push_str("\n\nMethod: equilibrium charge balance for a polyprotic acid titrated with a strong base; ideal aqueous conditions, 25 °C. Table values can be plotted externally.");
                    Ok(output)
                }
                _ => Err("This is not an acid/base calculation input".into()),
            }
        })();
        self.info_overlay = Some(match result {
            Ok(report) => localize_acid_base_report(self.language, report),
            Err(error) => format!("Acid/base calculation unavailable\n{error}"),
        });
        if self
            .info_overlay
            .as_deref()
            .is_some_and(|report| report.starts_with("Acid/base calculation unavailable"))
        {
            self.acid_base_chart = None;
        }
        self.info_scroll = 0;
    }
    fn apply_molecule_property_edit(&mut self, input: &str) {
        let input = input.trim();
        if let Some(key) = input.strip_prefix('-') {
            let key = key.trim();
            if !valid_property_key(key) {
                self.status = "Use -KEY to remove a valid data key".into();
                return;
            }
            let found = self
                .mol
                .properties
                .iter()
                .any(|(existing, _)| existing.eq_ignore_ascii_case(key))
                || (key.eq_ignore_ascii_case("name") && !self.mol.name.is_empty());
            if !found {
                self.status = format!("No molecule data field named {key}");
                return;
            }
            self.checkpoint();
            self.mol
                .properties
                .retain(|(existing, _)| !existing.eq_ignore_ascii_case(key));
            if key.eq_ignore_ascii_case("name") {
                self.mol.name.clear();
            }
            self.status = format!("Removed molecule data field {key}");
            return;
        }
        let Some((key, value)) = input.split_once('=') else {
            self.status = "Enter KEY=VALUE to save a field, or -KEY to remove it".into();
            return;
        };
        let key = key.trim();
        let value = value.trim();
        if !valid_property_key(key)
            || value.is_empty()
            || value.len() > 2048
            || value.chars().any(char::is_control)
        {
            self.status = "Data key must use letters, digits, '_' or '-'; value must be 1–2048 printable characters".into();
            return;
        }
        self.checkpoint();
        if key.eq_ignore_ascii_case("name") {
            self.mol.name = value.to_string();
        }
        self.mol
            .properties
            .retain(|(existing, _)| !existing.eq_ignore_ascii_case(key));
        self.mol
            .properties
            .push((key.to_string(), value.to_string()));
        self.status = format!("Saved molecule data field {key}");
    }
    fn start_auto_atom_mapping(&mut self, mode: &str) {
        if self.refuse_overlapping_worker() {
            return;
        }
        let selected: BTreeSet<usize> = self.selected_atoms.iter().copied().collect();
        if selected.is_empty() {
            self.info_overlay = Some("Select every complete reactant and product component in one reaction step, then run atom mapping. Leave agents unselected.".into());
            return;
        }
        let chosen: Vec<Vec<usize>> = molecule_components(&self.mol)
            .into_iter()
            .filter(|component| component.iter().any(|atom| selected.contains(atom)))
            .collect();
        if chosen.is_empty()
            || chosen
                .iter()
                .any(|component| component.iter().any(|atom| !selected.contains(atom)))
        {
            self.info_overlay = Some("Automatic mapping requires whole connected components. Select every reactant/product component you want mapped, and leave agents unselected.".into());
            return;
        }
        let mut source = Vec::new();
        let mut target = Vec::new();
        let mut reaction_steps = BTreeSet::new();
        for component in &chosen {
            let roles: BTreeSet<&str> = component
                .iter()
                .map(|&i| self.mol.atoms[i].reaction_role.as_deref().unwrap_or(""))
                .collect();
            let steps: BTreeSet<u16> = component
                .iter()
                .map(|&i| self.mol.atoms[i].reaction_step)
                .collect();
            if roles == BTreeSet::from(["reactant"]) {
                source.extend(component.iter().copied());
            } else if roles == BTreeSet::from(["product"]) {
                target.extend(component.iter().copied());
            } else {
                self.info_overlay = Some("Every selected component must have one explicit role: reactant or product. Leave agents unselected and set roles in Insert → Reaction component.".into());
                return;
            }
            if steps.len() != 1 {
                self.info_overlay =
                    Some("A connected component cannot span multiple reaction steps.".into());
                return;
            }
            reaction_steps.extend(steps);
        }
        if source.is_empty() || target.is_empty() {
            self.info_overlay = Some(
                "Select at least one explicitly assigned reactant and one product component."
                    .into(),
            );
            return;
        }
        if reaction_steps.len() != 1 {
            self.info_overlay = Some(
                "All selected reactant and product components must belong to the same reaction step.".into(),
            );
            return;
        }
        let invalid = |indices: &[usize]| {
            indices.iter().any(|&i| {
                let a = &self.mol.atoms[i];
                a.alias.is_some()
                    || a.query.is_some()
                    || a.r_group.is_some()
                    || !chemistry::valid_element(&a.element)
            }) || self.mol.bonds.iter().any(|b| {
                indices.contains(&b.a)
                    && indices.contains(&b.b)
                    && (b.query.is_some()
                        || b.style != BondStyle::Normal && b.style != BondStyle::Aromatic)
            })
        };
        if source.len() > 999 || target.len() > 999 || invalid(&source) || invalid(&target) {
            self.info_overlay = Some("This optional mapping bridge accepts up to 999 atoms per reaction side. Query atoms, aliases, R-groups and special stereobonds are refused so they are not silently lost.".into());
            return;
        }
        let component_mol = |indices: &[usize]| {
            let mut component = Molecule::default();
            let mut remap = BTreeMap::new();
            for &old in indices {
                remap.insert(old, component.atoms.len());
                component.atoms.push(self.mol.atoms[old].clone());
            }
            component.bonds = self
                .mol
                .bonds
                .iter()
                .filter_map(|bond| {
                    Some(Bond {
                        a: *remap.get(&bond.a)?,
                        b: *remap.get(&bond.b)?,
                        ..bond.clone()
                    })
                })
                .collect();
            component
        };
        let source_mol = component_mol(&source);
        let target_mol = component_mol(&target);
        let source_block = molfile_for(&source_mol);
        let target_block = molfile_for(&target_mol);
        let captured = self.mol.clone();
        let source_indices = source;
        let target_indices = target;
        let operation_mode = mode.to_string();
        let (tx, rx) = mpsc::channel();
        std::thread::spawn(move || {
            let _ = tx.send(rdkit::map_reaction_pair(&source_block, &target_block));
        });
        self.calc_rx = Some((captured, rx));
        self.auto_map_pending = Some((source_indices, target_indices, operation_mode.clone()));
        self.info_overlay = Some(format!("Running optional AtomMapper for all selected components in this reaction step ({operation_mode} maps)…"));
        self.info_scroll = 0;
        self.status = "Mapping selected reaction components…".into();
    }

    fn apply_auto_atom_mapping(
        &mut self,
        report: &str,
        source: &[usize],
        target: &[usize],
        mode: &str,
    ) {
        let mut pairs = Vec::new();
        let mut optimality = "unknown";
        let mut engine = "AtomMapper";
        let mut complete_status = false;
        let mut malformed = false;
        for line in report.lines() {
            let fields: Vec<&str> = line.split('\t').collect();
            match fields.as_slice() {
                ["STATUS", "complete"] if !complete_status => complete_status = true,
                ["OPTIMALITY", value] => optimality = value,
                ["ENGINE", value] => engine = value,
                ["PAIR", a, b, changed] => {
                    match a
                        .parse::<usize>()
                        .ok()
                        .zip(b.parse::<usize>().ok())
                        .zip(changed.parse::<u8>().ok().filter(|value| *value <= 1))
                    {
                        Some(((a, b), changed)) => pairs.push((a, b, changed != 0)),
                        None => malformed = true,
                    }
                }
                [""] => {}
                _ => malformed = true,
            }
        }
        if !complete_status
            || malformed
            || pairs.is_empty()
            || pairs
                .iter()
                .any(|(a, b, _)| *a >= source.len() || *b >= target.len())
        {
            self.info_overlay = Some(
                "Atom mapping returned an incomplete or invalid result. No atom maps were changed."
                    .into(),
            );
            return;
        }
        let mut seen_source = BTreeSet::new();
        let mut seen_target = BTreeSet::new();
        if pairs
            .iter()
            .any(|(a, b, _)| !seen_source.insert(*a) || !seen_target.insert(*b))
        {
            self.info_overlay = Some(
                "Atom mapping returned duplicate atom correspondences. No atom maps were changed."
                    .into(),
            );
            return;
        }
        let selected_atoms: BTreeSet<usize> = source.iter().chain(target).copied().collect();
        let mut map_ids: BTreeSet<u16> = self
            .mol
            .atoms
            .iter()
            .enumerate()
            .filter(|(i, _)| !selected_atoms.contains(i))
            .filter_map(|(_, atom)| atom.atom_map)
            .collect();
        let mut next = 1u16;
        let mut assignments = Vec::new();
        for (a, b, changed) in &pairs {
            if mode == "Changing" && !changed {
                continue;
            }
            while map_ids.contains(&next) && next < 32767 {
                next += 1;
            }
            if next == 0 || map_ids.contains(&next) {
                self.info_overlay = Some("No free atom-map numbers remain (maximum is 32767). No atom maps were changed.".into());
                return;
            }
            let id = next;
            map_ids.insert(id);
            assignments.push((source[*a], id));
            assignments.push((target[*b], id));
        }
        if mode == "Complete" {
            let mapped_source: BTreeSet<usize> = assignments.iter().map(|(i, _)| *i).collect();
            for &atom in source.iter().chain(target) {
                if mapped_source.contains(&atom) {
                    continue;
                }
                while map_ids.contains(&next) && next < 32767 {
                    next += 1;
                }
                if next == 0 || map_ids.contains(&next) {
                    self.info_overlay = Some("No free atom-map numbers remain (maximum is 32767). No atom maps were changed.".into());
                    return;
                }
                map_ids.insert(next);
                assignments.push((atom, next));
            }
        }
        if assignments.is_empty() {
            self.info_overlay = Some(
                "The selected mode produced no atom-map assignments. Existing maps were left unchanged."
                    .into(),
            );
            return;
        }
        self.checkpoint();
        for &atom in &selected_atoms {
            self.mol.atoms[atom].atom_map = None;
        }
        for (atom, id) in &assignments {
            self.mol.atoms[*atom].atom_map = Some(*id);
        }
        self.dirty = true;
        let mapped_pairs = assignments.len() / 2;
        self.info_overlay = Some(format!("Automatic atom mapping · {mode}\n{engine}\nSearch status: complete · optimality: {optimality}\nMapped pairs: {}\n\nAtomMapper maps shared connected cores across the selected reaction sides and may return a seed-conditioned optimum; a completed search is not proof of the globally best mapping. Changing mode marks atoms at engine-reported atom/bond changes. Unmatched atoms are assigned one-sided IDs only in Complete mode. All selected reactant and product components in this step were processed.", if mode == "Complete" { pairs.len() } else { mapped_pairs }));
        self.info_scroll = 0;
        self.status = format!(
            "Applied {} automatic atom-map pair(s)",
            if mode == "Complete" {
                pairs.len()
            } else {
                mapped_pairs
            }
        );
    }

    fn menu_action(&mut self, cmd: &str) {
        match cmd{
        "New"=>self.new_document(),
        "Open"=>self.file_dialog=Some(false),"Save"=>self.file_dialog=Some(true),"Save all"=>self.save_all(),"Recent files"=>self.recent_open=true,"Close tab"=>self.request_close_document(self.active_document),"Undo"=>self.undo(),"Redo"=>self.redo(),
        "Toggle dark theme"=>self.set_dark_theme(!self.dark_theme),
        "Show implicit hydrogens"=>self.set_implicit_hydrogen_visibility(!self.show_implicit_hydrogens),
        "Delete selection"=>self.delete_selection(),
        "Duplicate selection"=>self.duplicate_selection(),
        "Rotate selection 90°"=>self.rotate_selection(std::f32::consts::FRAC_PI_2),
        "Select all"=>self.select_all(),
        "Fit"=>{let w=screen_width();let h=screen_height();let mut points:Vec<Vec2>=self.mol.atoms.iter().map(|a|a.pos).collect();for g in &self.graphics{points.push(match g{Graphic::Arrow{start,end,..}|Graphic::Shape{start,end,..}=>(*start+*end)/2.,Graphic::ElectronArrow{control,..}=>*control,Graphic::Text{pos,..}|Graphic::Plus{pos}=>*pos});}if !points.is_empty(){let(min_x,max_x,min_y,max_y)=points.iter().fold((f32::MAX,f32::MIN,f32::MAX,f32::MIN),|(lx,hx,ly,hy),p|(lx.min(p.x),hx.max(p.x),ly.min(p.y),hy.max(p.y)));let sx=(max_x-min_x).max(BOND);let sy=(max_y-min_y).max(BOND);self.zoom=((w-LEFT-RIGHT-60.)/sx).min((h-TOP-BOTTOM-60.)/sy).clamp(0.2,4.);self.preferred_zoom=self.zoom;self.zoom_preference_dirty=true;self.pan=Vec2::ZERO;}},
        "Zoom in"=>self.adjust_zoom(1.2,vec2(screen_width()/2.,screen_height()/2.),screen_width(),screen_height()),
        "Zoom out"=>self.adjust_zoom(1./1.2,vec2(screen_width()/2.,screen_height()/2.),screen_width(),screen_height()),
        "Elements"=>self.periodic_open=true,"Templates"=>{self.templates_open=true;self.template_tab_custom=false;},
        "Save selection as template"=>{if self.selected.is_none()&&self.selected_atoms.is_empty(){self.status="Select a connected fragment before saving a template".into();}else{self.templates_open=true;self.template_tab_custom=true;self.text_buffer.clear();self.text_target=Some(TextTarget::TemplateName);}},
        "Edit molecule data"=>{self.text_buffer.clear();self.text_target=Some(TextTarget::MoleculeProperty);},
        "Arrow"=>{self.tool=Tool::Arrow;self.arrow_kind="reaction".into()},
        "Line"=>self.tool=Tool::Shape(ShapeKind::Line),
        "Rectangle"=>self.tool=Tool::Shape(ShapeKind::Rectangle),
        "Ellipse"=>self.tool=Tool::Shape(ShapeKind::Ellipse),
        "Equilibrium arrow"=>{self.tool=Tool::Arrow;self.arrow_kind="equilibrium".into()},
        "Resonance arrow"=>{self.tool=Tool::Arrow;self.arrow_kind="resonance".into()},
        "Retrosynthesis arrow"=>{self.tool=Tool::Arrow;self.arrow_kind="retrosynthesis".into()},
        "Electron pair arrow"=>self.tool=Tool::Electron(2),
        "Single electron arrow"=>self.tool=Tool::Electron(1),
        "Plus"=>self.tool=Tool::Plus,"Text"=>self.tool=Tool::Text,
        "Charge +"=>self.tool=Tool::Charge(1),"Charge -"=>self.tool=Tool::Charge(-1),
        "Single"=>self.tool=Tool::Bond(1),"Double"=>self.tool=Tool::Bond(2),"Triple"=>self.tool=Tool::Bond(3),
        "Aromatic"=>self.tool=Tool::StyledBond(BondStyle::Aromatic),"Wedge"=>self.tool=Tool::StyledBond(BondStyle::Wedge),"Hash"=>self.tool=Tool::StyledBond(BondStyle::Hash),"Wavy"=>self.tool=Tool::StyledBond(BondStyle::Wavy),
        "Validate"=>{let a=chemistry::analyze(&self.mol);self.info_overlay=Some(if a.warnings.is_empty(){"No basic valence warnings found. Advanced stereochemistry and structure checks require a chemistry engine.".into()}else{a.warnings.join("\n")});},
        "Find substructure matches"=>{
            let query_atoms=if self.selected_atoms.is_empty(){self.selected.into_iter().collect::<Vec<_>>()}else{self.selected_atoms.clone()};
            match query::find_selected_matches(&self.mol,&query_atoms){
                Ok(report)=>{
                    if report.matches.is_empty(){self.info_overlay=Some(format!("No matches. Query: {} atoms; target: {} atoms. The selected query component is excluded from the target.",report.query_atoms,report.target_atoms));}
                    else{
                        let mut atoms=BTreeSet::new();
                        for found in &report.matches{atoms.extend(found.iter().copied());}
                        self.selected_atoms=atoms.into_iter().collect();self.selected=self.selected_atoms.first().copied();
                        let lines=report.matches.iter().take(30).enumerate().map(|(i,found)|format!("{}. target atoms {}",i+1,found.iter().map(|n|(n+1).to_string()).collect::<Vec<_>>().join(", "))).collect::<Vec<_>>().join("\n");
                        self.info_overlay=Some(format!("{} match(es) · query {} atoms · target {} atoms\n\n{}{}\n\nSelected atoms highlight the matching occurrences. This local matcher supports element/list, aromaticity, hydrogen, connectivity, degree, valence, unsaturation, substitution count, ring-bond count, minimum ring size, deterministic minimum-cycle-basis ring count, bond-order predicates and exact stored R/S and E/Z/cis/trans labels. Ring counts may differ from other toolkits when cycle bases tie; wedge/hash stereochemistry and stereo groups are unsupported.",report.matches.len(),report.query_atoms,report.target_atoms,lines,if report.matches.len()>30{"\n…showing first 30"}else{""}));
                    }
                },
                Err(error)=>self.info_overlay=Some(error),
            }
        },
        "Name structure"=>{self.naming_highlight.clear();self.info_overlay=Some(match naming::name(&self.mol){Ok(n)=>{self.naming_highlight=n.highlighted_atoms.clone();format!("Teaching name: {}\n\n{}",n.name,n.explanation)},Err(e)=>format!("Name unavailable\n{e}")});},
        "Clean 2D"=>self.clean_2d(),
        "Generate 3D coordinates"=>{
            if self.refuse_overlapping_worker() { return; }
            if self.mol.atoms.is_empty() {
                self.info_overlay=Some("Draw or open a molecule before generating 3D coordinates.".into());
            } else if self.mol.atoms.iter().any(|a|a.query.is_some()||a.r_group.is_some()||!chemistry::valid_element(&a.element)) {
                self.info_overlay=Some("This structure contains query atoms or nonstandard elements that cannot be sent to the RDKit 3D worker.".into());
            } else {
                let block=self.molfile();
                let captured=self.mol.clone();
                let (tx,rx)=mpsc::channel();
                std::thread::spawn(move||{let _=tx.send(rdkit::generate_3d(&block));});
                self.calc_rx=Some((captured,rx));
                self.info_overlay=Some("Generating 3D coordinates with optional RDKit…".into());
            }
        },
        "3D viewer"=>{if self.mol.atoms.iter().all(|a|a.coord3d.is_some())&&!self.mol.atoms.is_empty(){self.view_3d_open=true;}else{self.info_overlay=Some("No 3D coordinates are stored. Choose Structure → Generate 3D coordinates; that action requires optional RDKit.".into());}},
        "Add explicit hydrogens"=>self.add_explicit_hydrogens(),
        "Remove explicit hydrogens"=>self.remove_explicit_hydrogens(),
        "Rotate 90 CW"=>self.transform_document(std::f32::consts::FRAC_PI_2,false,false),
        "Rotate 90 CCW"=>self.transform_document(-std::f32::consts::FRAC_PI_2,false,false),
        "Rotate 180"=>self.transform_document(std::f32::consts::PI,false,false),
        "Flip left-right"=>self.transform_document(0.,true,false),
        "Flip up-down"=>self.transform_document(0.,false,true),
        "Benzene template"=>self.activate_builtin_template("Benzene"),
        "Formula & mass"=>{let a=chemistry::analyze(&self.mol);self.info_overlay=Some(format!("Formula: {}\nMolar mass: {} g/mol\nExact mass: {} Da\nMethod: local valence/H count; CIAAW standard atomic weights (2024), interval midpoint estimates; listed isotope masses and ionic electron correction.\nUnusual valence, query chemistry and unsupported elements may be unresolved.",a.formula,a.molar_mass.map(|x|format!("{x:.3}")).unwrap_or_else(||"unavailable".into()),a.exact_mass.map(|x|format!("{x:.4}")).unwrap_or_else(||"unavailable".into())));},
        "Molecular properties"=>{
            let a=chemistry::analyze(&self.mol);
            let comp=a.composition.iter().map(|(el,n)|format!("{el}: {n}")).collect::<Vec<_>>().join(", ");
            self.info_overlay=Some(format!("Local graph properties\nFormula: {}\nElements: {}\nMolar mass: {} g/mol\nExact mass: {} Da\nAtoms: {} (heavy: {})\nBonds: {}\nConnected components: {}\nCycle rank (E − V + C): {}\nAverage graph degree: {:.2}\nFormal charge: {}\nMethod: graph counts; CIAAW 2024 standard atomic weights (interval midpoints) and explicit isotope masses. Exact ionic mass applies electron mass correction. Cycle rank includes every stored bond, including coordination links.\nlogP, TPSA, pKa, NMR and 3D require a verified chemistry engine.",a.formula,comp,a.molar_mass.map(|x|format!("{x:.3}")).unwrap_or_else(||"unavailable".into()),a.exact_mass.map(|x|format!("{x:.4}")).unwrap_or_else(||"unavailable".into()),a.atom_count,a.heavy_atom_count,a.bond_count,a.connected_components,a.cycle_rank,a.average_degree,a.formal_charge));
        },
        "Local topological descriptors"=>{
            self.info_overlay=Some(match chemistry::topological_descriptors(&self.mol) {
                Ok(d)=>format!("Classical heavy-atom graph descriptors\nHeavy atoms: {}\nHeavy-atom bonds: {}\nWiener index: {}\nGraph diameter: {} bonds\nFirst Zagreb index (M₁): {}\nSecond Zagreb index (M₂): {}\nRandić connectivity index (α = −½): {:.4}\nMethod: shortest-path sum and degree products on the connected heavy-atom graph; explicit hydrogens excluded. These are graph descriptors, not physical property predictions.",d.heavy_atoms,d.heavy_atom_bonds,d.wiener_index,d.diameter,d.first_zagreb,d.second_zagreb,d.randic_connectivity),
                Err(error)=>format!("Topological descriptors unavailable\n{error}"),
            });
        },
        "Reaction analysis"=>self.analyze_reaction(),
        "Auto-map reaction (Complete)"=>self.start_auto_atom_mapping("Complete"),
        "Auto-map reaction (Matching)"=>self.start_auto_atom_mapping("Matching"),
        "Auto-map reaction (Changing)"=>self.start_auto_atom_mapping("Changing"),
        "Gasteiger partial charges"=>self.start_rdkit_partial_charges(),
        "Auto-balance reaction"=>self.auto_balance_reaction(),
        "Auto-balance all reaction steps"=>self.auto_balance_all_reactions(),
        "View all reaction steps"=>self.reaction_step_view=None,
        "View previous reaction step"=>{
            let count=self.graphics.iter().filter(|g|matches!(g,Graphic::Arrow{kind,..} if kind=="reaction"||kind=="equilibrium"||kind=="retrosynthesis")).count() as u16;
            if count>0 { let current=self.reaction_step_view.unwrap_or(1).min(count); self.reaction_step_view=Some(current.saturating_sub(1).max(1)); }
        },
        "View next reaction step"=>{
            let count=self.graphics.iter().filter(|g|matches!(g,Graphic::Arrow{kind,..} if kind=="reaction"||kind=="equilibrium"||kind=="retrosynthesis")).count() as u16;
            if count>0 { let current=self.reaction_step_view.unwrap_or(1).min(count); self.reaction_step_view=Some((current+1).min(count)); }
        },
        "View reaction step 1"=>self.reaction_step_view=Some(1),
        "View reaction step 2"=>self.reaction_step_view=Some(2),
        "View reaction step 3"=>self.reaction_step_view=Some(3),
        "Reaction component: reactant"=>self.set_reaction_component(Some("reactant"),None,None),
        "Reaction component: product"=>self.set_reaction_component(Some("product"),None,None),
        "Reaction component: agent"=>self.set_reaction_component(Some("agent"),None,None),
        "Reaction component: automatic role"=>self.set_reaction_component(Some(""),None,None),
        "Reaction coefficient: 2"=>self.set_reaction_component(None,None,Some(2)),
        "Reaction coefficient: 3"=>self.set_reaction_component(None,None,Some(3)),
        "Reaction coefficient: custom"=>{self.text_buffer.clear();self.text_target=Some(TextTarget::ReactionCoefficient);},
        "Chemistry A acid/base tools"=>self.acid_base_open=true,
        "Reaction step: 1"=>self.set_reaction_component(None,Some(1),None),
        "Reaction step: 2"=>self.set_reaction_component(None,Some(2),None),
        "Reaction step: 3"=>self.set_reaction_component(None,Some(3),None),
        "RDKit descriptors"=>{
            if self.refuse_overlapping_worker(){return;}
            if self.mol.atoms.is_empty(){self.info_overlay=Some("Draw a molecule first".into());}
            else if self.mol.atoms.len()>999||self.mol.bonds.len()>999 {self.info_overlay=Some("RDKit bridge currently accepts MOL V2000 (maximum 999 atoms/bonds)".into());}
            else if self.mol.atoms.iter().any(|a|a.query.is_some()||a.r_group.is_some()||a.alias.is_some()||!chemistry::valid_element(&a.element))||self.mol.bonds.iter().any(|b|b.style!=BondStyle::Normal||b.query.is_some()) {self.info_overlay=Some("This structure contains query, R-group, alias, or special bond data that the RDKit MOL bridge cannot preserve. No descriptor was calculated.".into());}
            else {let block=self.molfile();let captured=self.mol.clone();let (tx,rx)=mpsc::channel();std::thread::spawn(move||{let _=tx.send(rdkit::descriptors(&block));});self.calc_rx=Some((captured,rx));self.info_overlay=Some("Calculating RDKit descriptors…".into());self.status="Calculating with RDKit...".into();}
        },
        "Enumerate tautomers" => self.start_rdkit_chemical_report("tautomer-report"),
        "Enumerate stereoisomers" => self.start_rdkit_chemical_report("stereoisomer-report"),
        "Fetch PubChem properties" => self.start_pubchem_lookup(),
        "Settings"=>self.settings_open=true,
        "Command palette"=>{self.palette_open=true;self.palette_query.clear();self.palette_selection=0;},
        "Search PubChem in browser"=>{let query=self.mol.name.trim().to_string();let query=if query.is_empty(){self.mol.properties.iter().find(|(key,_)|key.eq_ignore_ascii_case("name")).map(|(_,value)|value.clone()).unwrap_or_else(||chemistry::analyze(&self.mol).formula)}else{query};let url=format!("https://pubchem.ncbi.nlm.nih.gov/#query={}",url_encode(&query));self.status=match open_external_url(&url){Ok(())=>format!("Opened PubChem browser search for {query} (browser search only; no data was retrieved)"),Err(error)=>error};},
        "Naming lesson"=>self.start_naming_lesson(),
        "Reaction balancing lesson"=>self.start_reaction_lesson(),
        "Aromatize"=>self.start_rdkit_graph_operation("aromatize"),
        "Kekulé form"=>self.start_rdkit_graph_operation("kekulize"),
        "Stereochemistry report"=>self.start_stereo_report(),
        "Assign R/S labels with RDKit"=>self.start_stereo_assignment(),
        "Copy as SMILES"=>self.start_copy_smiles(),
        "Copy as SMARTS"=>self.start_copy_smarts(),
        "Copy drawing as image"=>{
            self.status=match render_png(&self.mol,&self.graphics,1920,1080)
                .and_then(|png|copy_png_to_clipboard(&png))
            {
                Ok(())=>"Copied the 1920×1080 drawing image to the clipboard".into(),
                Err(error)=>format!("Could not copy drawing image: {error}"),
            };
        },
        "Resume lesson"=>{if self.lesson.is_some(){self.lesson_open=true;}else{self.status="Start a lesson from Help first".into();}},
        "Integration status"=>self.info_overlay=Some("PubChem direct structure lookup uses the official NCBI PUG REST API and requires an internet connection plus curl. It sends a locally representable exact structure only after you select Fetch PubChem properties; the drawing is not changed.\n\nRDKit is an optional external Python integration. Without it, 3D generation, local advanced descriptors, stereochemistry assignment and broad format conversion are unavailable. Local drawing, MOL and MOLS still work.".into()),
        "Shortcuts"=>{let m=if cfg!(target_os="macos"){"Cmd"}else{"Ctrl"};self.info_overlay=Some(format!("{m}+N: New | {m}+O: Open | {m}+S: Save | {m}+Z: Undo | {m}+Y: Redo\n{m}+L: Copy as SMILES | {m}+Shift+K: Copy as SMARTS | {m}+Shift+P: Command palette | Delete: delete selection\nEsc: select tool | Shift+R: six-membered ring\nDrawing buffer: 1/2/3 orders, 4 aromatic, 5 wedge, 6 hash, 7 wavy, 12/14/24 bond queries\nC, N, O, Cl... select element; + and - change selected charge\n@13 sets isotope mass 13 on the selected atom; @0 clears it\nlabel:Ph assigns alias Ph; label:- clears the alias\nR1...R32767 assigns R-group; M1... assigns a map; M0 clears; M+/M= assign unique maps\n{m}+E opens the periodic table."))},
        "Quick start guide"=>self.info_overlay=Some(quick_start_guide(self.language)),
        "Feature status"=>self.info_overlay=Some("See STATUS.md in the source folder for implemented, partial and unavailable functions.".into()),
        _=>{}
    }
    }
    fn apply_quick(&mut self, a: shortcuts::QuickAction) {
        use shortcuts::QuickAction::*;
        match a {
            Bond(n) => {
                self.bond_query = None;
                self.tool = match n {
                    1..=3 => Tool::Bond(n),
                    4 => Tool::StyledBond(BondStyle::Aromatic),
                    5 => Tool::StyledBond(BondStyle::Wedge),
                    6 => Tool::StyledBond(BondStyle::Hash),
                    7 => Tool::StyledBond(BondStyle::Wavy),
                    _ => Tool::StyledBond(BondStyle::Any),
                };
            }
            BondQuery(q) => {
                self.tool = Tool::StyledBond(BondStyle::Any);
                self.bond_query = Some(q.into());
                self.status = format!("Query bond: {q}");
            }
            Element(s) => {
                self.element = s.clone();
                if let Some(i) = self.selected {
                    self.checkpoint();
                    self.mol.atoms[i].element = s;
                } else {
                    self.tool = Tool::Atom;
                }
            }
            Charge(n) => {
                if let Some(i) = self.selected {
                    self.checkpoint();
                    self.mol.atoms[i].charge = (self.mol.atoms[i].charge + n).clamp(-9, 9);
                } else {
                    self.tool = Tool::Charge(n);
                }
            }
            RGroup(n) => {
                if let Some(i) = self.selected {
                    self.checkpoint();
                    self.mol.atoms[i].r_group = Some(n);
                    self.mol.atoms[i].element = "R".into();
                } else {
                    self.status = "Select an atom before assigning an R-group".into();
                }
            }
            AtomMap(n) => {
                if let Some(i) = self.selected {
                    self.checkpoint();
                    self.mol.atoms[i].atom_map = if n == 0 { None } else { Some(n) };
                } else {
                    self.status = "Select an atom before assigning a map number".into();
                }
            }
            AssignUniqueAtomMaps => self.assign_unique_atom_maps(),
            Isotope(mass) => {
                if let Some(i) = self.selected.filter(|i| *i < self.mol.atoms.len()) {
                    if self.mol.atoms[i].isotope != mass {
                        self.checkpoint();
                        self.mol.atoms[i].isotope = mass;
                    }
                } else {
                    self.status =
                        "Select an atom before assigning an isotope (type @0 to clear)".into();
                }
            }
            Alias(alias) => {
                if let Some(i) = self.selected.filter(|i| *i < self.mol.atoms.len()) {
                    if self.mol.atoms[i].alias != alias {
                        self.checkpoint();
                        self.mol.atoms[i].alias = alias;
                    }
                } else {
                    self.status = "Select an atom before assigning an alias (type label:Ph)".into();
                }
            }
            Query(q) => {
                if let Some(i) = self.selected {
                    match shortcuts::merge_query_terms(self.mol.atoms[i].query.as_deref(), &q) {
                        Ok(merged) => {
                            if self.mol.atoms[i].query.as_deref() != Some(merged.as_str()) {
                                self.checkpoint();
                                self.mol.atoms[i].query = Some(merged);
                            }
                        }
                        Err(error) => {
                            self.status = error;
                        }
                    }
                } else {
                    self.status = "Select an atom before assigning a query".into();
                }
            }
        }
    }
    fn molfile(&self) -> String {
        molfile_for(&self.mol)
    }
    fn svg_with_size(&self, width: f32, height: f32) -> String {
        let mut s = format!(
            "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\" width=\"{}\" height=\"{}\"><rect width=\"100%\" height=\"100%\" fill=\"white\"/>",
            width as i32, height as i32, width as i32, height as i32
        );
        for b in &self.mol.bonds {
            let a = self.mol.atoms[b.a].pos;
            let c = self.mol.atoms[b.b].pos;
            for primitive in export_bond_primitives(a, c, b, 1.) {
                match primitive {
                    BondPrimitive::Line(p, q) => {
                        s += &format!("<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#222\" stroke-width=\"2\"/>", p.x, p.y, q.x, q.y);
                    }
                    BondPrimitive::Triangle(p1, p2, p3) => {
                        s += &format!(
                            "<polygon points=\"{},{} {},{} {},{}\" fill=\"#222\"/>",
                            p1.x, p1.y, p2.x, p2.y, p3.x, p3.y
                        );
                    }
                }
            }
        }
        for (i, a) in self.mol.atoms.iter().enumerate() {
            if show_label(&self.mol, i) {
                let label = format!(
                    "{}{}{}",
                    a.isotope.map(|n| n.to_string()).unwrap_or_default(),
                    a.alias.as_ref().unwrap_or(&a.element),
                    charge_label(a.charge)
                );
                s+=&format!("<text x=\"{}\" y=\"{}\" font-family=\"Arial\" font-size=\"17\" text-anchor=\"middle\" fill=\"#222\">{}</text>",a.pos.x,a.pos.y+5.,xml_escape(&label));
            }
        }
        for g in &self.graphics {
            match g {
                Graphic::Arrow {
                    start,
                    end,
                    kind,
                    above,
                    below,
                } => {
                    s+=&format!("<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#222\" stroke-width=\"2\"/>",start.x,start.y,end.x,end.y);
                    let d = (*end - *start).normalize_or_zero();
                    let q = vec2(-d.y, d.x);
                    s += &format!(
                        "<polygon points=\"{},{} {},{} {},{}\" fill=\"#222\"/>",
                        end.x,
                        end.y,
                        end.x - d.x * 12. + q.x * 5.,
                        end.y - d.y * 12. + q.y * 5.,
                        end.x - d.x * 12. - q.x * 5.,
                        end.y - d.y * 12. - q.y * 5.
                    );
                    if kind == "equilibrium" {
                        s+=&format!("<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#222\" stroke-width=\"1.5\"/>",start.x+q.x*6.,start.y+q.y*6.,end.x+q.x*6.,end.y+q.y*6.);
                    }
                    if kind == "resonance" {
                        s += &format!(
                            "<polygon points=\"{},{} {},{} {},{}\" fill=\"#222\"/>",
                            start.x,
                            start.y,
                            start.x + d.x * 12. + q.x * 5.,
                            start.y + d.y * 12. + q.y * 5.,
                            start.x + d.x * 12. - q.x * 5.,
                            start.y + d.y * 12. - q.y * 5.
                        );
                    }
                    if kind == "retrosynthesis" {
                        s+=&format!("<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#222\" stroke-width=\"2\"/>",end.x+q.x*7.,end.y+q.y*7.,end.x-q.x*7.,end.y-q.y*7.);
                    }
                    let mid = (*start + *end) / 2.;
                    s += &format!(
                        "<text x=\"{}\" y=\"{}\">{}</text><text x=\"{}\" y=\"{}\">{}</text>",
                        mid.x,
                        mid.y - 12.,
                        xml_escape(above),
                        mid.x,
                        mid.y + 25.,
                        xml_escape(below)
                    );
                }
                Graphic::Text {
                    pos, text, size, ..
                } => {
                    s += &format!(
                        "<text x=\"{}\" y=\"{}\" font-size=\"{}\">{}</text>",
                        pos.x,
                        pos.y,
                        size,
                        xml_escape(text)
                    )
                }
                Graphic::Plus { pos } => {
                    s += &format!(
                        "<text x=\"{}\" y=\"{}\" font-size=\"30\">+</text>",
                        pos.x, pos.y
                    )
                }
                Graphic::Shape { kind, start, end } => {
                    let x = start.x.min(end.x);
                    let y = start.y.min(end.y);
                    let w = (start.x - end.x).abs();
                    let h = (start.y - end.y).abs();
                    match kind {
                        ShapeKind::Line=>s+=&format!("<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#222\" stroke-width=\"2\"/>",start.x,start.y,end.x,end.y),
                        ShapeKind::Rectangle=>s+=&format!("<rect x=\"{x}\" y=\"{y}\" width=\"{w}\" height=\"{h}\" fill=\"none\" stroke=\"#222\" stroke-width=\"2\"/>"),
                        ShapeKind::Ellipse=>s+=&format!("<ellipse cx=\"{}\" cy=\"{}\" rx=\"{}\" ry=\"{}\" fill=\"none\" stroke=\"#222\" stroke-width=\"2\"/>",(start.x+end.x)/2.,(start.y+end.y)/2.,w/2.,h/2.),
                    }
                }
                Graphic::ElectronArrow {
                    start,
                    control,
                    end,
                    electrons,
                } => {
                    s+=&format!("<path d=\"M {} {} Q {} {} {} {}\" fill=\"none\" stroke=\"#222\" stroke-width=\"2\"/>",start.x,start.y,control.x,control.y,end.x,end.y);
                    let d = (*end - *control).normalize_or_zero();
                    let q = vec2(-d.y, d.x);
                    s+=&format!("<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#222\" stroke-width=\"2\"/>",end.x,end.y,end.x-d.x*12.+q.x*5.,end.y-d.y*12.+q.y*5.);
                    if *electrons == 2 {
                        s+=&format!("<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"#222\" stroke-width=\"2\"/>",end.x,end.y,end.x-d.x*12.-q.x*5.,end.y-d.y*12.-q.y*5.);
                    }
                }
            }
        }
        s += "</svg>";
        s
    }
    fn write_svg(&self, path: &str, width: f32, height: f32) -> std::io::Result<()> {
        fs::write(path, self.svg_with_size(width, height))
    }
}
fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}
fn export_bounds(mol: &Molecule, graphics: &[Graphic]) -> (f32, f32, f32, f32) {
    let mut points: Vec<Vec2> = mol.atoms.iter().map(|a| a.pos).collect();
    for g in graphics {
        match g {
            Graphic::Arrow { start, end, .. } | Graphic::Shape { start, end, .. } => {
                points.push(*start);
                points.push(*end)
            }
            Graphic::ElectronArrow {
                start,
                control,
                end,
                ..
            } => {
                points.push(*start);
                points.push(*control);
                points.push(*end)
            }
            Graphic::Text { pos, .. } | Graphic::Plus { pos } => points.push(*pos),
        }
    }
    if points.is_empty() {
        return (0., 100., 0., 100.);
    }
    let bounds = points.iter().fold(
        (f32::MAX, f32::MIN, f32::MAX, f32::MIN),
        |(lx, hx, ly, hy), p| (lx.min(p.x), hx.max(p.x), ly.min(p.y), hy.max(p.y)),
    );
    let pad = 50.;
    (
        bounds.0 - pad,
        bounds.1 + pad,
        bounds.2 - pad,
        bounds.3 + pad,
    )
}
fn export_transform(p: Vec2, b: (f32, f32, f32, f32), width: f32, height: f32) -> Vec2 {
    let scale = ((width - 80.) / (b.1 - b.0).max(1.)).min((height - 80.) / (b.3 - b.2).max(1.));
    vec2((p.x - b.0) * scale + 40., (p.y - b.2) * scale + 40.)
}
#[derive(Clone, Copy, Debug, PartialEq)]
enum BondPrimitive {
    Line(Vec2, Vec2),
    Triangle(Vec2, Vec2, Vec2),
}
fn export_bond_primitives(a: Vec2, b: Vec2, bond: &Bond, scale: f32) -> Vec<BondPrimitive> {
    let d = (b - a).normalize_or_zero();
    let q = vec2(-d.y, d.x);
    let unit = scale.max(0.01);
    let line = |start, end| BondPrimitive::Line(start, end);
    match bond.style {
        BondStyle::Normal => {
            let count = bond.order.clamp(1, 3);
            (0..count)
                .map(|k| {
                    let offset = (k as f32 - (count as f32 - 1.) / 2.) * 5. * unit;
                    line(a + q * offset, b + q * offset)
                })
                .collect()
        }
        BondStyle::Aromatic => {
            let mut primitives = vec![line(a, b)];
            let offset = 5. * unit;
            for k in 0..8 {
                let start = a.lerp(b, k as f32 / 8.) + q * offset;
                let end = a.lerp(b, (k as f32 + 0.55) / 8.) + q * offset;
                primitives.push(line(start, end));
            }
            primitives
        }
        BondStyle::Wedge => vec![BondPrimitive::Triangle(
            a,
            b + q * (6. * unit),
            b - q * (6. * unit),
        )],
        BondStyle::Hash => (0..11)
            .map(|k| {
                let t = (k as f32 + 1.) / 11.;
                let center = a.lerp(b, t);
                let half_width = 6. * unit * t;
                line(center - q * half_width, center + q * half_width)
            })
            .collect(),
        BondStyle::Wavy => {
            let mut primitives = Vec::with_capacity(24);
            let mut previous = a;
            for k in 1..=24 {
                let t = k as f32 / 24.;
                let current = a.lerp(b, t) + q * (t * std::f32::consts::TAU * 4.).sin() * 3. * unit;
                primitives.push(line(previous, current));
                previous = current;
            }
            primitives
        }
        BondStyle::Coordinate => {
            let head = 10. * unit;
            let half = 4. * unit;
            vec![
                line(a, b),
                line(b, b - d * head + q * half),
                line(b, b - d * head - q * half),
            ]
        }
        BondStyle::Any => (0..8)
            .map(|k| {
                let start = a.lerp(b, k as f32 / 8.);
                let end = a.lerp(b, (k as f32 + 0.55) / 8.);
                line(start, end)
            })
            .collect(),
    }
}
fn raster_line(image: &mut image::RgbaImage, a: Vec2, b: Vec2, color: image::Rgba<u8>, thick: i32) {
    let (mut x0, mut y0, x1, y1) = (
        a.x.round() as i32,
        a.y.round() as i32,
        b.x.round() as i32,
        b.y.round() as i32,
    );
    let dx = (x1 - x0).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let dy = -(y1 - y0).abs();
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy;
    loop {
        for ox in -thick..=thick {
            for oy in -thick..=thick {
                let x = x0 + ox;
                let y = y0 + oy;
                if x >= 0 && y >= 0 && (x as u32) < image.width() && (y as u32) < image.height() {
                    image.put_pixel(x as u32, y as u32, color);
                }
            }
        }
        if x0 == x1 && y0 == y1 {
            break;
        }
        let e2 = 2 * err;
        if e2 >= dy {
            err += dy;
            x0 += sx
        }
        if e2 <= dx {
            err += dx;
            y0 += sy
        }
    }
}
fn raster_text(image: &mut image::RgbaImage, text: &str, at: Vec2, size: f32) {
    let Ok(font) = fontdue::Font::from_bytes(
        include_bytes!("../assets/DejaVuSans.ttf") as &[u8],
        fontdue::FontSettings::default(),
    ) else {
        return;
    };
    let color = [34, 43, 38, 255];
    let mut x = at.x.round() as i32;
    let baseline = at.y.round() as i32;
    for ch in text.chars() {
        let (metrics, bitmap) = font.rasterize(ch, size.max(7.));
        for gy in 0..metrics.height {
            for gx in 0..metrics.width {
                let alpha = bitmap[gy * metrics.width + gx];
                let px = x + metrics.xmin + gx as i32;
                let py = baseline - metrics.ymin - metrics.height as i32 + gy as i32;
                if alpha > 0
                    && px >= 0
                    && py >= 0
                    && (px as u32) < image.width()
                    && (py as u32) < image.height()
                {
                    image.put_pixel(
                        px as u32,
                        py as u32,
                        image::Rgba([color[0], color[1], color[2], alpha]),
                    );
                }
            }
        }
        x += metrics.advance_width.round() as i32;
    }
}
fn render_png(
    mol: &Molecule,
    graphics: &[Graphic],
    width: u32,
    height: u32,
) -> Result<Vec<u8>, String> {
    let (width, height) = (width.clamp(320, 4096), height.clamp(240, 4096));
    let mut image = image::RgbaImage::from_pixel(width, height, image::Rgba([255, 255, 255, 255]));
    let bounds = export_bounds(mol, graphics);
    let sx = width as f32 / 1920.;
    let sy = height as f32 / 1080.;
    let point = |p: Vec2| export_transform(p, bounds, width as f32, height as f32);
    for bond in &mol.bonds {
        let Some(a) = mol.atoms.get(bond.a) else {
            continue;
        };
        let Some(c) = mol.atoms.get(bond.b) else {
            continue;
        };
        let aa = point(a.pos);
        let cc = point(c.pos);
        for primitive in export_bond_primitives(aa, cc, bond, sx.min(sy)) {
            match primitive {
                BondPrimitive::Line(start, end) => {
                    raster_line(&mut image, start, end, image::Rgba([45, 53, 48, 255]), 1);
                }
                BondPrimitive::Triangle(p1, p2, p3) => {
                    let steps = ((p2 - p1).length().max((p3 - p1).length()) * 1.5).ceil() as usize;
                    for step in 0..=steps {
                        let t = step as f32 / steps.max(1) as f32;
                        raster_line(
                            &mut image,
                            p1.lerp(p2, t),
                            p1.lerp(p3, t),
                            image::Rgba([45, 53, 48, 255]),
                            0,
                        );
                    }
                }
            }
        }
    }
    for (i, a) in mol.atoms.iter().enumerate() {
        if show_label(mol, i) {
            let label = format!(
                "{}{}{}",
                a.isotope.map(|n| n.to_string()).unwrap_or_default(),
                a.alias.as_ref().unwrap_or(&a.element),
                charge_label(a.charge)
            );
            let at = point(a.pos);
            raster_text(
                &mut image,
                &label,
                at - vec2(label.len() as f32 * 3. * sx, 0.),
                18. * sx.max(sy),
            );
        }
    }
    for g in graphics {
        match g {
            Graphic::Arrow {
                start,
                end,
                kind,
                above,
                below,
            } => {
                let a = point(*start);
                let b = point(*end);
                raster_line(&mut image, a, b, image::Rgba([42, 52, 46, 255]), 1);
                if kind != "resonance" {
                    let d = (b - a).normalize_or_zero();
                    let q = vec2(-d.y, d.x);
                    raster_line(
                        &mut image,
                        b,
                        b - d * 12. * sx + q * 6. * sy,
                        image::Rgba([42, 52, 46, 255]),
                        1,
                    );
                    raster_line(
                        &mut image,
                        b,
                        b - d * 12. * sx - q * 6. * sy,
                        image::Rgba([42, 52, 46, 255]),
                        1,
                    );
                }
                let mid = (a + b) / 2.;
                raster_text(
                    &mut image,
                    above,
                    mid - vec2(20. * sx, 12. * sy),
                    15. * sx.max(sy),
                );
                raster_text(
                    &mut image,
                    below,
                    mid + vec2(-12. * sx, 20. * sy),
                    15. * sx.max(sy),
                );
            }
            Graphic::Text {
                pos, text, size, ..
            } => raster_text(&mut image, text, point(*pos), *size * sx),
            Graphic::Plus { pos } => raster_text(&mut image, "+", point(*pos), 30. * sx),
            Graphic::Shape { kind, start, end } => {
                let a = point(*start);
                let b = point(*end);
                let color = image::Rgba([42, 52, 46, 255]);
                match kind {
                    ShapeKind::Line => raster_line(&mut image, a, b, color, 1),
                    ShapeKind::Rectangle => {
                        let c = vec2(b.x, a.y);
                        let d = vec2(a.x, b.y);
                        for (u, v) in [(a, c), (c, b), (b, d), (d, a)] {
                            raster_line(&mut image, u, v, color, 1);
                        }
                    }
                    ShapeKind::Ellipse => {
                        let mid = (a + b) / 2.;
                        let radius = (b - a) / 2.;
                        let mut prev = mid + vec2(radius.x, 0.);
                        for i in 1..=48 {
                            let t = i as f32 / 48. * std::f32::consts::TAU;
                            let cur = mid + vec2(radius.x * t.cos(), radius.y * t.sin());
                            raster_line(&mut image, prev, cur, color, 1);
                            prev = cur;
                        }
                    }
                }
            }
            Graphic::ElectronArrow {
                start,
                control,
                end,
                ..
            } => {
                let mut prev = point(*start);
                for i in 1..=24 {
                    let t = i as f32 / 24.;
                    let cur = point(*start) * (1. - t).powi(2)
                        + point(*control) * (2. * t * (1. - t))
                        + point(*end) * t * t;
                    raster_line(&mut image, prev, cur, image::Rgba([42, 52, 46, 255]), 1);
                    prev = cur;
                }
            }
        }
    }
    let mut bytes = std::io::Cursor::new(Vec::new());
    image::DynamicImage::ImageRgba8(image)
        .write_to(&mut bytes, image::ImageFormat::Png)
        .map_err(|e| e.to_string())?;
    Ok(bytes.into_inner())
}
fn render_pdf(mol: &Molecule, graphics: &[Graphic], width: f32, height: f32) -> Vec<u8> {
    let bounds = export_bounds(mol, graphics);
    let point = |p: Vec2| {
        let q = export_transform(p, bounds, width, height);
        vec2(q.x, height - q.y)
    };
    let mut stream = String::from("0.12 0.16 0.13 RG 1.5 w\n");
    for b in &mol.bonds {
        if let (Some(a), Some(c)) = (mol.atoms.get(b.a), mol.atoms.get(b.b)) {
            let aa = point(a.pos);
            let cc = point(c.pos);
            for primitive in export_bond_primitives(aa, cc, b, 1.) {
                match primitive {
                    BondPrimitive::Line(start, end) => {
                        stream += &format!(
                            "{:.2} {:.2} m {:.2} {:.2} l S\n",
                            start.x, start.y, end.x, end.y
                        );
                    }
                    BondPrimitive::Triangle(p1, p2, p3) => {
                        stream += &format!("0.12 0.16 0.13 rg\n{:.2} {:.2} m {:.2} {:.2} l {:.2} {:.2} l h f\n0.12 0.16 0.13 RG\n", p1.x, p1.y, p2.x, p2.y, p3.x, p3.y);
                    }
                }
            }
        }
    }
    for (i, a) in mol.atoms.iter().enumerate() {
        if show_label(mol, i) {
            let label = format!(
                "{}{}{}",
                a.isotope.map(|n| n.to_string()).unwrap_or_default(),
                a.alias.as_ref().unwrap_or(&a.element),
                charge_label(a.charge)
            );
            let p = point(a.pos);
            stream += &format!(
                "BT /F1 12 Tf {:.2} {:.2} Td ({}) Tj ET\n",
                p.x - 4. * label.len() as f32,
                p.y - 3.,
                pdf_escape(&label)
            );
        }
    }
    for g in graphics {
        match g {
            Graphic::Arrow { start, end, .. } => {
                let a = point(*start);
                let b = point(*end);
                stream += &format!("{:.2} {:.2} m {:.2} {:.2} l S\n", a.x, a.y, b.x, b.y);
                let d = (b - a).normalize_or_zero();
                let q = vec2(-d.y, d.x);
                for tip in [b - d * 12. + q * 5., b - d * 12. - q * 5.] {
                    stream += &format!("{:.2} {:.2} m {:.2} {:.2} l S\n", b.x, b.y, tip.x, tip.y);
                }
            }
            Graphic::Shape { kind, start, end } => {
                let a = point(*start);
                let b = point(*end);
                match kind {
                    ShapeKind::Line => {
                        stream += &format!("{:.2} {:.2} m {:.2} {:.2} l S\n", a.x, a.y, b.x, b.y)
                    }
                    ShapeKind::Rectangle => {
                        let c = vec2(b.x, a.y);
                        let d = vec2(a.x, b.y);
                        stream += &format!(
                            "{:.2} {:.2} m {:.2} {:.2} l {:.2} {:.2} l {:.2} {:.2} l h S\n",
                            a.x, a.y, c.x, c.y, b.x, b.y, d.x, d.y
                        );
                    }
                    ShapeKind::Ellipse => {
                        let mid = (a + b) / 2.;
                        let radius = (b - a) / 2.;
                        let mut prev = mid + vec2(radius.x, 0.);
                        for i in 1..=48 {
                            let t = i as f32 / 48. * std::f32::consts::TAU;
                            let cur = mid + vec2(radius.x * t.cos(), radius.y * t.sin());
                            stream += &format!(
                                "{:.2} {:.2} m {:.2} {:.2} l S\n",
                                prev.x, prev.y, cur.x, cur.y
                            );
                            prev = cur;
                        }
                    }
                }
            }
            Graphic::Plus { pos } => {
                let p = point(*pos);
                stream += &format!("BT /F1 18 Tf {:.2} {:.2} Td (+) Tj ET\n", p.x, p.y);
            }
            Graphic::Text {
                pos, text, size, ..
            } => {
                let p = point(*pos);
                stream += &format!(
                    "BT /F1 {:.1} Tf {:.2} {:.2} Td ({}) Tj ET\n",
                    size,
                    p.x,
                    p.y,
                    pdf_escape(text)
                );
            }
            Graphic::ElectronArrow { start, end, .. } => {
                let a = point(*start);
                let b = point(*end);
                stream += &format!("{:.2} {:.2} m {:.2} {:.2} l S\n", a.x, a.y, b.x, b.y);
            }
        }
    }
    let objects=["<< /Type /Catalog /Pages 2 0 R >>".to_string(),"<< /Type /Pages /Kids [3 0 R] /Count 1 >>".to_string(),format!("<< /Type /Page /Parent 2 0 R /MediaBox [0 0 {:.0} {:.0}] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>",width,height),"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>".to_string(),format!("<< /Length {} >>\nstream\n{}\nendstream",stream.len(),stream)];
    let mut pdf = String::from("%PDF-1.4\n");
    let mut offsets = vec![0usize];
    for (i, obj) in objects.iter().enumerate() {
        offsets.push(pdf.len());
        pdf += &format!("{} 0 obj\n{}\nendobj\n", i + 1, obj);
    }
    let xref = pdf.len();
    pdf += &format!("xref\n0 {}\n0000000000 65535 f \n", objects.len() + 1);
    for offset in offsets.iter().skip(1) {
        pdf += &format!("{offset:010} 00000 n \n");
    }
    pdf += &format!(
        "trailer\n<< /Size {} /Root 1 0 R >>\nstartxref\n{xref}\n%%EOF\n",
        objects.len() + 1
    );
    pdf.into_bytes()
}
fn pdf_escape(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            '⁺' => '+',
            '⁻' => '-',
            '²' => '2',
            '³' => '3',
            '·' => '.',
            c if c.is_ascii() => c,
            _ => '?',
        })
        .collect::<String>()
        .replace('\\', "\\\\")
        .replace('(', "\\(")
        .replace(')', "\\)")
}
fn charge_label(c: i8) -> String {
    match c {
        0 => String::new(),
        1 => "+".into(),
        -1 => "-".into(),
        x if x > 0 => format!("{}+", x),
        x => format!("{}-", -x),
    }
}
fn charge_suffix(charge: usize) -> &'static str {
    match charge {
        0 | 1 => "",
        2 => "²",
        3 => "³",
        4 => "⁴",
        5 => "⁵",
        6 => "⁶",
        _ => "ⁿ",
    }
}
fn show_label(m: &Molecule, i: usize) -> bool {
    let a = &m.atoms[i];
    a.element != "C"
        || a.alias.is_some()
        || a.charge != 0
        || m.bonds.iter().filter(|b| b.a == i || b.b == i).count() == 0
}
fn atom_display_label(atom: &Atom, implicit_hydrogens: Option<u8>, show_implicit: bool) -> String {
    let mut label = format!(
        "{}{}",
        atom.isotope.map(|n| n.to_string()).unwrap_or_default(),
        atom.alias.as_ref().unwrap_or(&atom.element)
    );
    if show_implicit && atom.alias.is_none() {
        if let Some(count) = implicit_hydrogens.filter(|count| *count > 0) {
            label.push('H');
            if count > 1 {
                label.push_str(&count.to_string());
            }
        }
    }
    label.push_str(&charge_label(atom.charge));
    label
}
fn distance_segment(p: Vec2, a: Vec2, b: Vec2) -> f32 {
    let d = b - a;
    let t = ((p - a).dot(d) / d.length_squared().max(0.001)).clamp(0., 1.);
    p.distance(a + d * t)
}
fn project_3d(point: Vec3, center: Vec3, rotation: Vec2, scale: f32, origin: Vec2) -> Vec2 {
    let p = rotate_3d(point - center, rotation);
    origin + vec2(p.x, -p.y) * scale
}
fn rotate_3d(p: Vec3, rotation: Vec2) -> Vec3 {
    let (sy, cy) = rotation.x.sin_cos();
    let (sx, cx) = rotation.y.sin_cos();
    let x = p.x * cy - p.z * sy;
    let z = p.x * sy + p.z * cy;
    let y = p.y * cx - z * sx;
    let rotated_z = p.y * sx + z * cx;
    vec3(x, y, rotated_z)
}
const THREE_D_STYLE_COUNT: u8 = 5;
fn three_d_style_name(style: u8, language: Language) -> &'static str {
    let names = if language == Language::Danish {
        [
            "Trådmodel",
            "Pinde",
            "Kugle og pind",
            "Rumfyldning",
            "Monokrom",
        ]
    } else {
        [
            "Wireframe",
            "Sticks",
            "Ball & stick",
            "Spacefill",
            "Monochrome",
        ]
    };
    names[usize::from(style.min(THREE_D_STYLE_COUNT - 1))]
}
fn three_d_atom_radius(style: u8, element: &str, scale: f32) -> f32 {
    match style {
        0 => 2.5,
        1 => {
            if element == "C" {
                3.5
            } else {
                7.
            }
        }
        3 => (van_der_waals_radius(element) * scale * 0.55).clamp(7., 90.),
        _ => 13.,
    }
}
fn three_d_atom_color(style: u8, element: &str) -> Color {
    if style == 4 {
        Color::from_rgba(157, 190, 169, 255)
    } else {
        element_color(element)
    }
}
fn three_d_bond_color(style: u8) -> Color {
    match style {
        0 => Color::from_rgba(120, 145, 132, 180),
        4 => Color::from_rgba(98, 142, 116, 220),
        _ => Color::from_rgba(75, 84, 78, 255),
    }
}
fn element_color(element: &str) -> Color {
    match element {
        "H" => Color::from_rgba(255, 255, 255, 255),
        "C" => Color::from_rgba(144, 144, 144, 255),
        "N" => Color::from_rgba(48, 80, 248, 255),
        "O" => Color::from_rgba(255, 13, 13, 255),
        "F" | "Cl" => Color::from_rgba(31, 240, 31, 255),
        "Br" => Color::from_rgba(166, 41, 41, 255),
        "I" => Color::from_rgba(148, 0, 148, 255),
        "S" => Color::from_rgba(255, 200, 50, 255),
        "P" => Color::from_rgba(255, 128, 0, 255),
        "B" => Color::from_rgba(255, 181, 181, 255),
        "Si" => Color::from_rgba(240, 200, 160, 255),
        "Na" | "K" | "Li" | "Rb" | "Cs" => Color::from_rgba(171, 92, 242, 255),
        "Ca" | "Mg" | "Ba" | "Sr" => Color::from_rgba(61, 255, 0, 255),
        "Fe" | "Cu" | "Co" | "Ni" | "Mn" | "Cr" => Color::from_rgba(224, 102, 51, 255),
        "Zn" | "Cd" | "Hg" => Color::from_rgba(125, 128, 176, 255),
        "Al" | "Ga" | "In" | "Sn" | "Pb" => Color::from_rgba(191, 166, 166, 255),
        _ => Color::from_rgba(140, 145, 142, 255),
    }
}
fn van_der_waals_radius(element: &str) -> f32 {
    match element {
        "H" => 1.20,
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "F" => 1.47,
        "P" | "S" => 1.80,
        "Cl" => 1.75,
        "Br" => 1.85,
        "I" => 1.98,
        "He" => 1.40,
        "Ne" => 1.54,
        "Ar" => 1.88,
        "Kr" => 2.02,
        "Xe" => 2.16,
        _ => 1.80,
    }
}
fn v3000_query_attributes(query: Option<&str>) -> (Vec<String>, bool) {
    let mut attributes = BTreeMap::<&'static str, i16>::new();
    let mut unsupported = false;
    for term in query
        .unwrap_or_default()
        .split(';')
        .filter(|term| !term.is_empty())
    {
        let Some(predicate) = shortcuts::parse_query_predicate(term) else {
            unsupported = true;
            continue;
        };
        let value = match predicate {
            shortcuts::QueryPredicate::Substitutions(Some(n)) => {
                Some(("SUBST", if n == 0 { -1 } else { n as i16 }))
            }
            shortcuts::QueryPredicate::Unsaturated => Some(("UNSAT", 1)),
            shortcuts::QueryPredicate::RingBonds(Some(n)) => {
                Some(("RBCNT", if n == 0 { -1 } else { n as i16 }))
            }
            shortcuts::QueryPredicate::AnyAtom
            | shortcuts::QueryPredicate::Heteroatom
            | shortcuts::QueryPredicate::ElementList(_)
            | shortcuts::QueryPredicate::ExcludedElementList(_) => None,
            _ => {
                unsupported = true;
                None
            }
        };
        if let Some((key, value)) = value {
            if attributes.insert(key, value).is_some() {
                unsupported = true;
            }
        }
    }
    (
        attributes
            .into_iter()
            .map(|(key, value)| format!("{key}={value}"))
            .collect(),
        unsupported,
    )
}

/// Encode query element predicates in the V3000 atom-type field. The current
/// graph model stores these as composable query terms, so intersect them into
/// one equivalent finite element set before writing the CTfile list syntax.
fn v3000_atom_type(atom: &Atom) -> (String, bool) {
    let mut allowed: BTreeSet<String> = chemistry::SYMBOLS
        .iter()
        .map(|symbol| (*symbol).to_string())
        .collect();
    let mut has_element_query = false;
    let unsupported = false;
    for term in atom
        .query
        .as_deref()
        .unwrap_or_default()
        .split(';')
        .filter(|t| !t.is_empty())
    {
        let Some(predicate) = shortcuts::parse_query_predicate(term) else {
            // Other query fields are handled by V3000 attributes when possible.
            continue;
        };
        match predicate {
            shortcuts::QueryPredicate::AnyAtom => has_element_query = true,
            shortcuts::QueryPredicate::Heteroatom => {
                has_element_query = true;
                allowed.retain(|element| element != "C" && element != "H");
            }
            shortcuts::QueryPredicate::ElementList(elements) => {
                has_element_query = true;
                allowed.retain(|element| elements.contains(element));
            }
            shortcuts::QueryPredicate::ExcludedElementList(elements) => {
                has_element_query = true;
                allowed.retain(|element| !elements.contains(element));
            }
            _ => {}
        }
    }
    if !has_element_query {
        if !chemistry::valid_element(&atom.element) {
            return (atom.element.clone(), true);
        }
        return (atom.element.clone(), unsupported);
    }
    if allowed.is_empty() {
        return ("A".into(), true);
    }
    let all: BTreeSet<String> = chemistry::SYMBOLS
        .iter()
        .map(|s| (*s).to_string())
        .collect();
    if allowed == all {
        return ("A".into(), unsupported);
    }
    let hetero: BTreeSet<String> = all
        .iter()
        .filter(|element| element.as_str() != "C" && element.as_str() != "H")
        .cloned()
        .collect();
    if allowed == hetero {
        return ("Q".into(), unsupported);
    }
    let excluded: Vec<String> = all.difference(&allowed).cloned().collect();
    if excluded.len() < allowed.len() {
        (format!("NOT [{}]", excluded.join(",")), unsupported)
    } else {
        (
            format!("[{}]", allowed.into_iter().collect::<Vec<_>>().join(",")),
            unsupported,
        )
    }
}

/// Return the CTfile bond type and attributes supported by this document model.
/// V3000 bond types 5–8 are query types; TOPO=1/2 restrict bonds to rings/chains.
fn v3000_bond_attributes(bond: &Bond) -> (u8, Vec<String>, bool) {
    let mut order_type = None;
    let mut topology = None;
    let mut unsupported = false;
    if let Some(query) = bond.query.as_deref() {
        for term in query.split(';').filter(|term| !term.is_empty()) {
            match term {
                "single-or-double" => {
                    if order_type.replace(5).is_some() {
                        unsupported = true;
                    }
                }
                "single-or-aromatic" => {
                    if order_type.replace(6).is_some() {
                        unsupported = true;
                    }
                }
                "double-or-aromatic" => {
                    if order_type.replace(7).is_some() {
                        unsupported = true;
                    }
                }
                "any" => {
                    if order_type.replace(8).is_some() {
                        unsupported = true;
                    }
                }
                "topology-ring" => {
                    if topology.replace(1).is_some() {
                        unsupported = true;
                    }
                }
                "topology-chain" => {
                    if topology.replace(2).is_some() {
                        unsupported = true;
                    }
                }
                _ => unsupported = true,
            }
        }
    }
    let order = order_type.unwrap_or_else(|| {
        if bond.style == BondStyle::Any {
            8
        } else if bond.style == BondStyle::Aromatic || bond.order == 4 {
            4
        } else {
            bond.order
        }
    });
    if !(1..=8).contains(&order) {
        unsupported = true;
    }
    let mut attributes = Vec::new();
    if let Some(topology) = topology {
        attributes.push(format!("TOPO={topology}"));
    }
    let cfg = match bond.style {
        BondStyle::Wedge => Some(1),
        BondStyle::Wavy => Some(2),
        BondStyle::Hash => Some(3),
        _ => None,
    };
    if let Some(cfg) = cfg {
        attributes.push(format!("CFG={cfg}"));
    }
    (order, attributes, unsupported)
}

fn mol_uses_v3000(mol: &Molecule) -> bool {
    mol.atoms.len() > 999
        || mol.bonds.len() > 999
        || mol.atoms.iter().any(|atom| {
            atom.valence_override.is_some_and(|value| value > 14)
                || atom.implicit_h_override == Some(5)
                || atom.query.as_deref().is_some_and(|query| {
                    let (attributes, _) = v3000_query_attributes(Some(query));
                    !attributes.is_empty()
                })
                || v3000_atom_type(atom).0 != atom.element
        })
        || mol.bonds.iter().any(|bond| {
            bond.style == BondStyle::Any
                || bond.query.as_deref().is_some_and(|query| {
                    query.split(';').any(|term| {
                        matches!(
                            term,
                            "any"
                                | "single-or-double"
                                | "single-or-aromatic"
                                | "double-or-aromatic"
                                | "topology-ring"
                                | "topology-chain"
                        )
                    })
                })
        })
}

fn parse_v3000_atom_type(atom_type: &str) -> Result<(String, Option<String>), String> {
    match atom_type {
        "A" | "*" => return Ok(("C".into(), Some("any atom".into()))),
        "Q" => return Ok(("C".into(), Some("heteroatom".into()))),
        _ => {}
    }
    let (excluded, list_text) = if let Some(list) = atom_type.strip_prefix("NOT ") {
        (true, list)
    } else {
        (false, atom_type)
    };
    if let Some(list) = list_text
        .strip_prefix('[')
        .and_then(|s| s.strip_suffix(']'))
    {
        let elements: Vec<&str> = list.split(',').map(str::trim).collect();
        if elements.is_empty()
            || elements
                .iter()
                .any(|element| !chemistry::valid_element(element))
        {
            return Err("V3000 atom list contains an invalid or empty element".into());
        }
        let query = format!("{}{}", if excluded { "!" } else { "" }, elements.join(","));
        if shortcuts::parse_query_predicate(&query).is_none() {
            return Err("V3000 atom list is invalid or contains duplicate elements".into());
        }
        return Ok(("C".into(), Some(query)));
    }
    if excluded {
        return Err("V3000 NOT atom type must contain a bracketed element list".into());
    }
    if chemistry::valid_element(atom_type) {
        Ok((atom_type.into(), None))
    } else {
        Err(format!("Unsupported V3000 atom type: {atom_type}"))
    }
}

fn split_v3000_tokens(record: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let mut current = String::new();
    let mut quoted = false;
    let mut chars = record.chars().peekable();
    while let Some(ch) = chars.next() {
        match ch {
            '"' if quoted && chars.peek() == Some(&'"') => {
                current.push('"');
                chars.next();
            }
            '"' => quoted = !quoted,
            ch if ch.is_whitespace() && !quoted => {
                if !current.is_empty() {
                    tokens.push(std::mem::take(&mut current));
                }
            }
            _ => current.push(ch),
        }
    }
    if !current.is_empty() {
        tokens.push(current);
    }
    tokens
}

fn molfile_for(mol: &Molecule) -> String {
    if mol_uses_v3000(mol) {
        let has_3d = !mol.atoms.is_empty() && mol.atoms.iter().all(|a| a.coord3d.is_some());
        let mut s = format!("{}\n  Molecula\n\n  0  0  0  0  0  0  0  0  0  0999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS {} {} 0 0 0\nM  V30 BEGIN ATOM\n", if has_3d {"MOLECULA_3D"} else {"Molecula"}, mol.atoms.len(), mol.bonds.len());
        for (i, a) in mol.atoms.iter().enumerate() {
            let (atom_type, _) = v3000_atom_type(a);
            let atom_type = if atom_type.chars().any(char::is_whitespace) {
                format!("\"{atom_type}\"")
            } else {
                atom_type
            };
            let (x, y, z) = if has_3d {
                let p = a.coord3d.unwrap();
                (p.x, p.y, p.z)
            } else {
                (a.pos.x, a.pos.y, 0.)
            };
            let mut attrs = String::new();
            if a.charge != 0 {
                attrs += &format!(" CHG={}", a.charge);
            }
            if let Some(n) = a.isotope {
                attrs += &format!(" MASS={n}");
            }
            if a.radical > 0 {
                attrs += &format!(" RAD={}", a.radical);
            }
            if let Some(valence) = a.valence_override {
                attrs += &format!(" VAL={}", if valence == 0 { -1 } else { valence as i16 });
            }
            if let Some(h_count) = a.implicit_h_override.filter(|count| *count <= 5) {
                attrs += &format!(" HCOUNT={}", if h_count == 0 { -1 } else { h_count as i16 });
            }
            let (query_attrs, _) = v3000_query_attributes(a.query.as_deref());
            for query_attr in query_attrs {
                attrs += &format!(" {query_attr}");
            }
            if let Some(cfg) = match a.stereo.as_deref() {
                Some("@") => Some(1),
                Some("@@") => Some(2),
                Some("either") => Some(3),
                _ => None,
            } {
                attrs += &format!(" CFG={cfg}");
            }
            s += &format!(
                "M  V30 {} {} {:.4} {:.4} {:.4} {}{}\n",
                i + 1,
                atom_type,
                x,
                y,
                z,
                a.atom_map.unwrap_or(0),
                attrs
            );
        }
        s += "M  V30 END ATOM\n";
        if !mol.bonds.is_empty() {
            s += "M  V30 BEGIN BOND\n";
        }
        for (i, b) in mol.bonds.iter().enumerate() {
            let (bond_type, attrs, _) = v3000_bond_attributes(b);
            s += &format!(
                "M  V30 {} {} {} {}{}{}\n",
                i + 1,
                bond_type,
                b.a + 1,
                b.b + 1,
                if attrs.is_empty() { "" } else { " " },
                attrs.join(" ")
            );
        }
        if !mol.bonds.is_empty() {
            s += "M  V30 END BOND\n";
        }
        s += "M  V30 END CTAB\nM  END\n";
        return s;
    }
    let has_3d = !mol.atoms.is_empty() && mol.atoms.iter().all(|a| a.coord3d.is_some());
    let mut s = String::from(if has_3d {
        "MOLECULA_3D\n  Rust\n\n"
    } else {
        "Molecular Sketch\n  Rust\n\n"
    });
    s += &format!(
        "{:>3}{:>3}  0  0  0  0            999 V2000\n",
        mol.atoms.len(),
        mol.bonds.len()
    );
    for a in &mol.atoms {
        let (x, y, z) = if has_3d {
            let p = a.coord3d.unwrap();
            (p.x, p.y, p.z)
        } else {
            (
                (a.pos.x - LEFT - 100.) / BOND * 1.5,
                -(a.pos.y - TOP - 100.) / BOND * 1.5,
                0.,
            )
        };
        let h_count = a
            .implicit_h_override
            .filter(|count| *count <= 4 && a.explicit_h == 0)
            .map(|count| count + 1)
            .unwrap_or(0);
        let valence = a
            .valence_override
            .filter(|v| *v <= 14)
            .map(|v| if v == 0 { 15 } else { v });
        // In V2000, an unset atom map is represented by a blank field. Writing
        // zero creates an explicit map number 0 in RDKit and exported SMILES.
        let atom_map = a
            .atom_map
            .map(|map| format!("{map:>3}"))
            .unwrap_or_else(|| "   ".to_string());
        s += &format!(
            "{:>10.4}{:>10.4}{:>10.4} {:<3} 0  0 {:>1} {:>1}  0  0 {:>2}  0  0  {}  0  0\n",
            x,
            y,
            z,
            a.element,
            match a.stereo.as_deref() {
                Some("@") => 1,
                Some("@@") => 2,
                Some("either") => 3,
                _ => 0,
            },
            h_count,
            valence.unwrap_or(0),
            atom_map
        );
    }
    for b in &mol.bonds {
        let stereo = match b.style {
            BondStyle::Wedge => 1,
            BondStyle::Wavy => 4,
            BondStyle::Hash => 6,
            _ => 0,
        };
        s += &format!(
            "{:>3}{:>3}{:>3}{:>3}  0  0  0\n",
            b.a + 1,
            b.b + 1,
            b.order,
            stereo
        );
    }
    for (i, a) in mol.atoms.iter().enumerate() {
        if a.charge != 0 {
            s += &format!("M  CHG  1 {:>3} {:>3}\n", i + 1, a.charge);
        }
        if let Some(isotope) = a.isotope {
            s += &format!("M  ISO  1 {:>3} {:>3}\n", i + 1, isotope);
        }
        if a.radical > 0 {
            s += &format!("M  RAD  1 {:>3} {:>3}\n", i + 1, a.radical);
        }
    }
    s.push_str("M  END\n");
    s
}

fn parse_mol(s: &str) -> Result<Molecule, String> {
    parse_mol_at(
        s,
        vec2(screen_width() / 2., (TOP + screen_height() - BOTTOM) / 2.),
    )
}
fn parse_sdf(input: &str) -> Result<Vec<Molecule>, String> {
    let mut records = Vec::new();
    for raw in input.split("$$$$") {
        if raw.trim().is_empty() {
            continue;
        }
        let mut mol = parse_mol_at(raw, vec2(0., 0.))?;
        let mut ended = false;
        let mut current: Option<String> = None;
        let mut value = String::new();
        for line in raw.lines() {
            if !ended {
                if line.trim() == "M  END" {
                    ended = true;
                }
                continue;
            }
            if let Some(tag) = line
                .trim()
                .strip_prefix(">  <")
                .and_then(|s| s.strip_suffix('>'))
            {
                if let Some(key) = current.replace(tag.to_string()) {
                    mol.properties.push((key, value.trim_end().to_string()));
                }
                value.clear();
            } else if current.is_some() {
                if line.is_empty() {
                    if let Some(key) = current.take() {
                        mol.properties.push((key, value.trim_end().to_string()));
                    }
                    value.clear();
                } else {
                    if !value.is_empty() {
                        value.push('\n');
                    }
                    value.push_str(line);
                }
            }
        }
        if let Some(key) = current {
            mol.properties.push((key, value.trim_end().to_string()));
        }
        records.push(mol);
    }
    Ok(records)
}
fn molecule_components(m: &Molecule) -> Vec<Vec<usize>> {
    let mut seen = vec![false; m.atoms.len()];
    let mut groups = Vec::new();
    for root in 0..m.atoms.len() {
        if seen[root] {
            continue;
        }
        let mut stack = vec![root];
        let mut group = Vec::new();
        seen[root] = true;
        while let Some(i) = stack.pop() {
            group.push(i);
            for b in &m.bonds {
                let next = if b.a == i {
                    Some(b.b)
                } else if b.b == i {
                    Some(b.a)
                } else {
                    None
                };
                if let Some(j) = next {
                    if !seen[j] {
                        seen[j] = true;
                        stack.push(j);
                    }
                }
            }
        }
        groups.push(group);
    }
    groups
}
fn url_encode(value: &str) -> String {
    let mut out = String::new();
    for byte in value.bytes() {
        if byte.is_ascii_alphanumeric() || b"-_.~".contains(&byte) {
            out.push(byte as char)
        } else {
            out += &format!("%{byte:02X}");
        }
    }
    out
}
fn open_external_url(url: &str) -> Result<(), String> {
    #[cfg(target_os = "macos")]
    let result = std::process::Command::new("open").arg(url).status();
    #[cfg(target_os = "windows")]
    let result = std::process::Command::new("rundll32")
        .arg("url.dll,FileProtocolHandler")
        .arg(url)
        .status();
    #[cfg(all(unix, not(target_os = "macos")))]
    let result = std::process::Command::new("xdg-open").arg(url).status();
    result.map_err(|e| format!("Could not launch a browser: {e}"))?;
    Ok(())
}
fn rxn_for(m: &Molecule) -> Result<String, String> {
    if m.atoms.is_empty() {
        return Err("RXN export needs at least one atom".into());
    }
    let mean_x = m.atoms.iter().map(|a| a.pos.x).sum::<f32>() / m.atoms.len() as f32;
    let mut sides: [Vec<Molecule>; 3] = std::array::from_fn(|_| Vec::new());
    for group in molecule_components(m) {
        let explicit: BTreeSet<_> = group
            .iter()
            .filter_map(|&i| m.atoms[i].reaction_role.as_deref())
            .collect();
        if explicit.len() > 1 {
            return Err("A connected component has conflicting reaction roles".into());
        }
        let center = group.iter().map(|&i| m.atoms[i].pos.x).sum::<f32>() / group.len() as f32;
        let role = match explicit
            .iter()
            .next()
            .copied()
            .unwrap_or(if center < mean_x {
                "reactant"
            } else {
                "product"
            }) {
            "reactant" => 0,
            "product" => 1,
            "agent" => 2,
            _ => return Err("Unknown reaction component role".into()),
        };
        let mut sub = Molecule::default();
        let mut map = BTreeMap::new();
        for &i in &group {
            map.insert(i, sub.atoms.len());
            sub.atoms.push(m.atoms[i].clone());
        }
        for b in &m.bonds {
            if let (Some(&a), Some(&c)) = (map.get(&b.a), map.get(&b.b)) {
                let mut copy = b.clone();
                copy.a = a;
                copy.b = c;
                sub.bonds.push(copy);
            }
        }
        sides[role].push(sub);
    }
    if sides[0].is_empty() || sides[1].is_empty() {
        return Err("RXN export needs at least one reactant and one product; assign component roles in the inspector".into());
    }
    let mut out = format!(
        "$RXN\n{}\n  Molecula\n\n{:>3}{:>3}{:>3}\n",
        m.name.replace(['\n', '\r'], " "),
        sides[0].len(),
        sides[1].len(),
        sides[2].len()
    );
    for side in sides {
        for mol in side {
            out += "$MOL\n";
            out += &molfile_for(&mol);
        }
    }
    Ok(out)
}
fn parse_rxn(input: &str, center: Vec2) -> Result<Molecule, String> {
    let lines: Vec<_> = input.lines().collect();
    if lines.first() != Some(&"$RXN") {
        return Err("Invalid RXN header".into());
    }
    let counts_i = lines
        .iter()
        .position(|line| line.trim().starts_with("$RXN"))
        .ok_or("Missing RXN header")?
        + 4;
    let counts: Vec<_> = lines
        .get(counts_i)
        .ok_or("RXN counts are missing")?
        .split_whitespace()
        .collect();
    if counts.len() < 2 {
        return Err("Invalid RXN counts".into());
    }
    let nr = counts[0]
        .parse::<usize>()
        .map_err(|_| "Invalid reactant count")?;
    let np = counts[1]
        .parse::<usize>()
        .map_err(|_| "Invalid product count")?;
    let ng = counts.get(2).map_or(Ok(0), |x| {
        x.parse::<usize>().map_err(|_| "Invalid agent count")
    })?;
    if nr + np + ng > 10000 {
        return Err("RXN contains too many components".into());
    }
    let blocks: Vec<_> = input.split("$MOL").skip(1).collect();
    if blocks.len() != nr + np + ng {
        return Err("RXN component count does not match its header".into());
    }
    let mut out = Molecule {
        name: lines.get(1).copied().unwrap_or("").to_string(),
        ..Molecule::default()
    };
    for (ci, block) in blocks.iter().enumerate() {
        let role = if ci < nr {
            "reactant"
        } else if ci < nr + np {
            "product"
        } else {
            "agent"
        };
        let offset = vec2(center.x + (ci as f32) * BOND * 2., center.y);
        let sub = parse_mol_at(block, offset)?;
        let base = out.atoms.len();
        for mut a in sub.atoms {
            a.reaction_role = Some(role.into());
            out.atoms.push(a);
        }
        for mut b in sub.bonds {
            b.a += base;
            b.b += base;
            out.bonds.push(b);
        }
    }
    Ok(out)
}
fn parse_mol_at(s: &str, center: Vec2) -> Result<Molecule, String> {
    let lines: Vec<&str> = s.lines().collect();
    if lines.iter().any(|l| l.contains("V3000")) {
        let counts = lines
            .iter()
            .find(|l| l.contains("V30 COUNTS"))
            .ok_or("V3000 file is missing counts")?;
        let counts: Vec<_> = counts.split_whitespace().collect();
        if counts.len() < 7 {
            return Err("Invalid V3000 counts line".into());
        }
        let na = counts[3]
            .parse::<usize>()
            .map_err(|_| "Invalid V3000 atom count")?;
        let nb = counts[4]
            .parse::<usize>()
            .map_err(|_| "Invalid V3000 bond count")?;
        if na > 100_000 || nb > 200_000 {
            return Err("V3000 structure exceeds safety limits".into());
        }
        let is_3d = lines.first().is_some_and(|l| l.trim() == "MOLECULA_3D");
        let mut m = Molecule::default();
        let mut coords = Vec::new();
        let mut in_atoms = false;
        let mut in_bonds = false;
        let mut ended = false;
        for line in &lines {
            let line = line.trim();
            if line == "M  V30 BEGIN ATOM" {
                in_atoms = true;
                continue;
            }
            if line == "M  V30 END ATOM" {
                in_atoms = false;
                continue;
            }
            if line == "M  V30 BEGIN BOND" {
                in_bonds = true;
                continue;
            }
            if line == "M  V30 END BOND" {
                in_bonds = false;
                continue;
            }
            if line == "M  END" {
                ended = true;
                break;
            }
            if in_atoms && line.starts_with("M  V30 ") {
                let f = split_v3000_tokens(&line[7..]);
                if f.len() < 6 {
                    return Err("Invalid V3000 atom record".into());
                }
                let index = f[0]
                    .parse::<usize>()
                    .map_err(|_| "Invalid V3000 atom index")?;
                if index != m.atoms.len() + 1 {
                    return Err("V3000 atom indices must be sequential".into());
                }
                let not_list = f[1] == "NOT";
                let type_offset = usize::from(not_list);
                if f.len() < 6 + type_offset {
                    return Err("Invalid V3000 atom-list record".into());
                }
                let raw_type = if not_list {
                    if !f[2].starts_with('[') {
                        return Err(
                            "V3000 NOT atom type must be followed by a bracketed list".into()
                        );
                    }
                    format!("NOT {}", f[2])
                } else {
                    f[1].clone()
                };
                let (element, element_query) = parse_v3000_atom_type(&raw_type)?;
                let x = f[2 + type_offset]
                    .parse::<f32>()
                    .map_err(|_| "Invalid coordinate")?;
                let y = f[3 + type_offset]
                    .parse::<f32>()
                    .map_err(|_| "Invalid coordinate")?;
                let z = f[4 + type_offset]
                    .parse::<f32>()
                    .map_err(|_| "Invalid coordinate")?;
                if !x.is_finite() || !y.is_finite() || !z.is_finite() {
                    return Err("MOL coordinates must be finite numbers".into());
                }
                coords.push(vec2(x, y));
                let mut a = Atom::new(vec2(0., 0.), &element);
                a.query = element_query;
                if is_3d {
                    a.coord3d = Some(vec3(x, y, z));
                }
                if let Ok(map) = f[5 + type_offset].parse::<u16>() {
                    if map > 0 {
                        a.atom_map = Some(map);
                    }
                }
                for attr in f.iter().skip(6 + type_offset) {
                    if let Some(v) = attr.strip_prefix("CHG=") {
                        a.charge = v.parse().map_err(|_| "Invalid V3000 charge")?;
                    } else if let Some(v) = attr.strip_prefix("MASS=") {
                        a.isotope = Some(v.parse().map_err(|_| "Invalid V3000 isotope")?);
                    } else if let Some(v) = attr.strip_prefix("RAD=") {
                        a.radical = v.parse().map_err(|_| "Invalid V3000 radical")?;
                        if !(1..=3).contains(&a.radical) {
                            return Err(
                                "V3000 radical must be singlet (1), doublet (2), or triplet (3)"
                                    .into(),
                            );
                        }
                    } else if let Some(v) = attr.strip_prefix("VAL=") {
                        let valence = v.parse::<i16>().map_err(|_| "Invalid V3000 valence")?;
                        if !(-1..=255).contains(&valence) {
                            return Err(
                                "V3000 valence is outside the supported integer range".into()
                            );
                        }
                        if valence != 0 {
                            a.valence_override =
                                Some(if valence == -1 { 0 } else { valence as u8 });
                        }
                    } else if let Some(v) = attr.strip_prefix("HCOUNT=") {
                        let count = v.parse::<i16>().map_err(|_| "Invalid V3000 HCOUNT")?;
                        if !(-1..=5).contains(&count) {
                            return Err("V3000 HCOUNT must be -1, 0, or 1 through 5".into());
                        }
                        a.implicit_h_override = match count {
                            -1 => Some(0),
                            0 => None,
                            value => Some(value as u8),
                        };
                    } else if let Some(v) = attr.strip_prefix("SUBST=") {
                        let count = v.parse::<i16>().map_err(|_| "Invalid V3000 SUBST")?;
                        if !(-1..=255).contains(&count) {
                            return Err("V3000 SUBST must be -1, 0, or 1 through 255".into());
                        }
                        let term = match count {
                            -1 => Some(".s0".to_string()),
                            0 => None,
                            value => Some(format!(".s{value}")),
                        };
                        if let Some(term) = term {
                            a.query =
                                Some(shortcuts::merge_query_terms(a.query.as_deref(), &term)?);
                        }
                    } else if let Some(v) = attr.strip_prefix("UNSAT=") {
                        let value = v.parse::<u8>().map_err(|_| "Invalid V3000 UNSAT")?;
                        if value > 1 {
                            return Err("V3000 UNSAT must be 0 or 1".into());
                        }
                        if value == 1 {
                            a.query = Some(shortcuts::merge_query_terms(a.query.as_deref(), ".u")?);
                        }
                    } else if let Some(v) = attr.strip_prefix("RBCNT=") {
                        let count = v.parse::<i16>().map_err(|_| "Invalid V3000 RBCNT")?;
                        if !(-1..=255).contains(&count) {
                            return Err("V3000 RBCNT must be -1, 0, or 1 through 255".into());
                        }
                        let term = match count {
                            -1 => Some(".rb0".to_string()),
                            0 => None,
                            value => Some(format!(".rb{value}")),
                        };
                        if let Some(term) = term {
                            a.query =
                                Some(shortcuts::merge_query_terms(a.query.as_deref(), &term)?);
                        }
                    } else if let Some(v) = attr.strip_prefix("CFG=") {
                        a.stereo = match v {
                            "1" => Some("@".into()),
                            "2" => Some("@@".into()),
                            "3" => Some("either".into()),
                            _ => None,
                        };
                    }
                }
                m.atoms.push(a);
            }
            if in_bonds && line.starts_with("M  V30 ") {
                let f: Vec<_> = line[7..].split_whitespace().collect();
                if f.len() < 4 {
                    return Err("Invalid V3000 bond record".into());
                }
                let v3000_type = f[1].parse::<u8>().map_err(|_| "Invalid V3000 bond order")?;
                let a = f[2]
                    .parse::<usize>()
                    .map_err(|_| "Invalid V3000 bond atom")?;
                let b = f[3]
                    .parse::<usize>()
                    .map_err(|_| "Invalid V3000 bond atom")?;
                if a == 0 || b == 0 || a > na || b > na || a == b || !(1..=8).contains(&v3000_type)
                {
                    return Err("V3000 bond references invalid atoms or order".into());
                }
                let mut bond_query = match v3000_type {
                    5 => Some("single-or-double".to_string()),
                    6 => Some("single-or-aromatic".to_string()),
                    7 => Some("double-or-aromatic".to_string()),
                    8 => Some("any".to_string()),
                    _ => None,
                };
                let mut topology_seen = None;
                for attr in f.iter().skip(4) {
                    if let Some(value) = attr.strip_prefix("TOPO=") {
                        match value {
                            "0" | "1" | "2" => {
                                if topology_seen.is_some_and(|seen| seen != value) {
                                    return Err("Conflicting V3000 TOPO values on one bond".into());
                                }
                                topology_seen = Some(value);
                                if value == "0" {
                                    continue;
                                }
                                let query_term = if value == "1" {
                                    "topology-ring"
                                } else {
                                    "topology-chain"
                                };
                                if let Some(query) = bond_query.as_mut() {
                                    if !query.split(';').any(|term| term == query_term) {
                                        query.push(';');
                                        query.push_str(query_term);
                                    }
                                } else {
                                    bond_query = Some(query_term.into());
                                }
                            }
                            _ => {
                                return Err(
                                    "Unsupported V3000 bond topology; expected TOPO=0, TOPO=1 or TOPO=2"
                                        .into(),
                                )
                            }
                        }
                    }
                }
                let style = f
                    .iter()
                    .skip(4)
                    .find_map(|x| x.strip_prefix("CFG="))
                    .and_then(|v| v.parse::<u8>().ok())
                    .map(|cfg| {
                        if cfg == 1 {
                            BondStyle::Wedge
                        } else if cfg == 2 {
                            BondStyle::Wavy
                        } else if cfg == 3 {
                            BondStyle::Hash
                        } else {
                            BondStyle::Normal
                        }
                    })
                    .unwrap_or(if v3000_type == 4 {
                        BondStyle::Aromatic
                    } else {
                        BondStyle::Normal
                    });
                m.bonds.push(Bond {
                    a: a - 1,
                    b: b - 1,
                    order: if v3000_type >= 5 { 1 } else { v3000_type },
                    style,
                    stereo: None,
                    query: bond_query,
                });
            }
        }
        if !ended || m.atoms.len() != na || m.bonds.len() != nb {
            return Err("V3000 file is incomplete or record counts do not match".into());
        }
        for bond in m.bonds.iter().filter(|bond| bond.order == 4) {
            m.atoms[bond.a].aromatic = true;
            m.atoms[bond.b].aromatic = true;
        }
        let mid = if coords.is_empty() {
            vec2(0., 0.)
        } else {
            coords.iter().copied().sum::<Vec2>() / coords.len() as f32
        };
        for (atom, p) in m.atoms.iter_mut().zip(coords) {
            atom.pos = center + vec2(p.x - mid.x, -p.y + mid.y) * (BOND / 1.5);
        }
        return Ok(m);
    }
    let ci = lines
        .iter()
        .position(|l| l.contains("V2000"))
        .ok_or("Only MOL V2000 files are supported")?;
    let count = lines[ci];
    let fields: Vec<&str> = count.split_whitespace().collect();
    if fields.len() < 2 {
        return Err("Invalid MOL counts line".into());
    }
    let na = fields[0]
        .parse::<usize>()
        .map_err(|_| "Invalid atom count")?;
    let nb = fields[1]
        .parse::<usize>()
        .map_err(|_| "Invalid bond count")?;
    if na > 999 || nb > 999 || lines.len() < ci + 1 + na + nb {
        return Err("MOL file is incomplete or too large".into());
    }
    if !lines
        .iter()
        .skip(ci + 1 + na + nb)
        .any(|l| l.trim() == "M  END")
    {
        return Err("MOL file is missing M  END".into());
    }
    let mut m = Molecule::default();
    let mut coords = Vec::new();
    let is_3d = lines
        .first()
        .is_some_and(|line| line.trim() == "MOLECULA_3D");
    for i in 0..na {
        let f: Vec<&str> = lines[ci + 1 + i].split_whitespace().collect();
        if f.len() < 4 {
            return Err("Invalid atom line".into());
        }
        let x = f[0].parse::<f32>().map_err(|_| "Invalid coordinate")?;
        let y = f[1].parse::<f32>().map_err(|_| "Invalid coordinate")?;
        if !x.is_finite() || !y.is_finite() {
            return Err("MOL coordinates must be finite numbers".into());
        }
        let z = f
            .get(2)
            .and_then(|field| field.parse::<f32>().ok())
            .filter(|z| z.is_finite())
            .unwrap_or(0.);
        coords.push(vec2(x, y));
        let mut atom = Atom::new(vec2(0., 0.), f[3]);
        if let Some(map) = f.get(13).and_then(|value| value.parse::<u16>().ok()) {
            if map > 0 {
                atom.atom_map = Some(map);
            }
        }
        if let Some(parity) = f.get(6).and_then(|value| value.parse::<u8>().ok()) {
            atom.stereo = match parity {
                1 => Some("@".into()),
                2 => Some("@@".into()),
                3 => Some("either".into()),
                _ => None,
            };
        }
        if let Some(h_count) = f.get(7).and_then(|value| value.parse::<u8>().ok()) {
            if h_count > 5 {
                return Err("Invalid MOL hydrogen count (expected 0 through 5)".into());
            }
            if h_count > 0 {
                atom.implicit_h_override = Some(h_count - 1);
            }
        }
        if let Some(valence) = f.get(10).and_then(|value| value.parse::<u8>().ok()) {
            if valence > 15 {
                return Err("Invalid MOL valence (expected 0 through 15)".into());
            }
            if valence != 0 {
                atom.valence_override = Some(if valence == 15 { 0 } else { valence });
            }
        }
        if is_3d {
            atom.coord3d = Some(vec3(x, y, z));
        }
        m.atoms.push(atom);
    }
    let mid = if coords.is_empty() {
        vec2(0., 0.)
    } else {
        coords.iter().copied().sum::<Vec2>() / na as f32
    };
    for (a, p) in m.atoms.iter_mut().zip(coords) {
        a.pos = center + vec2(p.x - mid.x, -p.y + mid.y) * (BOND / 1.5);
    }
    for i in 0..nb {
        let f: Vec<&str> = lines[ci + 1 + na + i].split_whitespace().collect();
        if f.len() < 3 {
            return Err("Invalid bond line".into());
        }
        let a = f[0].parse::<usize>().map_err(|_| "Invalid bond")?;
        let b = f[1].parse::<usize>().map_err(|_| "Invalid bond")?;
        let order = f[2].parse::<u8>().map_err(|_| "Invalid bond order")?;
        if a == 0 || b == 0 || a > na || b > na || a == b {
            return Err("Bond references missing atom".into());
        }
        if !(1..=4).contains(&order) {
            return Err("Unsupported MOL bond order".into());
        }
        let stereo = f
            .get(3)
            .and_then(|value| value.parse::<u8>().ok())
            .unwrap_or(0);
        m.bonds.push(Bond {
            a: a - 1,
            b: b - 1,
            order,
            style: match stereo {
                1 => BondStyle::Wedge,
                4 => BondStyle::Wavy,
                6 => BondStyle::Hash,
                _ if order == 4 => BondStyle::Aromatic,
                _ => BondStyle::Normal,
            },
            stereo: None,
            query: None,
        });
    }
    for l in lines.iter().skip(ci + 1 + na + nb) {
        let f: Vec<&str> = l.split_whitespace().collect();
        if f.len() >= 5 && f[0] == "M" && f[1] == "CHG" {
            if let (Ok(i), Ok(c)) = (f[3].parse::<usize>(), f[4].parse::<i8>()) {
                if i > 0 && i <= na {
                    m.atoms[i - 1].charge = c;
                }
            }
        }
        if f.len() >= 5 && f[0] == "M" && f[1] == "ISO" {
            if let (Ok(i), Ok(n)) = (f[3].parse::<usize>(), f[4].parse::<u16>()) {
                if i > 0 && i <= na {
                    m.atoms[i - 1].isotope = Some(n);
                }
            }
        }
        if f.len() >= 3 && f[0] == "M" && f[1] == "RAD" {
            let pair_count = f[2].parse::<usize>().map_err(|_| "Invalid radical count")?;
            if pair_count > na || f.len() < 3 + pair_count * 2 {
                return Err("MOL radical record is incomplete".into());
            }
            for pair in 0..pair_count {
                let atom_index = f[3 + pair * 2]
                    .parse::<usize>()
                    .map_err(|_| "Invalid radical atom index")?;
                let radical = f[4 + pair * 2]
                    .parse::<u8>()
                    .map_err(|_| "Invalid radical type")?;
                if atom_index == 0 || atom_index > na || !(1..=3).contains(&radical) {
                    return Err("MOL radical must reference an atom and use type 1, 2, or 3".into());
                }
                m.atoms[atom_index - 1].radical = radical;
            }
        }
    }
    for bond in m.bonds.iter().filter(|bond| bond.order == 4) {
        m.atoms[bond.a].aromatic = true;
        m.atoms[bond.b].aromatic = true;
    }
    Ok(m)
}

fn button(x: f32, y: f32, w: f32, h: f32, label: &str, active: bool) -> bool {
    let label = UI_LANGUAGE.with(|language| localized(language.get(), label));
    let p = vec2(mouse_position().0, mouse_position().1);
    let hover = p.x >= x && p.x < x + w && p.y >= y && p.y < y + h;
    let top = y < 40.;
    let toolbar = (40. ..TOP).contains(&y);
    let palette = y >= screen_height() - BOTTOM && y < screen_height() - BOTTOM + 52.;
    let fill = if top {
        if active {
            Color::from_rgba(45, 102, 76, 255)
        } else if hover {
            Color::from_rgba(32, 74, 56, 255)
        } else {
            navy()
        }
    } else if active {
        if palette {
            accent()
        } else {
            Color::from_rgba(229, 237, 255, 255)
        }
    } else if hover {
        Color::from_rgba(235, 240, 250, 255)
    } else if toolbar {
        Color::from_rgba(255, 255, 255, 18)
    } else {
        WHITE
    };
    round_rect(x, y, w, h, 8., fill);
    if !top && !toolbar && !active && !hover {
        round_border(x, y, w, h, 8., Color::from_rgba(220, 228, 239, 255));
    }
    let mut size = if top {
        14.
    } else if toolbar {
        15.
    } else {
        16.
    };
    while size > 11. && measure_text(label.as_ref(), None, size as u16, 1.).width > w - 8. {
        size -= 1.;
    }
    let d = measure_text(label.as_ref(), None, size as u16, 1.);
    draw_text(
        label.as_ref(),
        x + (w - d.width) / 2.,
        y + (h + d.height) / 2. - 3.,
        size,
        if top || palette && active {
            WHITE
        } else if active {
            accent()
        } else {
            ink()
        },
    );
    hover && is_mouse_button_pressed(MouseButton::Left) && UI_INPUT_ENABLED.with(Cell::get)
}
fn draw_bond(a: Vec2, b: Vec2, order: u8, color: Color) {
    let d = (b - a).normalize_or_zero();
    let q = vec2(-d.y, d.x);
    for k in 0..order {
        let o = (k as f32 - (order as f32 - 1.) / 2.) * 5.;
        let start = a + q * o;
        let end = b + q * o;
        draw_line(start.x, start.y, end.x, end.y, 2., color);
    }
}
fn draw_chem_bond(a: Vec2, b: Vec2, bond: &Bond, color: Color) {
    let d = (b - a).normalize_or_zero();
    let q = vec2(-d.y, d.x);
    match bond.style {
        BondStyle::Normal => draw_bond(a, b, bond.order.min(3), color),
        BondStyle::Aromatic => {
            draw_bond(a, b, 1, color);
            for k in 0..6 {
                let t = k as f32 / 6.;
                let u = (t + 0.08).min(1.);
                let p = a.lerp(b, t) + q * 5.;
                let r = a.lerp(b, u) + q * 5.;
                draw_line(p.x, p.y, r.x, r.y, 1.5, color);
            }
        }
        BondStyle::Wedge => draw_triangle(a, b + q * 6., b - q * 6., color),
        BondStyle::Hash => {
            for k in 0..11 {
                let t = k as f32 / 11.;
                let p = a.lerp(b, t);
                let width = 1. + t * 6.;
                draw_line(
                    p.x - q.x * width,
                    p.y - q.y * width,
                    p.x + q.x * width,
                    p.y + q.y * width,
                    1.5,
                    color,
                );
            }
        }
        BondStyle::Wavy => {
            let mut prev = a;
            for k in 1..=24 {
                let t = k as f32 / 24.;
                let p = a.lerp(b, t) + q * (t * std::f32::consts::TAU * 4.).sin() * 3.;
                draw_line(prev.x, prev.y, p.x, p.y, 1.7, color);
                prev = p;
            }
        }
        BondStyle::Coordinate => {
            draw_line(a.x, a.y, b.x, b.y, 1.7, color);
            draw_triangle(b, b - d * 10. + q * 4., b - d * 10. - q * 4., color);
        }
        BondStyle::Any => {
            for k in 0..8 {
                let p = a.lerp(b, k as f32 / 8.);
                let r = a.lerp(b, (k as f32 + 0.55) / 8.);
                draw_line(p.x, p.y, r.x, r.y, 2., color);
            }
        }
    }
}
fn canvas(p: Vec2) -> bool {
    p.x > LEFT && p.x < screen_width() - RIGHT && p.y > TOP && p.y < screen_height() - BOTTOM
}
fn electron_control(start: Vec2, end: Vec2) -> Vec2 {
    let d = end - start;
    (start + end) / 2. + vec2(d.y, -d.x) * 0.34
}
fn draw_electron_arrow(start: Vec2, control: Vec2, end: Vec2, electrons: u8, color: Color) {
    let mut prev = start;
    for k in 1..=32 {
        let t = k as f32 / 32.;
        let p = start * (1. - t).powi(2) + control * (2. * t * (1. - t)) + end * t * t;
        draw_line(prev.x, prev.y, p.x, p.y, 2., color);
        prev = p;
    }
    let d = (end - control).normalize_or_zero();
    let q = vec2(-d.y, d.x);
    draw_line(
        end.x,
        end.y,
        end.x - d.x * 12. + q.x * 5.,
        end.y - d.y * 12. + q.y * 5.,
        2.,
        color,
    );
    if electrons == 2 {
        draw_line(
            end.x,
            end.y,
            end.x - d.x * 12. - q.x * 5.,
            end.y - d.y * 12. - q.y * 5.,
            2.,
            color,
        );
    }
}
fn draw_reaction_arrow(start: Vec2, end: Vec2, kind: &str, color: Color) {
    let d = (end - start).normalize_or_zero();
    let q = vec2(-d.y, d.x);
    if kind == "equilibrium" {
        let a = start - q * 4.;
        let b = end - q * 4.;
        draw_line(a.x, a.y, b.x, b.y, 2., color);
        draw_line(
            b.x,
            b.y,
            b.x - d.x * 10. - q.x * 5.,
            b.y - d.y * 10. - q.y * 5.,
            2.,
            color,
        );
        let a = end + q * 4.;
        let b = start + q * 4.;
        draw_line(a.x, a.y, b.x, b.y, 2., color);
        draw_line(
            b.x,
            b.y,
            b.x + d.x * 10. + q.x * 5.,
            b.y + d.y * 10. + q.y * 5.,
            2.,
            color,
        );
    } else {
        draw_line(start.x, start.y, end.x, end.y, 2., color);
        let head = |at: Vec2, dir: Vec2| {
            draw_line(
                at.x,
                at.y,
                at.x - dir.x * 12. + q.x * 6.,
                at.y - dir.y * 12. + q.y * 6.,
                2.,
                color,
            );
            draw_line(
                at.x,
                at.y,
                at.x - dir.x * 12. - q.x * 6.,
                at.y - dir.y * 12. - q.y * 6.,
                2.,
                color,
            );
        };
        head(end, d);
        if kind == "resonance" {
            head(start, -d);
        }
        if kind == "retrosynthesis" {
            draw_line(
                end.x - q.x * 8.,
                end.y - q.y * 8.,
                end.x + q.x * 8.,
                end.y + q.y * 8.,
                2.,
                color,
            );
        }
    }
}
fn periodic_cell(z: usize) -> (usize, usize) {
    let n = z + 1;
    match n {
        1 => (0, 0),
        2 => (17, 0),
        3..=10 => {
            let j = n - 3;
            (if j < 2 { j } else { j + 10 }, 1)
        }
        11..=18 => {
            let j = n - 11;
            (if j < 2 { j } else { j + 10 }, 2)
        }
        19..=36 => (n - 19, 3),
        37..=54 => (n - 37, 4),
        55..=56 => (n - 55, 5),
        57..=71 => (n - 57 + 2, 7),
        72..=86 => (n - 72 + 3, 5),
        87..=88 => (n - 87, 6),
        89..=103 => (n - 89 + 2, 8),
        104..=118 => (n - 104 + 3, 6),
        _ => (0, 0),
    }
}
const MENU_NAMES: [&str; 10] = [
    "File",
    "Edit",
    "View",
    "Insert",
    "Atom",
    "Bond",
    "Structure",
    "Calculations",
    "Services",
    "Help",
];
fn menu_items(i: usize) -> &'static [&'static str] {
    match i {
        0 => &[
            "New",
            "Open",
            "Save",
            "Save all",
            "Recent files",
            "Close tab",
        ],
        1 => &[
            "Undo",
            "Redo",
            "Copy as SMILES",
            "Copy as SMARTS",
            "Copy drawing as image",
            "Delete selection",
            "Duplicate selection",
            "Rotate selection 90°",
            "Select all",
            "Settings",
        ],
        2 => &[
            "Fit",
            "Zoom in",
            "Zoom out",
            "Elements",
            "Toggle dark theme",
            "Show implicit hydrogens",
        ],
        3 => &[
            "Templates",
            "Save selection as template",
            "Arrow",
            "Equilibrium arrow",
            "Resonance arrow",
            "Retrosynthesis arrow",
            "Electron pair arrow",
            "Single electron arrow",
            "Plus",
            "Text",
            "Edit molecule data",
            "Line",
            "Rectangle",
            "Ellipse",
        ],
        4 => &["Elements", "Charge +", "Charge -"],
        5 => &[
            "Single", "Double", "Triple", "Aromatic", "Wedge", "Hash", "Wavy",
        ],
        6 => &[
            "Validate",
            "Find substructure matches",
            "Name structure",
            "Auto-map reaction (Complete)",
            "Auto-map reaction (Matching)",
            "Auto-map reaction (Changing)",
            "Clean 2D",
            "Aromatize",
            "Kekulé form",
            "Stereochemistry report",
            "Assign R/S labels with RDKit",
            "Generate 3D coordinates",
            "3D viewer",
            "Add explicit hydrogens",
            "Remove explicit hydrogens",
            "Rotate 90 CW",
            "Rotate 90 CCW",
            "Rotate 180",
            "Flip left-right",
            "Flip up-down",
            "Benzene template",
        ],
        7 => &[
            "Formula & mass",
            "Molecular properties",
            "Local topological descriptors",
            "Chemistry A acid/base tools",
            "Reaction analysis",
            "Auto-balance reaction",
            "Auto-balance all reaction steps",
            "View all reaction steps",
            "View previous reaction step",
            "View next reaction step",
            "View reaction step 1",
            "View reaction step 2",
            "View reaction step 3",
            "Reaction component: reactant",
            "Reaction component: product",
            "Reaction component: agent",
            "Reaction component: automatic role",
            "Reaction coefficient: 2",
            "Reaction coefficient: 3",
            "Reaction coefficient: custom",
            "Reaction step: 1",
            "Reaction step: 2",
            "Reaction step: 3",
            "RDKit descriptors",
            "Gasteiger partial charges",
            "Enumerate tautomers",
            "Enumerate stereoisomers",
        ],
        8 => &[
            "Fetch PubChem properties",
            "Search PubChem in browser",
            "Integration status",
        ],
        9 => &[
            "Quick start guide",
            "Naming lesson",
            "Reaction balancing lesson",
            "Resume lesson",
            "Shortcuts",
            "Feature status",
            "Command palette",
        ],
        _ => &[],
    }
}
fn all_commands() -> Vec<&'static str> {
    let mut commands = Vec::new();
    for menu in 0..MENU_NAMES.len() {
        for command in menu_items(menu) {
            if !commands.contains(command) {
                commands.push(*command);
            }
        }
    }
    commands
}

fn wrap_info_lines(text: &str, width: usize) -> Vec<String> {
    let width = width.max(1);
    let mut wrapped = Vec::new();
    for line in text.lines() {
        let chars: Vec<char> = line.chars().collect();
        if chars.is_empty() {
            wrapped.push(String::new());
        } else {
            wrapped.extend(chars.chunks(width).map(|part| part.iter().collect()));
        }
    }
    wrapped
}

fn background_calculation_error(message: &str) -> (String, String) {
    if let Some(detail) = message.strip_prefix("PUBCHEM_LOOKUP\n") {
        (
            format!("PubChem lookup failed. The drawing was not changed.\n\n{detail}"),
            "PubChem lookup failed".into(),
        )
    } else {
        (
            format!("RDKit unavailable: {message}\nInstall RDKit in a Python environment and set MOLECULA_PYTHON to its interpreter."),
            "Calculation failed".into(),
        )
    }
}

#[macroquad::main(window_conf)]
async fn main() {
    if let Ok(font) = load_ttf_font_from_bytes(include_bytes!("../assets/DejaVuSans.ttf")) {
        UI_FONT.with(|f| *f.borrow_mut() = Some(font));
    }
    let mut e = Editor::new();
    let mut autosave_clock = 0.0f32;
    #[cfg(target_os = "windows")]
    record_startup_note("Editor initialized; ready.");
    UI_LANGUAGE.with(|language| language.set(e.language));
    UI_DARK_MODE.with(|mode| mode.set(e.dark_theme));
    loop {
        UI_LANGUAGE.with(|language| language.set(e.language));
        UI_DARK_MODE.with(|mode| mode.set(e.dark_theme));
        UI_INPUT_ENABLED.with(|enabled| enabled.set(!e.settings_open));
        e.persist_zoom_preference();
        if is_quit_requested() {
            prevent_quit();
            if e.documents.iter().any(|doc| doc.dirty) {
                e.pending_quit = true;
            } else {
                e.quit_now = true;
            }
        }
        if e.quit_when_clean && !e.documents.iter().any(|doc| doc.dirty) {
            e.quit_now = true;
        }
        if e.quit_now {
            break;
        }
        let quit_shortcut = is_key_down(KeyCode::LeftSuper)
            || is_key_down(KeyCode::RightSuper)
            || is_key_down(KeyCode::LeftControl)
            || is_key_down(KeyCode::RightControl);
        if quit_shortcut && is_key_pressed(KeyCode::Q) {
            if e.documents.iter().any(|doc| doc.dirty) {
                e.pending_quit = true;
            } else {
                e.quit_now = true;
            }
        }
        if let Some((current, result)) = e.calc_rx.as_ref().and_then(|(captured, rx)| {
            rx.try_recv()
                .ok()
                .map(|result| (captured == &e.mol, result))
        }) {
            e.calc_rx = None;
            let auto_mapping = e.auto_map_pending.take();
            let auto_mapping_result = auto_mapping.is_some();
            if current {
                let calculation_failed = result.is_err();
                let cip_assignment = result
                    .as_ref()
                    .is_ok_and(|report| report.starts_with("ASSIGN_STEREO\n"));
                let clipboard_copy = result
                    .as_ref()
                    .is_ok_and(|report| report.starts_with("COPIED_SMILES\n"));
                match result {
                    Ok(report) if auto_mapping.is_some() => {
                        let (source, target, mode) = auto_mapping.as_ref().unwrap();
                        e.apply_auto_atom_mapping(&report, source, target, mode);
                    }
                    Ok(report) if report.starts_with("COPIED_SMILES\n") => {
                        let smiles = report.trim_start_matches("COPIED_SMILES\n");
                        e.info_overlay = Some(format!("SMILES copied to the system clipboard\nStructure notation does not include 2D layout, drawing objects or molecule metadata.\n\n{smiles}"));
                        e.info_scroll = 0;
                        e.status = "Copied isomeric SMILES to clipboard".into();
                    }
                    Ok(report) if report.starts_with("ASSIGN_STEREO\n") => {
                        e.apply_stereo_labels(report.trim_start_matches("ASSIGN_STEREO\n"));
                    }
                    Ok(report) if report.starts_with("IMPORT\t") => {
                        if let Some((path, payload)) = report.split_once('\n') {
                            if let Some((_, block)) = payload.split_once("MOLBLOCK\n") {
                                match parse_mol(block) {
                                    Ok(mol) => {
                                        let filename = path.trim_start_matches("IMPORT\t");
                                        e.open_snapshot_tab(
                                            Snapshot {
                                                mol,
                                                graphics: vec![],
                                            },
                                            Some(filename.into()),
                                        );
                                        e.record_recent(filename);
                                        e.status =
                                            format!("Imported {filename} with optional RDKit");
                                    }
                                    Err(error) => {
                                        e.info_overlay = Some(format!(
                                            "RDKit produced an invalid MOL structure: {error}"
                                        ))
                                    }
                                }
                            }
                        }
                    }
                    Ok(report) if report.contains("MOLBLOCK\n") => {
                        if let Some((metadata, block)) = report.split_once("MOLBLOCK\n") {
                            e.apply_rdkit_molblock(metadata, block)
                        }
                    }
                    Ok(report) => {
                        e.info_overlay = Some(report);
                        e.info_scroll = 0;
                    }
                    Err(message) => {
                        let (notice, status) = background_calculation_error(&message);
                        e.info_overlay = Some(notice);
                        e.status = status;
                    }
                }
                if !calculation_failed && !cip_assignment && !clipboard_copy && !auto_mapping_result
                {
                    e.status = "Calculation finished".into();
                }
            } else {
                e.auto_map_pending = None;
                e.status = "Discarded calculation for an earlier molecule".into();
            }
        }
        let w = screen_width();
        let h = screen_height();
        autosave_clock += get_frame_time();
        if autosave_clock >= 20.0 {
            e.autosave();
            autosave_clock = 0.0;
        }
        let p = vec2(mouse_position().0, mouse_position().1);
        if e.view_3d_open {
            clear_background(surface());
            let w = screen_width();
            let h = screen_height();
            draw_rectangle(0., 0., w, 58., navy());
            let (view_title, rotate_hint, no_coordinates, partial_coordinates) =
                if e.language == Language::Danish {
                    (
                        "MOLECULA  ·  3D-VISNING",
                        "Træk/hjul · R nulstil · A symboler",
                        "Ingen 3D-koordinater. Vælg Beregninger → Generér 3D-koordinater.",
                        "Nogle atomer mangler 3D-koordinater.",
                    )
                } else {
                    (
                        "MOLECULA  ·  3D VIEW",
                        "Drag/wheel · R reset · A symbols",
                        "No 3D coordinates. Choose Calculations → Generate 3D coordinates.",
                        "Some atoms have no 3D coordinates.",
                    )
                };
            draw_text(view_title, 24., 37., 22., WHITE);
            draw_text(
                rotate_hint,
                275.,
                36.,
                15.,
                Color::from_rgba(185, 204, 194, 255),
            );
            if button(
                w - 500.,
                12.,
                140.,
                34.,
                three_d_style_name(e.view_3d_style, e.language),
                false,
            ) {
                e.view_3d_style = (e.view_3d_style + 1) % THREE_D_STYLE_COUNT;
            }
            let symbols_label = if e.language == Language::Danish {
                if e.show_3d_atom_symbols {
                    "Symboler: Til"
                } else {
                    "Symboler: Fra"
                }
            } else if e.show_3d_atom_symbols {
                "Symbols: On"
            } else {
                "Symbols: Off"
            };
            if button(w - 352., 12., 118., 34., symbols_label, false) {
                e.toggle_3d_atom_symbols();
            }
            let reset_label = if e.language == Language::Danish {
                "Nulstil"
            } else {
                "Reset view"
            };
            if button(w - 226., 12., 96., 34., reset_label, false) || is_key_pressed(KeyCode::R) {
                e.rotation_3d = Vec2::ZERO;
                e.zoom_3d = 1.;
                e.last_3d_pointer = None;
                e.status = if e.language == Language::Danish {
                    "3D-visning nulstillet".into()
                } else {
                    "3D view reset".into()
                };
            }
            if is_key_pressed(KeyCode::A) {
                e.toggle_3d_atom_symbols();
            }
            let close_label = if e.language == Language::Danish {
                "Luk  ×"
            } else {
                "Close  ×"
            };
            if button(w - 125., 12., 105., 34., close_label, true)
                || is_key_pressed(KeyCode::Escape)
            {
                e.view_3d_open = false;
                e.last_3d_pointer = None;
                continue;
            }
            let atoms: Vec<Vec3> = e.mol.atoms.iter().filter_map(|a| a.coord3d).collect();
            if !atoms.is_empty() {
                let center =
                    atoms.iter().copied().fold(Vec3::ZERO, |sum, p| sum + p) / atoms.len() as f32;
                let scale = (70. * e.zoom_3d).clamp(12., 160.);
                let origin = vec2(w * 0.5, h * 0.53);
                if is_mouse_button_pressed(MouseButton::Left) && p.y > 60. {
                    e.last_3d_pointer = Some(p);
                }
                if is_mouse_button_down(MouseButton::Left) {
                    if let Some(last) = e.last_3d_pointer {
                        let delta = p - last;
                        e.rotation_3d += vec2(delta.x, delta.y) * 0.008;
                        e.last_3d_pointer = Some(p);
                    }
                } else {
                    e.last_3d_pointer = None;
                }
                let (_, wheel) = mouse_wheel();
                if wheel != 0. {
                    e.zoom_3d = (e.zoom_3d * 1.12_f32.powf(wheel)).clamp(0.2, 5.);
                }
                let projected: Vec<Option<(Vec2, f32)>> = e
                    .mol
                    .atoms
                    .iter()
                    .map(|a| {
                        a.coord3d.map(|q| {
                            let rotated = rotate_3d(q - center, e.rotation_3d);
                            (
                                project_3d(q, center, e.rotation_3d, scale, origin),
                                rotated.z,
                            )
                        })
                    })
                    .collect();
                if e.view_3d_style != 3 {
                    for bond in &e.mol.bonds {
                        if let (Some((a, _)), Some((b, _))) = (
                            projected.get(bond.a).copied().flatten(),
                            projected.get(bond.b).copied().flatten(),
                        ) {
                            draw_line(
                                a.x,
                                a.y,
                                b.x,
                                b.y,
                                if e.view_3d_style == 0 { 1. } else { 3. },
                                three_d_bond_color(e.view_3d_style),
                            );
                        }
                    }
                }
                let mut draw_order: Vec<(usize, Vec2, f32)> = projected
                    .iter()
                    .enumerate()
                    .filter_map(|(index, position)| {
                        position.map(|(point, depth)| (index, point, depth))
                    })
                    .collect();
                draw_order.sort_by(|left, right| left.2.total_cmp(&right.2));
                for (index, at, _) in draw_order {
                    let atom = &e.mol.atoms[index];
                    let radius = three_d_atom_radius(e.view_3d_style, &atom.element, scale);
                    draw_circle(
                        at.x,
                        at.y,
                        radius,
                        three_d_atom_color(e.view_3d_style, &atom.element),
                    );
                    if e.show_3d_atom_symbols {
                        draw_text(&atom.element, at.x + radius + 2., at.y + 5., 14., ink());
                    }
                }
                if e.mol.atoms.iter().any(|atom| atom.coord3d.is_none()) {
                    draw_text(
                        partial_coordinates,
                        24.,
                        82.,
                        16.,
                        Color::from_rgba(196, 112, 55, 255),
                    );
                }
            } else if !e.mol.atoms.is_empty() {
                draw_text(no_coordinates, 24., 88., 17., muted());
            }
            let coordinate_note = if e.language == Language::Danish {
                "3D-koordinater gemmes separat; kamerabevægelser ændrer ikke 2D-strukturen."
            } else {
                "3D coordinates are stored separately; camera movement does not change the 2D structure."
            };
            draw_text(coordinate_note, 24., h - 22., 14., muted());
            next_frame().await;
            continue;
        }
        if e.lesson_open {
            clear_background(surface());
            let w = screen_width();
            let h = screen_height();
            let x = (w - 760.) / 2.;
            let y = (h - 500.) / 2.;
            round_rect(
                x + 4.,
                y + 6.,
                760.,
                500.,
                18.,
                Color::from_rgba(0, 0, 0, 34),
            );
            round_rect(x, y, 760., 500., 18., card_surface());
            draw_rectangle(x, y, 760., 76., navy());
            let title = match e.lesson.as_ref().map(|l| &l.content) {
                Some(LessonContent::Naming { .. }) => "NAMING · GUIDED PRACTICE",
                Some(LessonContent::Reaction { .. }) => "REACTION · BALANCE PRACTICE",
                None => "TEACHING",
            };
            draw_text(title, x + 26., y + 34., 20., WHITE);
            draw_text(
                "Step by step · chemistry is checked against the implemented model",
                x + 26.,
                y + 58.,
                13.,
                Color::from_rgba(187, 211, 197, 255),
            );
            let mut action = None;
            if let Some(lesson) = e.lesson.as_ref() {
                match &lesson.content {
                    LessonContent::Naming { name, steps } => {
                        let total = steps.len();
                        let step = lesson.step.min(total.saturating_sub(1));
                        draw_text(
                            &format!("STEP {} OF {}", step + 1, total),
                            x + 28.,
                            y + 114.,
                            12.,
                            accent(),
                        );
                        draw_text(&steps[step], x + 28., y + 158., 19., ink());
                        draw_text("The highlighted atoms show the supported parent structure. The naming engine covers a limited teaching subset, not general IUPAC nomenclature.",x+28.,y+220.,15.,muted());
                        if step + 1 == total {
                            draw_text(
                                &format!("Assembled name: {name}"),
                                x + 28.,
                                y + 295.,
                                24.,
                                accent(),
                            );
                        }
                        if step > 0 && button(x + 28., y + 430., 105., 42., "← Previous", false) {
                            action = Some("previous");
                        }
                        if step + 1 < total && button(x + 145., y + 430., 105., 42., "Next →", true)
                        {
                            action = Some("next");
                        }
                    }
                    LessonContent::Reaction { feedback } => {
                        let steps=["1 · Separate the reaction sides. Put reactants before the arrow and products after it; set each connected component's role in the inspector.","2 · Set stoichiometric coefficients. Select an atom in a component and cycle its coefficient; all atoms in that connected component update together.","3 · Add a reaction arrow. The editor needs an explicit reaction arrow to identify the reaction step being analyzed.","4 · Check the result. Molecula compares element counts and formal charge. It reports ‘not evaluated’ when the graph is unsupported or incomplete."];
                        let step = lesson.step.min(steps.len() - 1);
                        draw_text(
                            &format!("STEP {} OF {}", step + 1, steps.len()),
                            x + 28.,
                            y + 114.,
                            12.,
                            accent(),
                        );
                        for (i, line) in steps[step].split('\n').enumerate() {
                            draw_text(
                                &line.chars().take(86).collect::<String>(),
                                x + 28.,
                                y + 158. + i as f32 * 28.,
                                17.,
                                ink(),
                            );
                        }
                        if let Some(result) = feedback {
                            draw_rectangle(x + 25., y + 270., 710., 130., panel_surface());
                            for (i, line) in result.lines().take(4).enumerate() {
                                draw_text(
                                    &line.chars().take(90).collect::<String>(),
                                    x + 38.,
                                    y + 298. + i as f32 * 23.,
                                    14.,
                                    if line.starts_with("Correct") {
                                        accent()
                                    } else {
                                        ink()
                                    },
                                );
                            }
                        }
                        if step > 0 && button(x + 28., y + 430., 105., 42., "← Previous", false) {
                            action = Some("previous");
                        }
                        if step + 1 < steps.len()
                            && button(x + 145., y + 430., 105., 42., "Next →", true)
                        {
                            action = Some("next");
                        }
                        if step + 1 == steps.len()
                            && button(x + 145., y + 430., 145., 42., "Check reaction", true)
                        {
                            action = Some("check");
                        }
                    }
                }
            }
            if button(x + 420., y + 430., 155., 42., "Edit in canvas", false) {
                action = Some("edit");
            }
            if button(x + 590., y + 430., 140., 42., "End lesson", false) {
                action = Some("end");
            }
            if is_key_pressed(KeyCode::Escape) {
                action = Some("edit");
            }
            match action {
                Some("previous") => {
                    if let Some(l) = &mut e.lesson {
                        l.step = l.step.saturating_sub(1)
                    }
                }
                Some("next") => {
                    if let Some(l) = &mut e.lesson {
                        l.step += 1
                    }
                }
                Some("check") => e.check_reaction_lesson(),
                Some("edit") => e.lesson_open = false,
                Some("end") => {
                    e.lesson = None;
                    e.lesson_open = false;
                    e.naming_highlight.clear()
                }
                _ => {}
            }
            next_frame().await;
            continue;
        }
        let world_p = e.to_world(p, w, h);
        let display = e.display_snapshot(w, h);
        let analysis = chemistry::analyze(&e.mol);
        if (is_mouse_button_pressed(MouseButton::Middle)
            || (is_key_down(KeyCode::Space) && is_mouse_button_pressed(MouseButton::Left)))
            && canvas(p)
        {
            e.panning = true;
            e.drag_start = Some(p);
        }
        if e.panning
            && (is_mouse_button_down(MouseButton::Middle)
                || (is_key_down(KeyCode::Space) && is_mouse_button_down(MouseButton::Left)))
        {
            if let Some(last) = e.drag_start {
                e.pan += p - last;
                e.drag_start = Some(p);
            }
        }
        if e.panning
            && !is_mouse_button_down(MouseButton::Middle)
            && !is_mouse_button_down(MouseButton::Left)
        {
            e.panning = false;
            e.drag_start = None;
        }
        let (wheel_x, wheel_y) = mouse_wheel();
        if e.info_overlay.is_some() && wheel_y != 0. {
            let line_count = wrap_info_lines(e.info_overlay.as_deref().unwrap_or(""), 75).len();
            let visible = ((h - 166.).max(24.) / 24.).floor() as usize;
            let max_scroll = line_count.saturating_sub(visible);
            e.info_scroll = if wheel_y < 0. {
                e.info_scroll
                    .saturating_add(wheel_y.abs().ceil() as usize)
                    .min(max_scroll)
            } else {
                e.info_scroll.saturating_sub(wheel_y.ceil() as usize)
            };
        } else if (wheel_x != 0. || wheel_y != 0.) && canvas(p) {
            e.adjust_zoom(1.12_f32.powf(wheel_y), p, w, h);
            e.pan.x -= wheel_x * 26.;
        }
        clear_background(surface());
        draw_rectangle(0., 0., w, 39., navy());
        draw_rectangle(0., 39., w, TOP - 39., card_surface());
        draw_line(0., TOP, w, TOP, 1., Color::from_rgba(218, 226, 238, 255));
        draw_line(
            16.,
            22.,
            27.,
            13.,
            2.5,
            Color::from_rgba(128, 219, 169, 255),
        );
        draw_line(
            27.,
            13.,
            38.,
            22.,
            2.5,
            Color::from_rgba(128, 219, 169, 255),
        );
        draw_line(
            38.,
            22.,
            27.,
            29.,
            2.5,
            Color::from_rgba(128, 219, 169, 255),
        );
        for (cx, cy) in [(16., 22.), (27., 13.), (38., 22.), (27., 29.)] {
            draw_circle(cx, cy, 3.5, WHITE);
        }
        draw_text("molecula", 52., 27., 21., WHITE);
        draw_text(
            "STRUCTURE STUDIO",
            165.,
            25.,
            10.,
            Color::from_rgba(150, 168, 194, 255),
        );
        for (i, name) in MENU_NAMES.iter().enumerate() {
            if button(
                304. + i as f32 * 94.,
                3.,
                90.,
                32.,
                name,
                e.menu_open == Some(i),
            ) {
                e.menu_open = if e.menu_open == Some(i) {
                    None
                } else {
                    Some(i)
                };
            }
        }
        let menu = [
            ("New", 0),
            ("Open", 1),
            ("Save", 2),
            ("Undo", 3),
            ("Redo", 4),
            ("Fit", 5),
            ("Elements", 6),
            ("Templates", 7),
        ];
        round_rect(
            14.,
            45.,
            610.,
            34.,
            10.,
            Color::from_rgba(246, 248, 252, 255),
        );
        let mut mx = 17.;
        for (label, id) in menu {
            if button(mx, 44., 70., 31., label, false) {
                match id {
                    0 => e.new_document(),
                    1 => e.file_dialog = Some(false),
                    2 => e.file_dialog = Some(true),
                    3 => e.undo(),
                    4 => e.redo(),
                    5 => {
                        let n = e.mol.atoms.len();
                        if n > 0 {
                            let center =
                                e.mol.atoms.iter().fold(vec2(0., 0.), |acc, a| acc + a.pos)
                                    / n as f32;
                            let dest = vec2((w + LEFT) / 2., (h + TOP - BOTTOM) / 2.);
                            for a in &mut e.mol.atoms {
                                a.pos += dest - center;
                            }
                        }
                    }
                    6 => e.periodic_open = true,
                    7 => e.templates_open = true,
                    _ => {}
                }
            }
            mx += 76.;
        }
        // Document tabs live in the dedicated strip below the toolbar.
        draw_rectangle(0., 80., w, TOP - 80., surface());
        let mut tx = LEFT + 8.;
        let tab_hits: Vec<(usize, f32, f32)> = e
            .documents
            .iter()
            .enumerate()
            .map(|(i, d)| {
                let label = format!("{}{}", if d.dirty { "● " } else { "" }, d.title);
                let width = (measure_text(&label, None, 15, 1.).width + 48.).clamp(130., 220.);
                round_rect(
                    tx,
                    84.,
                    width,
                    29.,
                    8.,
                    if i == e.active_document {
                        Color::from_rgba(222, 244, 234, 255)
                    } else {
                        card_surface()
                    },
                );
                draw_text(&label, tx + 12., 104., 15., ink());
                let hit = (i, tx, width);
                tx += width + 6.;
                hit
            })
            .collect();
        for (i, x, width) in tab_hits {
            if is_mouse_button_pressed(MouseButton::Left)
                && p.x >= x
                && p.x <= x + width
                && p.y >= 84.
                && p.y <= 113.
            {
                e.switch_document(i);
            }
            if button(x + width - 29., 84., 25., 29., "×", false) {
                e.request_close_document(i);
            }
        }
        if is_mouse_button_pressed(MouseButton::Left)
            && p.x >= tx
            && p.x < tx + 32.
            && p.y >= 84.
            && p.y <= 113.
        {
            e.new_document();
        }
        round_rect(tx, 84., 30., 29., 8., card_surface());
        draw_text("+", tx + 10., 105., 20., ink());
        draw_rectangle(0., TOP, LEFT, h - TOP - BOTTOM, navy());
        draw_text(
            "TOOLS",
            17.,
            TOP + 20.,
            10.,
            Color::from_rgba(132, 152, 183, 255),
        );
        let tools = [
            ("Sel", Tool::Select),
            ("Del", Tool::Erase),
            ("Chain", Tool::Chain),
            ("1", Tool::Bond(1)),
            ("=", Tool::Bond(2)),
            ("3", Tool::Bond(3)),
            ("Ar", Tool::StyledBond(BondStyle::Aromatic)),
            ("C", Tool::Atom),
            ("R6", Tool::Ring(6)),
            ("R5", Tool::Ring(5)),
            ("->", Tool::Arrow),
            ("e2", Tool::Electron(2)),
            ("e1", Tool::Electron(1)),
            ("+R", Tool::Plus),
            ("Txt", Tool::Text),
            ("+", Tool::Charge(1)),
            ("-", Tool::Charge(-1)),
            ("Wd", Tool::StyledBond(BondStyle::Wedge)),
            ("Hs", Tool::StyledBond(BondStyle::Hash)),
            ("Wv", Tool::StyledBond(BondStyle::Wavy)),
        ];
        for (i, (label, t)) in tools.iter().enumerate() {
            let y = TOP + 30. + i as f32 * 35.;
            if tool_button(18., y, 52., 33., *t, e.tool == *t) {
                e.pending_arrow = None;
                e.tool = *t;
                e.status = format!("Tool: {label}");
            }
        }
        if button(646., 44., 40., 31., "−", false) {
            e.adjust_zoom(1. / 1.2, p, w, h);
        }
        draw_text(
            &format!("{}%", (e.zoom * 100.).round() as i32),
            692.,
            65.,
            14.,
            ink(),
        );
        if button(738., 44., 40., 31., "+", false) {
            e.adjust_zoom(1.2, p, w, h);
        }
        if button(784., 44., 72., 31., "Fit view", false) {
            e.menu_action("Fit");
        }
        draw_text("Wheel: zoom · Space+drag: pan", 865., 64., 13., muted());
        let grid = Color::from_rgba(224, 231, 242, 255);
        let mut x = LEFT + 26.;
        while x < w - RIGHT {
            let mut y = TOP + 25.;
            while y < h - BOTTOM {
                draw_circle(x, y, 1., grid);
                y += 28.;
            }
            x += 28.;
        }
        if e.mol.atoms.is_empty() && e.graphics.is_empty() {
            let cx = (LEFT + w - RIGHT) / 2.;
            let cy = (TOP + h - BOTTOM) / 2.;
            round_rect(
                cx - 47.,
                cy - 115.,
                94.,
                94.,
                24.,
                Color::from_rgba(231, 238, 252, 255),
            );
            let ghost = Color::from_rgba(113, 148, 220, 255);
            for k in 0..6 {
                let a = std::f32::consts::TAU * k as f32 / 6.;
                let b = std::f32::consts::TAU * (k + 1) as f32 / 6.;
                draw_line(
                    cx + 26. * a.cos(),
                    cy - 68. + 26. * a.sin(),
                    cx + 26. * b.cos(),
                    cy - 68. + 26. * b.sin(),
                    2.,
                    ghost,
                );
            }
            draw_text("Start drawing chemistry", cx - 148., cy + 22., 25., ink());
            draw_text(
                "Choose a tool or drag on the canvas",
                cx - 154.,
                cy + 53.,
                16.,
                muted(),
            );
            round_rect(
                cx - 105.,
                cy + 74.,
                210.,
                30.,
                15.,
                Color::from_rgba(235, 241, 252, 255),
            );
            draw_text(
                "1  single bond    •    C  carbon",
                cx - 91.,
                cy + 95.,
                13.,
                accent(),
            );
        }
        for (i, b) in display.mol.bonds.iter().enumerate() {
            if e.reaction_step_view.is_some_and(|step| {
                display.mol.atoms[b.a].reaction_step != step
                    || display.mol.atoms[b.b].reaction_step != step
            }) {
                continue;
            }
            let a = display.mol.atoms[b.a].pos;
            let c = display.mol.atoms[b.b].pos;
            draw_chem_bond(
                a,
                c,
                b,
                if e.selected_bonds.contains(&i) || e.selected_bond == Some(i) {
                    accent()
                } else if e.tool == Tool::Erase && e.hit_bond(world_p) == Some(i) {
                    RED
                } else if e.reaction_highlight_bonds.contains(&i) {
                    Color::from_rgba(220, 76, 68, 255)
                } else {
                    DARKGRAY
                },
            );
        }
        let mut reaction_arrow_step = 0u16;
        for (gi, g) in display.graphics.iter().enumerate() {
            if let Graphic::Arrow { kind, .. } = g {
                if matches!(kind.as_str(), "reaction" | "equilibrium" | "retrosynthesis") {
                    reaction_arrow_step += 1;
                    if e.reaction_step_view
                        .is_some_and(|step| step != reaction_arrow_step)
                    {
                        continue;
                    }
                }
            }
            match g {
                Graphic::Arrow {
                    start,
                    end,
                    kind,
                    above,
                    below,
                } => {
                    draw_reaction_arrow(*start, *end, kind, DARKGRAY);
                    let mid = (*start + *end) / 2.;
                    draw_text(above, mid.x, mid.y - 12., 18., DARKGRAY);
                    draw_text(below, mid.x, mid.y + 25., 18., DARKGRAY);
                }
                Graphic::ElectronArrow {
                    start,
                    control,
                    end,
                    electrons,
                } => draw_electron_arrow(*start, *control, *end, *electrons, DARKGRAY),
                Graphic::Text {
                    pos, text, size, ..
                } => {
                    draw_text(text, pos.x, pos.y, *size, DARKGRAY);
                }
                Graphic::Plus { pos } => {
                    draw_text("+", pos.x, pos.y, 30., DARKGRAY);
                }
                Graphic::Shape { kind, start, end } => draw_shape(*kind, *start, *end, DARKGRAY),
            }
            if e.selected_graphics.contains(&gi) || e.selected_graphic == Some(gi) {
                let center = match g {
                    Graphic::Arrow { start, end, .. } => (*start + *end) / 2.,
                    Graphic::ElectronArrow { control, .. } => *control,
                    Graphic::Text { pos, .. } | Graphic::Plus { pos } => *pos,
                    Graphic::Shape { start, end, .. } => (*start + *end) / 2.,
                };
                draw_circle_lines(center.x, center.y, 16., 1.5, accent());
            }
        }
        for (i, a) in display.mol.atoms.iter().enumerate() {
            if e.reaction_step_view
                .is_some_and(|step| a.reaction_step != step)
            {
                continue;
            }
            if e.selected_atoms.contains(&i) || e.naming_highlight.contains(&i) {
                draw_circle(a.pos.x, a.pos.y, 17., Color::from_rgba(196, 223, 253, 255));
            }
            let implicit_h = analysis.implicit_h.get(i).copied().flatten();
            if show_label(&e.mol, i)
                || (e.show_implicit_hydrogens && implicit_h.is_some_and(|h| h > 0))
            {
                let label = atom_display_label(a, implicit_h, e.show_implicit_hydrogens);
                let d = measure_text(&label, None, 21, 1.);
                draw_rectangle(
                    a.pos.x - d.width / 2. - 3.,
                    a.pos.y - 16.,
                    d.width + 6.,
                    24.,
                    Color::from_rgba(248, 249, 251, 255),
                );
                draw_text(
                    &label,
                    a.pos.x - d.width / 2.,
                    a.pos.y + 5.,
                    21.,
                    if a.element == "O" {
                        RED
                    } else if a.element == "N" {
                        accent()
                    } else {
                        DARKGRAY
                    },
                );
            } else if e.tool == Tool::Select && e.hit_atom(world_p) == Some(i) {
                draw_circle_lines(a.pos.x, a.pos.y, 7., 1., accent());
            }
            if let Some(n) = a.atom_map {
                draw_text(&format!("{n}"), a.pos.x + 10., a.pos.y + 17., 13., GRAY);
            }
            if let Some(label @ ("R" | "S")) = a.stereo.as_deref() {
                draw_text(
                    label,
                    a.pos.x + 11.,
                    a.pos.y - 13.,
                    13.,
                    Color::from_rgba(149, 92, 24, 255),
                );
            }
        }
        if let Some(start) = e.selection_box {
            let start = e.to_screen(start, w, h);
            let left = start.x.min(p.x);
            let top = start.y.min(p.y);
            let width = (start.x - p.x).abs();
            let height = (start.y - p.y).abs();
            draw_rectangle(left, top, width, height, Color::from_rgba(36, 139, 94, 28));
            draw_rectangle_lines(left, top, width, height, 1.5, accent());
        }
        if let Some(start) = e.drag_start.or(e.pending_arrow) {
            if matches!(e.tool, Tool::Bond(_) | Tool::StyledBond(_)) && canvas(p) {
                let start_world = start;
                let origin = e
                    .hit_atom(start_world)
                    .map(|i| e.mol.atoms[i].pos)
                    .unwrap_or(start_world);
                let end_world = if world_p.distance(origin) < 14. {
                    e.hit_atom(start)
                        .map(|i| e.suggest_endpoint(i))
                        .unwrap_or(origin + vec2(BOND, 0.))
                } else {
                    world_p
                };
                let preview = Bond {
                    a: 0,
                    b: 1,
                    order: match e.tool {
                        Tool::Bond(o) => o,
                        _ => 1,
                    },
                    style: match e.tool {
                        Tool::StyledBond(s) => s,
                        _ => BondStyle::Normal,
                    },
                    stereo: None,
                    query: None,
                };
                let origin = e.to_screen(origin, w, h);
                let end = e.to_screen(end_world, w, h);
                draw_chem_bond(origin, end, &preview, accent());
                draw_circle_lines(end.x, end.y, 6., 1., accent());
            } else if e.tool == Tool::Chain && canvas(p) {
                let origin = e
                    .hit_atom(start)
                    .map(|i| e.mol.atoms[i].pos)
                    .unwrap_or(start);
                let delta = world_p - origin;
                let length = delta.length();
                let dir = if length > 1. {
                    delta / length
                } else {
                    vec2(1., 0.)
                };
                let count = if length < 18. {
                    3
                } else {
                    ((length / BOND).round() as usize).clamp(2, 10)
                };
                let end = if length < 18. {
                    origin + dir * BOND * count as f32
                } else {
                    world_p
                };
                let normal = vec2(-dir.y, dir.x);
                let mut prev = origin;
                for k in 1..=count {
                    let t = k as f32 / count as f32;
                    let zig = if k == count {
                        0.
                    } else if k % 2 == 1 {
                        -1.
                    } else {
                        1.
                    };
                    let next = origin + (end - origin) * t + normal * zig * (BOND * 0.22);
                    let prev_s = e.to_screen(prev, w, h);
                    let next_s = e.to_screen(next, w, h);
                    draw_line(prev_s.x, prev_s.y, next_s.x, next_s.y, 2., accent());
                    prev = next;
                }
            } else if e.tool == Tool::Arrow && canvas(p) {
                draw_reaction_arrow(
                    e.to_screen(start, w, h),
                    e.to_screen(world_p, w, h),
                    &e.arrow_kind,
                    accent(),
                );
            } else if let Tool::Shape(kind) = e.tool {
                if canvas(p) {
                    draw_shape(
                        kind,
                        e.to_screen(start, w, h),
                        e.to_screen(world_p, w, h),
                        accent(),
                    );
                }
            } else if let Tool::Electron(n) = e.tool {
                if canvas(p) {
                    let s = e.to_screen(start, w, h);
                    let ep = e.to_screen(world_p, w, h);
                    draw_electron_arrow(s, electron_control(s, ep), ep, n, accent());
                }
            }
        }
        if matches!(e.tool, Tool::Template(_)) && canvas(p) {
            if let Some(preview) = &e.template_preview {
                let placement_center = e.template_placement_center(world_p);
                let ghost = Color::new(accent().r, accent().g, accent().b, 0.58);
                let screen_pos = |atom: &Atom| {
                    let rotated =
                        rotate_template_point(atom.pos, vec2(0., 0.), e.template_rotation);
                    e.to_screen(placement_center + rotated, w, h)
                };
                for bond in &preview.bonds {
                    if let (Some(a), Some(b)) =
                        (preview.atoms.get(bond.a), preview.atoms.get(bond.b))
                    {
                        draw_chem_bond(screen_pos(a), screen_pos(b), bond, ghost);
                    }
                }
                for atom in &preview.atoms {
                    let at = screen_pos(atom);
                    let label = atom_display_label(atom, None, false);
                    if !label.is_empty() {
                        let dims = measure_text(&label, None, 16, 1.0);
                        draw_rectangle(
                            at.x - dims.width * 0.5 - 2.,
                            at.y - 10.,
                            dims.width + 4.,
                            14.,
                            Color::new(1., 1., 1., 0.78),
                        );
                        draw_text(&label, at.x - dims.width * 0.5, at.y + 2., 16., ghost);
                    } else if preview.bonds.is_empty() {
                        draw_circle_lines(at.x, at.y, 4., 1.5, ghost);
                    }
                }
                if let (Some(target), Some((anchor_index, _))) =
                    (e.template_attachment_target, e.template_anchor())
                {
                    if let (Some(target_atom), Some(anchor_atom)) =
                        (e.mol.atoms.get(target), preview.atoms.get(anchor_index))
                    {
                        let linker = Bond {
                            a: 0,
                            b: 1,
                            order: 1,
                            style: BondStyle::Normal,
                            stereo: None,
                            query: None,
                        };
                        draw_chem_bond(
                            e.to_screen(target_atom.pos, w, h),
                            screen_pos(anchor_atom),
                            &linker,
                            ghost,
                        );
                    }
                }
            }
        }
        let panel_x = w - RIGHT;
        draw_rectangle(panel_x, TOP, RIGHT, h - TOP - BOTTOM, panel_surface());
        draw_line(
            panel_x,
            TOP,
            panel_x,
            h - BOTTOM,
            1.,
            Color::from_rgba(218, 226, 238, 255),
        );
        draw_text("INSPECTOR", panel_x + 16., TOP + 23., 11., muted());
        round_rect(panel_x + 11., TOP + 38., 228., 545., 12., card_surface());
        draw_text("Properties", panel_x + 130., TOP + 24., 16., ink());
        if e.selected_atoms.len() > 1 {
            let mut elements = std::collections::BTreeMap::<String, usize>::new();
            let mut charge = 0i32;
            for &i in &e.selected_atoms {
                if let Some(atom) = e.mol.atoms.get(i) {
                    *elements.entry(atom.element.clone()).or_default() += 1;
                    charge += atom.charge as i32;
                }
            }
            let element_summary = elements
                .iter()
                .map(|(element, count)| format!("{element} × {count}"))
                .collect::<Vec<_>>()
                .join(", ");
            let internal_bonds = e
                .mol
                .bonds
                .iter()
                .filter(|b| e.selected_atoms.contains(&b.a) && e.selected_atoms.contains(&b.b))
                .count();
            draw_text("Canvas selection", panel_x + 16., TOP + 61., 18., DARKGRAY);
            draw_text(
                &format!("{} atoms", e.selected_atoms.len()),
                panel_x + 16.,
                TOP + 94.,
                16.,
                DARKGRAY,
            );
            draw_text(
                &element_summary.chars().take(28).collect::<String>(),
                panel_x + 16.,
                TOP + 124.,
                14.,
                DARKGRAY,
            );
            draw_text(
                &format!("Internal bonds: {internal_bonds}"),
                panel_x + 16.,
                TOP + 154.,
                15.,
                DARKGRAY,
            );
            draw_text(
                &format!("Formal charge: {charge:+}"),
                panel_x + 16.,
                TOP + 184.,
                15.,
                DARKGRAY,
            );
            draw_text(
                "Drag a box on the canvas to select",
                panel_x + 16.,
                TOP + 222.,
                12.,
                muted(),
            );
        } else if let Some(i) = e.selected.filter(|i| *i < e.mol.atoms.len()) {
            let atom = e.mol.atoms[i].clone();
            let role = atom.reaction_role.as_deref().unwrap_or("auto").to_string();
            let step = atom.reaction_step;
            let coefficient = atom.coefficient;
            draw_text(
                &format!(
                    "Atom {}: {}{}",
                    i + 1,
                    atom.alias.as_deref().unwrap_or(&atom.element),
                    atom.alias
                        .as_ref()
                        .map(|_| format!(" ({})", atom.element))
                        .unwrap_or_default()
                ),
                panel_x + 16.,
                TOP + 61.,
                18.,
                DARKGRAY,
            );
            draw_text(
                &format!("Charge: {}", atom.charge),
                panel_x + 16.,
                TOP + 91.,
                17.,
                DARKGRAY,
            );
            draw_text(
                &format!(
                    "Isotope: {}",
                    atom.isotope
                        .map(|n| n.to_string())
                        .unwrap_or_else(|| "natural".into())
                ),
                panel_x + 16.,
                TOP + 121.,
                17.,
                DARKGRAY,
            );
            draw_text(
                &format!(
                    "Map: {}",
                    atom.atom_map
                        .map(|n| n.to_string())
                        .unwrap_or_else(|| "none".into())
                ),
                panel_x + 16.,
                TOP + 151.,
                17.,
                DARKGRAY,
            );
            draw_text(
                &format!(
                    "Reaction: {} · step {} · ×{}",
                    atom.reaction_role.as_deref().unwrap_or("auto"),
                    atom.reaction_step,
                    atom.coefficient
                ),
                panel_x + 16.,
                TOP + 181.,
                14.,
                DARKGRAY,
            );
            if button(panel_x + 16., TOP + 203., 100., 33., "Charge +", false) {
                e.checkpoint();
                e.mol.atoms[i].charge += 1;
            }
            if button(panel_x + 126., TOP + 203., 100., 33., "Charge -", false) {
                e.checkpoint();
                e.mol.atoms[i].charge -= 1;
            }
            if button(panel_x + 16., TOP + 246., 100., 33., "13C", false) {
                e.checkpoint();
                e.mol.atoms[i].element = "C".into();
                e.mol.atoms[i].isotope = Some(13);
            }
            if button(panel_x + 126., TOP + 246., 100., 33., "Clear iso", false) {
                e.checkpoint();
                e.mol.atoms[i].isotope = None;
            }
            if button(
                panel_x + 16.,
                TOP + 280.,
                104.,
                29.,
                &format!("Role: {role}"),
                false,
            ) {
                let next = match role.as_str() {
                    "auto" => "reactant",
                    "reactant" => "product",
                    "product" => "agent",
                    _ => "",
                };
                e.set_reaction_component(Some(next), None, None);
            }
            if button(
                panel_x + 124.,
                TOP + 280.,
                55.,
                29.,
                &format!("Step {step}"),
                false,
            ) {
                e.set_reaction_component(None, Some(step % 3 + 1), None);
            }
            if button(
                panel_x + 183.,
                TOP + 280.,
                44.,
                29.,
                &format!("×{coefficient}"),
                false,
            ) {
                e.set_reaction_component(None, None, Some(coefficient % 5 + 1));
            }
            let valence_label = e.mol.atoms[i]
                .valence_override
                .map(|v| v.to_string())
                .unwrap_or_else(|| "auto".into());
            if button(
                panel_x + 16.,
                TOP + 314.,
                210.,
                29.,
                &format!("Valence: {valence_label}  ↻"),
                false,
            ) {
                e.checkpoint();
                e.mol.atoms[i].valence_override =
                    cycle_valence_override(e.mol.atoms[i].valence_override);
            }
            let radical_label = match atom.radical {
                1 => "singlet",
                2 => "doublet",
                3 => "triplet",
                _ => "none",
            };
            if button(
                panel_x + 16.,
                TOP + 348.,
                210.,
                29.,
                &format!(
                    "{}: {}  ↻",
                    localized(e.language, "Radical"),
                    localized(e.language, radical_label)
                ),
                false,
            ) {
                e.cycle_atom_radical(i);
            }
            let implicit_h_label = atom
                .implicit_h_override
                .map(|count| count.to_string())
                .unwrap_or_else(|| "auto".into());
            if button(
                panel_x + 16.,
                TOP + 382.,
                210.,
                29.,
                &format!(
                    "{}: {}  ↻",
                    localized(e.language, "Implicit H"),
                    localized(e.language, &implicit_h_label)
                ),
                false,
            ) {
                e.cycle_implicit_h_override(i);
            }
            draw_text(
                &localized(e.language, "Query predicates"),
                panel_x + 16.,
                TOP + 431.,
                15.,
                DARKGRAY,
            );
            if let Some(query) = &atom.query {
                draw_text(
                    &format!("{} {query}", localized(e.language, "Query:")),
                    panel_x + 16.,
                    TOP + 452.,
                    13.,
                    muted(),
                );
            }
            let query_y = TOP + 468.;
            for (n, (x, label, predicate)) in [
                (panel_x + 16., "Aromatic", ".a"),
                (panel_x + 126., "Aliphatic", ".A"),
                (panel_x + 16., "Unsaturated", ".u"),
                (panel_x + 126., "Any atom", "*"),
                (panel_x + 16., "Heteroatom", "Q"),
            ]
            .into_iter()
            .enumerate()
            {
                let y = TOP + 470. + (n / 2) as f32 * 34.;
                if button(x, y, 100., 29., &localized(e.language, label), false) {
                    e.apply_quick(shortcuts::QuickAction::Query(predicate.into()));
                }
            }
            if button(
                panel_x + 126.,
                query_y + 68.,
                100.,
                29.,
                &localized(e.language, "Clear query"),
                false,
            ) {
                e.clear_atom_query(i);
            }
        } else if let Some(i) = e.selected_bond.filter(|i| *i < e.mol.bonds.len()) {
            let bond = e.mol.bonds[i].clone();
            let a = &e.mol.atoms[bond.a];
            let b = &e.mol.atoms[bond.b];
            draw_text("Selected bond", panel_x + 16., TOP + 61., 18., DARKGRAY);
            draw_text(
                &format!("{} — {}", a.element, b.element),
                panel_x + 16.,
                TOP + 96.,
                16.,
                DARKGRAY,
            );
            draw_text(
                &format!("Bond order: {}", bond.order),
                panel_x + 16.,
                TOP + 128.,
                16.,
                DARKGRAY,
            );
            if button(
                panel_x + 16.,
                TOP + 138.,
                210.,
                31.,
                "Make double (=)",
                bond.order == 2,
            ) {
                e.set_bond_order(i, 2);
            }
            if button(panel_x + 16., TOP + 176., 100., 31., "Order −", false) {
                e.set_bond_order(i, bond.order.saturating_sub(1).max(1));
            }
            if button(panel_x + 126., TOP + 176., 100., 31., "Order +", false) {
                e.set_bond_order(i, if bond.order == 4 { 1 } else { bond.order + 1 });
            }
            if button(
                panel_x + 16.,
                TOP + 214.,
                210.,
                31.,
                &format!("Style: {:?}  ↻", bond.style),
                false,
            ) {
                e.cycle_bond_style(i);
            }
            if bond.order == 2 {
                if button(
                    panel_x + 16.,
                    TOP + 252.,
                    210.,
                    31.,
                    &format!(
                        "E/Z: {}  ↻",
                        bond.stereo.as_deref().unwrap_or("unspecified")
                    ),
                    false,
                ) {
                    e.cycle_double_bond_stereo(i);
                }
            } else {
                draw_text(
                    &localized(e.language, "E/Z applies to double bonds"),
                    panel_x + 16.,
                    TOP + 273.,
                    13.,
                    muted(),
                );
            }
            let topology_label = match bond
                .query
                .as_deref()
                .unwrap_or("")
                .split(';')
                .find(|term| matches!(*term, "topology-ring" | "topology-chain"))
            {
                Some("topology-ring") => "Ring only",
                Some("topology-chain") => "Chain only",
                _ => "Any topology",
            };
            if button(
                panel_x + 16.,
                TOP + 290.,
                210.,
                31.,
                &format!("Query topology: {topology_label}  ↻"),
                false,
            ) {
                e.cycle_bond_topology(i);
            }
            if button(panel_x + 16., TOP + 328., 210., 33., "Delete bond", false) {
                e.checkpoint();
                e.mol.bonds.remove(i);
                e.selected_bond = None;
            }
        } else if let Some(gi) = e.selected_graphic.filter(|i| *i < e.graphics.len()) {
            let is_arrow = matches!(e.graphics[gi], Graphic::Arrow { .. });
            let is_electron = matches!(e.graphics[gi], Graphic::ElectronArrow { .. });
            draw_text(
                if is_arrow {
                    "Reaction arrow"
                } else if is_electron {
                    "Electron-flow arrow"
                } else {
                    "Graphic object"
                },
                panel_x + 16.,
                TOP + 61.,
                18.,
                DARKGRAY,
            );
            if is_arrow {
                if button(
                    panel_x + 16.,
                    TOP + 95.,
                    210.,
                    35.,
                    "Edit text above",
                    false,
                ) {
                    if let Graphic::Arrow { above, .. } = &e.graphics[gi] {
                        e.text_buffer = above.clone();
                    }
                    e.text_target = Some(TextTarget::ArrowAbove(gi));
                }
                if button(
                    panel_x + 16.,
                    TOP + 140.,
                    210.,
                    35.,
                    "Edit text below",
                    false,
                ) {
                    if let Graphic::Arrow { below, .. } = &e.graphics[gi] {
                        e.text_buffer = below.clone();
                    }
                    e.text_target = Some(TextTarget::ArrowBelow(gi));
                }
                if button(
                    panel_x + 16.,
                    TOP + 185.,
                    210.,
                    35.,
                    "Change arrow type",
                    false,
                ) {
                    e.checkpoint();
                    if let Graphic::Arrow { kind, .. } = &mut e.graphics[gi] {
                        *kind = match kind.as_str() {
                            "reaction" => "equilibrium",
                            "equilibrium" => "resonance",
                            "resonance" => "retrosynthesis",
                            _ => "reaction",
                        }
                        .into();
                    }
                }
                if let Graphic::Arrow { kind, .. } = &e.graphics[gi] {
                    draw_text(
                        &format!("Type: {kind}"),
                        panel_x + 20.,
                        TOP + 245.,
                        14.,
                        muted(),
                    );
                }
            } else if let Graphic::ElectronArrow { electrons, .. } = &e.graphics[gi] {
                draw_text(
                    if *electrons == 1 {
                        "One electron (fishhook)"
                    } else {
                        "Electron pair"
                    },
                    panel_x + 16.,
                    TOP + 96.,
                    15.,
                    muted(),
                );
            }
        } else {
            let cx = panel_x + 125.;
            draw_circle(cx, TOP + 131., 29., Color::from_rgba(237, 242, 252, 255));
            draw_circle_lines(
                cx,
                TOP + 131.,
                12.,
                2.,
                Color::from_rgba(129, 155, 217, 255),
            );
            draw_line(
                cx + 9.,
                TOP + 141.,
                cx + 19.,
                TOP + 151.,
                2.,
                Color::from_rgba(129, 155, 217, 255),
            );
            draw_text("Nothing selected", panel_x + 49., TOP + 187., 17., ink());
            draw_text(
                "Drag a box to select atoms",
                panel_x + 40.,
                TOP + 215.,
                13.,
                muted(),
            );
            draw_text(
                "or click a bond/object",
                panel_x + 44.,
                TOP + 234.,
                13.,
                muted(),
            );
        }
        round_rect(panel_x + 11., TOP + 350., 228., 245., 12., card_surface());
        draw_text("Whole document", panel_x + 25., TOP + 378., 20., ink());
        draw_text(
            &format!("Formula: {}", analysis.formula),
            panel_x + 16.,
            TOP + 407.,
            17.,
            DARKGRAY,
        );
        draw_text(
            &format!(
                "Molar: {} g/mol",
                analysis
                    .molar_mass
                    .map(|x| format!("{x:.3}"))
                    .unwrap_or_else(|| "unavailable".into())
            ),
            panel_x + 16.,
            TOP + 433.,
            16.,
            DARKGRAY,
        );
        draw_text(
            &format!(
                "Exact: {} Da",
                analysis
                    .exact_mass
                    .map(|x| format!("{x:.4}"))
                    .unwrap_or_else(|| "unavailable".into())
            ),
            panel_x + 16.,
            TOP + 459.,
            16.,
            DARKGRAY,
        );
        draw_text(
            &format!(
                "Selected implicit H: {}",
                e.selected
                    .and_then(|i| analysis.implicit_h.get(i).copied().flatten())
                    .map(|n| n.to_string())
                    .unwrap_or_else(|| "-".into())
            ),
            panel_x + 16.,
            TOP + 486.,
            15.,
            DARKGRAY,
        );
        draw_text(
            &format!("Warnings: {}", analysis.warnings.len()),
            panel_x + 16.,
            TOP + 515.,
            16.,
            if analysis.warnings.is_empty() {
                DARKGRAY
            } else {
                RED
            },
        );
        for (j, msg) in analysis.warnings.iter().take(3).enumerate() {
            let short: String = msg.chars().take(28).collect();
            draw_text(&short, panel_x + 16., TOP + 543. + j as f32 * 21., 14., RED);
        }
        if h - BOTTOM > TOP + 678. {
            round_rect(panel_x + 11., TOP + 605., 228., 95., 12., card_surface());
            draw_text("GRAPH", panel_x + 25., TOP + 627., 11., muted());
            draw_text(
                &format!(
                    "{} atoms  ·  {} bonds",
                    analysis.atom_count, analysis.bond_count
                ),
                panel_x + 25.,
                TOP + 653.,
                16.,
                ink(),
            );
            draw_text(
                &format!("Formal charge  {:+}", analysis.formal_charge),
                panel_x + 25.,
                TOP + 681.,
                15.,
                muted(),
            );
        }
        draw_rectangle(0., h - BOTTOM, w, BOTTOM, card_surface());
        draw_line(
            0.,
            h - BOTTOM,
            w,
            h - BOTTOM,
            1.,
            Color::from_rgba(218, 226, 238, 255),
        );
        draw_text("ELEMENTS", 17., h - BOTTOM + 32., 11., muted());
        let elements = [
            "C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H", "B", "Si",
        ];
        let mut ex = 88.;
        for el in elements {
            if button(ex, h - BOTTOM + 10., 48., 38., el, e.element == el) {
                e.element = el.into();
                e.tool = Tool::Atom;
            }
            ex += 53.;
        }
        draw_line(
            0.,
            h - 31.,
            w,
            h - 31.,
            1.,
            Color::from_rgba(237, 241, 248, 255),
        );
        draw_circle(22., h - 15., 4., Color::from_rgba(54, 187, 143, 255));
        draw_text(&e.status, 34., h - 10., 14., muted());
        if !e.quick.pending.is_empty() {
            draw_text(
                &format!("Input: {}", e.quick.pending),
                w - RIGHT - 180.,
                70.,
                18.,
                accent(),
            );
        }
        if let Tool::Template(name) = e.tool {
            draw_text(
                &format!(
                    "{name} · {}° · R rotate · A attach · Esc cancel",
                    e.template_rotation
                ),
                w - RIGHT - 405.,
                TOP + 28.,
                15.,
                accent(),
            );
        }
        draw_text(
            "Keys: 1/2/3 bonds | Shift+R ring | E erase | Esc select | Cmd+Z undo | Cmd+S save",
            w - 660.,
            h - 10.,
            12.,
            muted(),
        );
        let cmd = is_key_down(KeyCode::LeftSuper)
            || is_key_down(KeyCode::RightSuper)
            || is_key_down(KeyCode::LeftControl)
            || is_key_down(KeyCode::RightControl);
        if e.file_dialog.is_none()
            && e.pending_export.is_none()
            && e.pending_close.is_none()
            && e.pending_recovery.is_none()
            && !e.pending_quit
            && !e.periodic_open
            && !e.templates_open
            && e.text_target.is_none()
            && e.menu_open.is_none()
            && e.info_overlay.is_none()
            && !e.settings_open
            && !e.recent_open
            && !e.palette_open
        {
            if cmd
                && (is_key_down(KeyCode::LeftShift) || is_key_down(KeyCode::RightShift))
                && is_key_pressed(KeyCode::P)
            {
                e.palette_open = true;
                e.palette_query.clear();
                e.palette_selection = 0;
            }
            if cmd && is_key_pressed(KeyCode::Z) {
                if is_key_down(KeyCode::LeftShift) || is_key_down(KeyCode::RightShift) {
                    e.redo()
                } else {
                    e.undo()
                }
            }
            if cmd && is_key_pressed(KeyCode::S) {
                e.file_dialog = Some(true);
            }
            if cmd && is_key_pressed(KeyCode::O) {
                e.file_dialog = Some(false);
            }
            if cmd && is_key_pressed(KeyCode::L) {
                e.start_copy_smiles();
            }
            if cmd
                && (is_key_down(KeyCode::LeftShift) || is_key_down(KeyCode::RightShift))
                && is_key_pressed(KeyCode::K)
            {
                e.start_copy_smarts();
            }
            if cmd && is_key_pressed(KeyCode::E) {
                e.periodic_open = true;
            }
            if cmd && is_key_pressed(KeyCode::Y) {
                e.redo();
            }
            if cmd && is_key_pressed(KeyCode::A) {
                e.select_all();
            }
            if cmd && is_key_down(KeyCode::LeftShift) && is_key_pressed(KeyCode::D) {
                e.duplicate_selection();
            }
            if !cmd && is_key_pressed(KeyCode::Q) {
                e.rotate_selection(std::f32::consts::FRAC_PI_2);
            }
            if !cmd && matches!(e.tool, Tool::Template(_)) && is_key_pressed(KeyCode::R) {
                e.template_rotation = (e.template_rotation + 60) % 360;
                e.quick.pending.clear();
                e.status = format!(
                    "Template rotation: {}° · click to place",
                    e.template_rotation
                );
            }
            if !cmd && matches!(e.tool, Tool::Template(_)) && is_key_pressed(KeyCode::A) {
                e.toggle_template_attachment();
            }
            if cmd && is_key_pressed(KeyCode::N) {
                e.new_document();
            }
            if is_key_pressed(KeyCode::Delete) || is_key_pressed(KeyCode::Backspace) {
                e.delete_selection();
            }
            if is_key_pressed(KeyCode::Escape) {
                e.cancel_template_placement();
                e.tool = Tool::Select;
                e.pending_arrow = None;
                e.selected = None;
                e.selected_atoms.clear();
                e.selected_bond = None;
                e.selected_graphic = None;
                e.selected_bonds.clear();
                e.selected_graphics.clear();
                e.selection_box = None;
                e.quick.pending.clear();
            }
            if is_key_down(KeyCode::LeftShift) && is_key_pressed(KeyCode::R) {
                e.tool = Tool::Ring(6);
            }
            if !cmd && is_key_pressed(KeyCode::E) {
                e.tool = Tool::Erase;
            }
            if !cmd && !matches!(e.tool, Tool::Template(_)) {
                while let Some(ch) = get_char_pressed() {
                    for a in e.quick.push(ch, get_time()) {
                        e.apply_quick(a);
                    }
                }
                for a in e.quick.tick(get_time()) {
                    e.apply_quick(a);
                }
            } else {
                e.quick.pending.clear();
            }
            if canvas(p) && !e.panning {
                if is_mouse_button_pressed(MouseButton::Left) {
                    e.drag_start = Some(world_p);
                    e.moved = false;
                    match e.tool {
                        Tool::Select => {
                            e.select_at(world_p, is_key_down(KeyCode::LeftShift));
                        }
                        Tool::Erase => {
                            if let Some(i) = e.hit_atom(world_p) {
                                e.checkpoint();
                                e.remove_atom(i);
                                e.selected = None;
                                e.selected_atoms.clear();
                                e.selected_bond = None;
                                e.selected_graphic = None;
                            } else if let Some(i) = e.hit_bond(world_p) {
                                e.checkpoint();
                                e.mol.bonds.remove(i);
                                e.selected_bond = None;
                            } else if let Some(i) = e.hit_graphic(world_p) {
                                e.checkpoint();
                                e.graphics.remove(i);
                                e.selected_graphic = None;
                            }
                        }
                        Tool::Atom => {
                            e.checkpoint();
                            if let Some(i) = e.hit_atom(world_p) {
                                e.mol.atoms[i].element = e.element.clone();
                            } else {
                                let el = e.element.clone();
                                e.add_atom(world_p, &el);
                            }
                        }
                        Tool::Ring(n) => {
                            e.checkpoint();
                            e.ring(world_p, n);
                        }
                        Tool::Template(kind) => e.insert_template(kind, world_p),
                        Tool::Charge(n) => {
                            if let Some(i) = e.hit_atom(world_p) {
                                e.checkpoint();
                                e.mol.atoms[i].charge = (e.mol.atoms[i].charge + n).clamp(-9, 9);
                            }
                        }
                        Tool::StyledBond(_) | Tool::Bond(_) => {
                            if e.apply_bond_tool_to_hit_bond(world_p) {
                                e.drag_start = None;
                            }
                        }
                        Tool::Arrow | Tool::Electron(_) | Tool::Chain | Tool::Shape(_) => {}
                        Tool::Plus => {
                            e.checkpoint();
                            e.graphics.push(Graphic::Plus { pos: world_p });
                        }
                        Tool::Text => {
                            e.text_target = Some(TextTarget::New(world_p));
                            e.text_buffer.clear();
                        }
                    }
                }
                if is_mouse_button_down(MouseButton::Left) && e.tool == Tool::Select {
                    if let Some(start) = e.selection_box {
                        if world_p.distance(start) > 3. {
                            e.select_box(start, world_p);
                        }
                    } else if let Some(start) = e.drag_start {
                        if e.has_selection() && world_p.distance(start) > 2. {
                            if !e.moved {
                                e.checkpoint();
                                e.moved = true;
                            }
                            e.move_selection(world_p - start);
                            e.drag_start = Some(world_p);
                        }
                    }
                }
                if is_mouse_button_released(MouseButton::Left) {
                    if let Some(start) = e.selection_box.take() {
                        e.select_box(start, world_p);
                    }
                    if matches!(e.tool, Tool::Arrow | Tool::Electron(_)) {
                        e.finish_arrow(world_p);
                    }
                    if let (Some(start), Tool::Shape(kind)) = (e.drag_start, e.tool) {
                        if canvas(e.to_screen(start, w, h)) {
                            e.checkpoint();
                            e.graphics.push(Graphic::Shape {
                                kind,
                                start,
                                end: world_p,
                            });
                            e.status = format!("Drew {}", kind.code());
                        }
                    }
                    if let Some(start) = e.drag_start.filter(|_| {
                        matches!(e.tool, Tool::Bond(_) | Tool::StyledBond(_) | Tool::Chain)
                    }) {
                        if canvas(e.to_screen(start, w, h)) {
                            if e.tool == Tool::Chain {
                                e.checkpoint();
                                e.draw_chain(start, world_p);
                            } else {
                                e.checkpoint();
                                let a = e.hit_atom(start).unwrap_or_else(|| e.add_atom(start, "C"));
                                let origin = e.mol.atoms[a].pos;
                                let end = if world_p.distance(origin) < 14. {
                                    e.suggest_endpoint(a)
                                } else {
                                    world_p
                                };
                                let b = e.hit_atom(end).unwrap_or_else(|| e.add_atom(end, "C"));
                                match e.tool {
                                    Tool::Bond(order) => e.add_bond(a, b, order),
                                    Tool::StyledBond(style) => {
                                        e.add_styled_bond(a, b, style);
                                        if let Some(q) = &e.bond_query {
                                            if let Some(bond) =
                                                e.mol.bonds.iter_mut().find(|bond| {
                                                    (bond.a == a && bond.b == b)
                                                        || (bond.a == b && bond.b == a)
                                                })
                                            {
                                                bond.query = Some(q.clone());
                                            }
                                        }
                                    }
                                    _ => {}
                                }
                            }
                        }
                    }
                    e.drag_start = None;
                }
            } else if is_mouse_button_released(MouseButton::Left) {
                e.drag_start = None;
            }
        } else if e.palette_open {
            draw_rectangle(0., 0., w, h, Color::from_rgba(0, 0, 0, 95));
            let commands = all_commands();
            let query = e.palette_query.to_lowercase();
            let matches: Vec<&str> = commands
                .into_iter()
                .filter(|name| *name != "Command palette" && name.to_lowercase().contains(&query))
                .collect();
            if !matches.is_empty() {
                e.palette_selection = e.palette_selection.min(matches.len() - 1);
            } else {
                e.palette_selection = 0;
            }
            let width = 540.;
            let row_height = 37.;
            let shown = matches.len().min(9);
            let height = 112. + shown as f32 * row_height;
            let x = (w - width) / 2.;
            let y = (h - height) / 2.;
            round_rect(
                x + 4.,
                y + 6.,
                width,
                height,
                15.,
                Color::from_rgba(0, 0, 0, 34),
            );
            round_rect(x, y, width, height, 15., WHITE);
            draw_text("COMMAND PALETTE", x + 20., y + 26., 12., muted());
            draw_rectangle_lines(x + 16., y + 39., width - 32., 40., 1., accent());
            draw_text(&e.palette_query, x + 27., y + 66., 20., ink());
            while let Some(ch) = get_char_pressed() {
                if !ch.is_control() {
                    e.palette_query.push(ch);
                    e.palette_selection = 0;
                }
            }
            if is_key_pressed(KeyCode::Backspace) {
                e.palette_query.pop();
                e.palette_selection = 0;
            }
            if is_key_pressed(KeyCode::Up) {
                e.palette_selection = e.palette_selection.saturating_sub(1);
            }
            if is_key_pressed(KeyCode::Down) && !matches.is_empty() {
                e.palette_selection = (e.palette_selection + 1).min(matches.len() - 1);
            }
            if is_key_pressed(KeyCode::Escape) {
                e.palette_open = false;
            }
            let mut selected_command = None;
            for (index, command) in matches.iter().take(9).enumerate() {
                if button(
                    x + 14.,
                    y + 84. + index as f32 * row_height,
                    width - 28.,
                    33.,
                    command,
                    index == e.palette_selection,
                ) {
                    selected_command = Some(*command);
                }
            }
            if is_key_pressed(KeyCode::Enter) {
                selected_command = matches.get(e.palette_selection).copied();
            }
            if let Some(command) = selected_command {
                e.palette_open = false;
                e.menu_action(command);
            }
            if matches.is_empty() {
                draw_text("No matching commands", x + 22., y + 100., 15., muted());
            }
        } else if let Some(idx) = e.menu_open {
            let x = (304. + idx as f32 * 94.).min(w - 265.);
            let y = 38.;
            let items = menu_items(idx);
            draw_rectangle(x, y, 250., items.len() as f32 * 39. + 12., WHITE);
            draw_rectangle_lines(x, y, 250., items.len() as f32 * 39. + 12., 1., GRAY);
            for (j, item) in items.iter().enumerate() {
                if button(x + 6., y + 6. + j as f32 * 39., 238., 34., item, false) {
                    e.menu_open = None;
                    e.menu_action(item);
                }
            }
            if is_key_pressed(KeyCode::Escape) {
                e.menu_open = None;
            }
            if is_mouse_button_pressed(MouseButton::Left)
                && !(p.x >= x
                    && p.x <= x + 250.
                    && p.y >= y
                    && p.y <= y + items.len() as f32 * 39. + 12.)
                && p.y > 36.
            {
                e.menu_open = None;
            }
        } else if let Some((docs, _)) = e.pending_recovery.as_ref() {
            draw_rectangle(0., 0., w, h, Color::from_rgba(0, 0, 0, 110));
            let dx = (w - 520.) / 2.;
            let dy = (h - 230.) / 2.;
            round_rect(dx, dy, 520., 230., 16., card_surface());
            draw_text(
                "Recover autosaved workspace?",
                dx + 24.,
                dy + 42.,
                23.,
                ink(),
            );
            draw_text(
                &format!(
                    "{} document tab(s) were autosaved. Original files stay untouched.",
                    docs.len()
                ),
                dx + 24.,
                dy + 82.,
                16.,
                muted(),
            );
            if button(dx + 24., dy + 150., 145., 42., "Recover", true) {
                e.apply_recovery(true);
            }
            if button(dx + 184., dy + 150., 145., 42., "Discard recovery", false) {
                e.apply_recovery(false);
            }
            if is_key_pressed(KeyCode::Escape) {
                e.apply_recovery(false);
            }
        } else if let Some(index) = e.pending_close {
            draw_rectangle(0., 0., w, h, Color::from_rgba(0, 0, 0, 110));
            let dx = (w - 520.) / 2.;
            let dy = (h - 230.) / 2.;
            round_rect(dx, dy, 520., 230., 16., card_surface());
            let title = e
                .documents
                .get(index)
                .map(|doc| doc.title.as_str())
                .unwrap_or("document");
            draw_text(
                "Save changes before closing?",
                dx + 24.,
                dy + 42.,
                23.,
                ink(),
            );
            draw_text(
                &format!("{} has unsaved changes.", title),
                dx + 24.,
                dy + 82.,
                16.,
                muted(),
            );
            if button(dx + 24., dy + 150., 130., 42., "Save", true) {
                e.close_pending_with_save();
            }
            if button(dx + 169., dy + 150., 145., 42., "Discard", false) {
                e.close_document(index);
            }
            if button(dx + 330., dy + 150., 150., 42., "Cancel", false)
                || is_key_pressed(KeyCode::Escape)
            {
                e.pending_close = None;
            }
        } else if e.pending_quit {
            draw_rectangle(0., 0., w, h, Color::from_rgba(0, 0, 0, 110));
            let dx = (w - 540.) / 2.;
            let dy = (h - 245.) / 2.;
            round_rect(dx, dy, 540., 245., 16., card_surface());
            draw_text("Unsaved documents", dx + 24., dy + 42., 24., ink());
            draw_text(
                "Save named documents, keep editing, or discard all unsaved work.",
                dx + 24.,
                dy + 84.,
                15.,
                muted(),
            );
            if button(dx + 24., dy + 150., 155., 42., "Save before quit", true) {
                e.save_before_quit();
            }
            if button(dx + 194., dy + 150., 150., 42., "Discard & quit", false) {
                e.discard_and_quit();
            }
            if button(dx + 360., dy + 150., 150., 42., "Keep editing", false)
                || is_key_pressed(KeyCode::Escape)
            {
                e.pending_quit = false;
            }
        } else if e.acid_base_open {
            draw_rectangle(0., 0., w, h, Color::from_rgba(0, 0, 0, 120));
            let width = 680.;
            let height = 540.;
            let x = (w - width) / 2.;
            let y = (h - height) / 2.;
            round_rect(
                x + 5.,
                y + 7.,
                width,
                height,
                18.,
                Color::from_rgba(0, 0, 0, 35),
            );
            round_rect(x, y, width, height, 18., card_surface());
            round_rect(x, y, width, 58., 18., Color::from_rgba(232, 245, 236, 255));
            draw_rectangle(x, y + 40., width, 18., Color::from_rgba(232, 245, 236, 255));
            draw_text(
                &localized(e.language, "Chemistry A acid-base workbench"),
                x + 22.,
                y + 36.,
                24.,
                ink(),
            );
            let tools = [
                TextTarget::BufferPh,
                TextTarget::PkaFromBuffer,
                TextTarget::WeakAcidPh,
                TextTarget::WeakBasePh,
                TextTarget::StrongMixturePh,
                TextTarget::PolyproticSpeciation,
                TextTarget::BjerrumDiagram,
                TextTarget::PolyproticTitrationPh,
                TextTarget::PolyproticTitrationCurve,
            ];
            let labels = [
                "Buffer pH from pKa",
                "pKa from buffer measurements",
                "Weak-acid pH",
                "Weak-base pH",
                "Strong acid/base mixture pH",
                "Polyprotic speciation",
                "Bjerrum diagram data",
                "Polyprotic titration pH",
                "Polyprotic titration curve",
            ];
            UI_INPUT_ENABLED.with(|enabled| enabled.set(true));
            for (index, (target, label)) in tools.iter().zip(labels).enumerate() {
                let column = index / 5;
                let row = index % 5;
                if button(
                    x + 22. + column as f32 * 326.,
                    y + 74. + row as f32 * 49.,
                    310.,
                    41.,
                    &localized(e.language, label),
                    false,
                ) {
                    e.start_acid_base_input(*target);
                }
            }
            draw_text(
                &localized(e.language, "pKa is entered from your class data or a cited table. No structure-based pKa prediction is made."),
                x + 24.,
                y + 355.,
                13.,
                muted(),
            );
            draw_text(
                "Aqueous ideal-solution model · 25 °C · Kw = 1.0×10⁻¹⁴ · concentrations in mol/L",
                x + 24.,
                y + 381.,
                13.,
                muted(),
            );
            draw_text(
                &localized(e.language, "Choose a calculation, then enter values in the order shown. Use ; between fields and | between successive pKa values."),
                x + 24.,
                y + 407.,
                12.,
                muted(),
            );
            draw_text(
                &localized(e.language, "Henderson–Hasselbalch: pH = pKa + log₁₀([A⁻]/[HA]). Charts can copy data as CSV."),
                x + 24.,
                y + 429.,
                12.,
                muted(),
            );
            if button(x + width - 122., y + height - 55., 98., 36., "Close", true)
                || is_key_pressed(KeyCode::Escape)
            {
                e.acid_base_open = false;
            }
            UI_INPUT_ENABLED.with(|enabled| enabled.set(false));
        } else if e.settings_open {
            let width = 440.;
            let height = 400.;
            let x = (w - width) / 2.;
            let y = (h - height) / 2.;
            round_rect(
                x + 5.,
                y + 7.,
                width,
                height,
                18.,
                Color::from_rgba(0, 0, 0, 30),
            );
            round_rect(x, y, width, height, 18., WHITE);
            round_rect(x, y, width, 58., 18., Color::from_rgba(238, 247, 241, 255));
            draw_rectangle(x, y + 40., width, 18., Color::from_rgba(238, 247, 241, 255));
            draw_text(
                &localized(e.language, "Settings"),
                x + 22.,
                y + 36.,
                25.,
                ink(),
            );
            draw_text(
                &localized(e.language, "Interface language"),
                x + 24.,
                y + 94.,
                17.,
                DARKGRAY,
            );
            let english_selected = e.language == Language::English;
            let danish_selected = e.language == Language::Danish;
            UI_INPUT_ENABLED.with(|enabled| enabled.set(true));
            if button(x + 24., y + 113., 185., 48., "English", english_selected) {
                e.set_language(Language::English);
            }
            if button(x + 231., y + 113., 185., 48., "Dansk", danish_selected) {
                e.set_language(Language::Danish);
            }
            draw_text(
                &localized(e.language, "Appearance"),
                x + 24.,
                y + 196.,
                17.,
                ink(),
            );
            if button(
                x + 24.,
                y + 214.,
                392.,
                42.,
                if e.dark_theme {
                    "Dark theme · on"
                } else {
                    "Light theme · on"
                },
                e.dark_theme,
            ) {
                e.set_dark_theme(!e.dark_theme);
            }
            draw_text(
                &localized(e.language, "Implicit hydrogen labels"),
                x + 24.,
                y + 285.,
                17.,
                ink(),
            );
            if button(
                x + 24.,
                y + 303.,
                392.,
                42.,
                &localized(
                    e.language,
                    if e.show_implicit_hydrogens {
                        "Shown"
                    } else {
                        "Hidden"
                    },
                ),
                e.show_implicit_hydrogens,
            ) {
                e.set_implicit_hydrogen_visibility(!e.show_implicit_hydrogens);
            }
            draw_text(
                &localized(
                    e.language,
                    "Language and appearance are saved automatically.",
                ),
                x + 24.,
                y + 366.,
                14.,
                muted(),
            );
            if button(x + width - 122., y + height - 55., 98., 36., "Close", true)
                || is_key_pressed(KeyCode::Escape)
            {
                e.settings_open = false;
            }
            UI_INPUT_ENABLED.with(|enabled| enabled.set(false));
        } else if e.recent_open {
            draw_rectangle(0., 0., w, h, Color::from_rgba(0, 0, 0, 100));
            let dx = (w - 700.) / 2.;
            let dy = (h - 500.).max(20.) / 2.;
            round_rect(dx, dy, 700., 500., 16., WHITE);
            draw_text("Recent files", dx + 22., dy + 38., 25., ink());
            if e.recent_files.is_empty() {
                draw_text("No recent files yet", dx + 24., dy + 80., 18., muted());
            }
            let recent_paths = e.recent_files.clone();
            for (i, path) in recent_paths.iter().take(10).enumerate() {
                let label = Path::new(path)
                    .file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or(path);
                if button(
                    dx + 20.,
                    dy + 55. + i as f32 * 38.,
                    660.,
                    33.,
                    &label.chars().take(72).collect::<String>(),
                    false,
                ) {
                    let path = path.clone();
                    e.recent_open = false;
                    e.load(&path);
                }
            }
            if button(dx + 570., dy + 450., 100., 34., "Close", true)
                || is_key_pressed(KeyCode::Escape)
            {
                e.recent_open = false;
            }
        } else if let Some(info) = e.info_overlay.clone() {
            draw_rectangle(0., 0., w, h, Color::from_rgba(0, 0, 0, 100));
            let chart = e.acid_base_chart.clone();
            let mw = w.min(if chart.is_some() { 980. } else { 650. });
            let dx = (w - mw) / 2.;
            let lines = wrap_info_lines(&info, if chart.is_some() { 116 } else { 75 });
            let chart_top = 88.;
            let chart_bottom = if chart.is_some() { 405. } else { 0. };
            let visible =
                ((h - if chart.is_some() { 575. } else { 166. }).max(24.) / 24.).floor() as usize;
            let max_scroll = lines.len().saturating_sub(visible);
            e.info_scroll = e.info_scroll.min(max_scroll);
            let end = (e.info_scroll + visible).min(lines.len());
            let mh = ((if chart.is_some() { 455. } else { 106. })
                + end.saturating_sub(e.info_scroll) as f32 * 24.)
                .max(if chart.is_some() { 550. } else { 300. })
                .min(h - 60.);
            let dy = (h - mh) / 2.;
            round_rect(dx, dy, mw, mh, 16., WHITE);
            draw_text(
                chart
                    .as_ref()
                    .map_or("Analysis & information", |chart| chart.title.as_str()),
                dx + 20.,
                dy + 35.,
                24.,
                ink(),
            );
            if let Some(chart) = chart.as_ref() {
                let plot_x = dx + 74.;
                let plot_y = dy + chart_top;
                let plot_w = mw - 112.;
                let plot_h = chart_bottom - chart_top - 22.;
                let x_min = chart.points.first().map_or(0.0, |point| point.0);
                let mut x_max = chart.points.last().map_or(1.0, |point| point.0);
                if x_max <= x_min {
                    x_max = x_min + 1.0;
                }
                let (y_min, mut y_max) =
                    if chart.y_label == "Fraction (%)" || chart.y_label == "Andel (%)" {
                        (0.0, 100.0)
                    } else {
                        let values = chart.points.iter().flat_map(|point| point.1.iter());
                        let lo = values.clone().copied().fold(f64::INFINITY, f64::min);
                        let hi = values.copied().fold(f64::NEG_INFINITY, f64::max);
                        ((lo.floor() - 1.0).max(0.0), (hi.ceil() + 1.0).min(14.0))
                    };
                if y_max <= y_min {
                    y_max = y_min + 1.0;
                }
                let grid = Color::from_rgba(223, 232, 226, 255);
                let palette = [
                    Color::from_rgba(26, 132, 87, 255),
                    Color::from_rgba(54, 120, 190, 255),
                    Color::from_rgba(203, 111, 40, 255),
                    Color::from_rgba(143, 83, 174, 255),
                    Color::from_rgba(195, 65, 77, 255),
                    Color::from_rgba(44, 150, 157, 255),
                    Color::from_rgba(132, 143, 50, 255),
                    Color::from_rgba(95, 100, 110, 255),
                    Color::from_rgba(194, 126, 152, 255),
                ];
                for tick in 0..=5 {
                    let fraction = tick as f32 / 5.0;
                    let y = plot_y + plot_h * (1.0 - fraction);
                    draw_line(plot_x, y, plot_x + plot_w, y, 1.0, grid);
                    let label = y_min + (y_max - y_min) * f64::from(fraction);
                    draw_text(&format!("{label:.0}"), dx + 31., y + 4., 13., muted());
                    let x = plot_x + plot_w * fraction;
                    draw_line(x, plot_y, x, plot_y + plot_h, 1.0, grid);
                    let x_label = x_min + (x_max - x_min) * f64::from(fraction);
                    draw_text(
                        &format!("{x_label:.1}"),
                        x - 13.,
                        plot_y + plot_h + 18.,
                        12.,
                        muted(),
                    );
                }
                draw_line(plot_x, plot_y, plot_x, plot_y + plot_h, 1.5, ink());
                draw_line(
                    plot_x,
                    plot_y + plot_h,
                    plot_x + plot_w,
                    plot_y + plot_h,
                    1.5,
                    ink(),
                );
                for (series_index, label) in chart.series.iter().enumerate() {
                    let color = palette[series_index % palette.len()];
                    let mut previous = None;
                    for (x_value, values) in &chart.points {
                        let Some(y_value) = values.get(series_index).copied() else {
                            continue;
                        };
                        let x = plot_x + ((*x_value - x_min) / (x_max - x_min)) as f32 * plot_w;
                        let y =
                            plot_y + (1.0 - ((y_value - y_min) / (y_max - y_min)) as f32) * plot_h;
                        if let Some((px, py)) = previous {
                            draw_line(px, py, x, y, 2.4, color);
                        }
                        previous = Some((x, y));
                    }
                    let legend_x = plot_x + (series_index % 5) as f32 * (plot_w / 5.0);
                    let legend_y = dy + chart_bottom - 1. + (series_index / 5) as f32 * 18.;
                    draw_line(
                        legend_x,
                        legend_y - 3.,
                        legend_x + 20.,
                        legend_y - 3.,
                        3.,
                        color,
                    );
                    draw_text(label, legend_x + 25., legend_y + 1., 13., ink());
                }
                draw_text(&chart.y_label, dx + 19., plot_y + plot_h / 2., 13., muted());
                draw_text(
                    &chart.x_label,
                    plot_x + plot_w / 2. - 70.,
                    dy + chart_bottom + 34.,
                    13.,
                    muted(),
                );
                draw_text(
                    "IDEAL AQUEOUS MODEL · 25 °C",
                    dx + mw - 260.,
                    dy + 35.,
                    11.,
                    muted(),
                );
            }
            for (j, line) in lines[e.info_scroll..end].iter().enumerate() {
                let text_y = dy + if chart.is_some() { 435. } else { 68. } + j as f32 * 24.;
                draw_text(line, dx + 20., text_y, 15., DARKGRAY);
            }
            if lines.len() > visible {
                draw_text(
                    &format!("Scroll · {}–{} of {}", e.info_scroll + 1, end, lines.len()),
                    dx + 20.,
                    dy + mh - 29.,
                    12.,
                    muted(),
                );
            }
            if let Some(chart) = chart.as_ref() {
                if button(
                    dx + 20.,
                    dy + mh - 50.,
                    124.,
                    34.,
                    &localized(e.language, "Copy CSV data"),
                    false,
                ) {
                    let csv = chemistry_chart_csv(chart);
                    match copy_text_to_clipboard(&csv) {
                        Ok(()) => e.status = "Copied plotted Chemistry A data as CSV".into(),
                        Err(error) => e.status = error,
                    }
                }
            }
            if button(
                dx + mw - 120.,
                dy + mh - 50.,
                100.,
                34.,
                &localized(e.language, "Close"),
                true,
            ) || is_key_pressed(KeyCode::Escape)
            {
                e.info_overlay = None;
                e.acid_base_chart = None;
                e.naming_highlight.clear();
            }
        } else if e.text_target.is_some() {
            draw_rectangle(0., 0., w, h, Color::from_rgba(0, 0, 0, 100));
            let text_target = e.text_target.unwrap();
            let acid_base_prompt_text = acid_base_prompt(text_target);
            let acid_base_input = acid_base_prompt_text.is_some();
            let width = if acid_base_input { 650. } else { 500. };
            let height = if acid_base_input { 250. } else { 195. };
            let dx = (w - width) / 2.;
            let template_name = matches!(e.text_target, Some(TextTarget::TemplateName));
            let template_rename = matches!(e.text_target, Some(TextTarget::RenameUserTemplate(_)));
            let molecule_property = matches!(e.text_target, Some(TextTarget::MoleculeProperty));
            let reaction_coefficient =
                matches!(e.text_target, Some(TextTarget::ReactionCoefficient));
            let template_edit = template_name || template_rename;
            let dy = (h - height) / 2.;
            round_rect(dx, dy, width, height, 14., card_surface());
            let dialog_title = if acid_base_input {
                localized(e.language, "Chemistry A acid-base workbench").into_owned()
            } else if molecule_property {
                localized(e.language, "Edit molecule data").into_owned()
            } else if reaction_coefficient {
                localized(e.language, "Reaction coefficient: custom").into_owned()
            } else if template_name {
                "Name custom template".into()
            } else if template_rename {
                "Rename custom template".into()
            } else {
                "Enter text".into()
            };
            draw_text(&dialog_title, dx + 20., dy + 35., 24., DARKGRAY);
            draw_rectangle_lines(dx + 20., dy + 55., width - 40., 39., 1., accent());
            draw_text(&e.text_buffer, dx + 27., dy + 81., 20., DARKGRAY);
            if acid_base_input {
                let prompt = localized(e.language, acid_base_prompt_text.unwrap());
                draw_text(&prompt, dx + 21., dy + 116., 14., muted());
                draw_text(
                    &localized(
                        e.language,
                        "Enter values separated by semicolons. Decimal comma or point is accepted.",
                    ),
                    dx + 21.,
                    dy + 139.,
                    13.,
                    muted(),
                );
                draw_text(
                    &localized(e.language, acid_base_example(text_target).unwrap_or("")),
                    dx + 21.,
                    dy + 165.,
                    13.,
                    accent(),
                );
            } else if template_edit {
                draw_text(
                    "1–48 letters, numbers, spaces, hyphens or underscores",
                    dx + 21.,
                    dy + 111.,
                    13.,
                    muted(),
                );
            } else if molecule_property {
                draw_text(
                    &localized(
                        e.language,
                        "KEY=VALUE saves or updates · -KEY removes · NAME sets molecule name",
                    ),
                    dx + 21.,
                    dy + 111.,
                    12.,
                    muted(),
                );
            } else if reaction_coefficient {
                draw_text(
                    &localized(
                        e.language,
                        "Set a whole-number coefficient from 1 to 99 for this component.",
                    ),
                    dx + 21.,
                    dy + 111.,
                    13.,
                    muted(),
                );
            }
            while let Some(ch) = get_char_pressed() {
                if !ch.is_control() {
                    e.text_buffer.push(ch);
                }
            }
            if is_key_pressed(KeyCode::Backspace) {
                e.text_buffer.pop();
            }
            if is_key_pressed(KeyCode::Escape) {
                e.text_target = None;
                if template_edit {
                    e.templates_open = false;
                }
            }
            if is_key_pressed(KeyCode::Enter) {
                e.commit_text();
            }
            if button(
                dx + width - 220.,
                dy + height - 52.,
                100.,
                35.,
                "Cancel",
                false,
            ) {
                e.text_target = None;
                if template_edit {
                    e.templates_open = false;
                }
            }
            if button(
                dx + width - 108.,
                dy + height - 52.,
                108.,
                35.,
                if template_rename {
                    "Rename"
                } else if template_name || molecule_property {
                    "Save"
                } else {
                    "Apply"
                },
                true,
            ) {
                e.commit_text();
            }
        } else if e.periodic_open {
            draw_rectangle(0., 0., w, h, Color::from_rgba(0, 0, 0, 100));
            let dx = (w - 800.) / 2.;
            let dy = (h - 420.) / 2.;
            draw_rectangle(dx, dy, 800., 420., WHITE);
            draw_text(
                "Periodic table - choose an element",
                dx + 15.,
                dy + 31.,
                23.,
                DARKGRAY,
            );
            for (z, sym) in chemistry::SYMBOLS.iter().enumerate() {
                let (col, row) = periodic_cell(z);
                let x = dx + 17. + col as f32 * 42.;
                let y = dy + 48. + row as f32 * 37.;
                if button(x, y, 38., 32., sym, false) {
                    e.element = (*sym).into();
                    e.tool = Tool::Atom;
                    e.periodic_open = false;
                    e.status = format!("Element: {sym}");
                }
            }
            if is_key_pressed(KeyCode::Escape) {
                e.periodic_open = false;
            }
        } else if e.templates_open {
            draw_rectangle(0., 0., w, h, Color::from_rgba(0, 0, 0, 100));
            let dx = (w - 500.) / 2.;
            let dy = (h - 465.) / 2.;
            round_rect(dx, dy, 500., 465., 16., card_surface());
            draw_text("Templates", dx + 20., dy + 35., 24., DARKGRAY);
            draw_rectangle_lines(dx + 18., dy + 48., 464., 36., 1., accent());
            draw_text(
                &format!("Search: {}", e.template_query),
                dx + 28.,
                dy + 72.,
                17.,
                ink(),
            );
            if button(
                dx + 18.,
                dy + 90.,
                216.,
                30.,
                "Built-in",
                !e.template_tab_custom,
            ) {
                e.template_tab_custom = false;
                e.pending_delete_template = None;
                e.template_page = 0;
            }
            if button(
                dx + 246.,
                dy + 90.,
                236.,
                30.,
                "My templates",
                e.template_tab_custom,
            ) {
                e.template_tab_custom = true;
                e.pending_delete_template = None;
                e.template_page = 0;
            }
            while let Some(ch) = get_char_pressed() {
                if !ch.is_control() {
                    e.template_query.push(ch);
                    e.template_page = 0;
                }
            }
            if is_key_pressed(KeyCode::Backspace) {
                e.template_query.pop();
                e.template_page = 0;
            }
            if !e.template_tab_custom {
                let all_choices: Vec<&str> = [
                    "Benzene",
                    "Cyclohexane",
                    "Naphthalene",
                    "Decalin",
                    "Indole",
                    "Quinoline",
                    "Anthracene",
                    "Phenanthrene",
                    "Norbornane",
                    "Pyridine",
                    "Ethanol",
                    "Carboxyl",
                    "Cyclopropane",
                    "Cyclobutane",
                    "Cyclopentane",
                    "Cycloheptane",
                    "Cyclooctane",
                    "Pyrimidine",
                    "Imidazole",
                    "Pyrrole",
                    "Furan",
                    "Thiophene",
                    "Hydroxyl",
                    "Amino",
                    "Carbonyl",
                    "Aldehyde",
                    "Nitrile",
                    "Ester",
                    "Amide",
                    "Glycine",
                    "Alanine",
                    "Valine",
                    "Leucine",
                    "Isoleucine",
                    "Serine",
                    "Threonine",
                    "Cysteine",
                    "Methionine",
                    "Aspartic acid",
                    "Asparagine",
                    "Glutamic acid",
                    "Glutamine",
                    "Lysine",
                    "Arginine",
                    "Histidine",
                    "Phenylalanine",
                    "Tyrosine",
                    "Tryptophan",
                    "Proline",
                ]
                .into_iter()
                .filter(|name| builtin_template_matches_query(name, &e.template_query))
                .collect();
                let page_count = all_choices.len().div_ceil(6).max(1);
                e.template_page = e.template_page.min(page_count - 1);
                let choices: Vec<&str> = all_choices
                    .iter()
                    .skip(e.template_page * 6)
                    .take(6)
                    .copied()
                    .collect();
                let ring_names = [
                    "Benzene",
                    "Cyclohexane",
                    "Naphthalene",
                    "Decalin",
                    "Indole",
                    "Quinoline",
                    "Anthracene",
                    "Phenanthrene",
                    "Norbornane",
                    "Pyridine",
                    "Cyclopropane",
                    "Cyclobutane",
                    "Cyclopentane",
                    "Cycloheptane",
                    "Cyclooctane",
                    "Pyrimidine",
                    "Imidazole",
                    "Pyrrole",
                    "Furan",
                    "Thiophene",
                ];
                draw_text("Rings", dx + 20., dy + 135., 12., muted());
                draw_text("Groups · amino acids", dx + 265., dy + 135., 12., muted());
                let mut category_rows = [0usize; 2];
                for kind in &choices {
                    let col = if ring_names.contains(kind) { 0 } else { 1 };
                    let row = category_rows[col];
                    category_rows[col] += 1;
                    if button(
                        dx + 18. + col as f32 * 232.,
                        dy + 142. + row as f32 * 44.,
                        224.,
                        36.,
                        kind,
                        false,
                    ) {
                        e.activate_builtin_template(kind);
                        e.template_query.clear();
                    }
                }
                if all_choices.is_empty() {
                    draw_text("No matching templates", dx + 20., dy + 165., 15., muted());
                }
                if page_count > 1 {
                    if button(dx + 18., dy + 421., 96., 30., "← Previous", false)
                        && e.template_page > 0
                    {
                        e.template_page -= 1;
                    }
                    draw_text(
                        &format!("{} / {page_count}", e.template_page + 1),
                        dx + 218.,
                        dy + 442.,
                        14.,
                        muted(),
                    );
                    if button(dx + 386., dy + 421., 96., 30., "Next →", false)
                        && e.template_page + 1 < page_count
                    {
                        e.template_page += 1;
                    }
                }
            } else {
                let matching: Vec<(usize, String)> = e
                    .user_templates
                    .iter()
                    .enumerate()
                    .filter(|(_, template)| {
                        template
                            .name
                            .to_lowercase()
                            .contains(&e.template_query.to_lowercase())
                    })
                    .map(|(index, template)| (index, template.name.clone()))
                    .collect();
                let page_count = matching.len().div_ceil(6).max(1);
                e.template_page = e.template_page.min(page_count - 1);
                for (row, (index, name)) in matching
                    .iter()
                    .skip(e.template_page * 6)
                    .take(6)
                    .enumerate()
                {
                    let y = dy + 132. + row as f32 * 44.;
                    if button(
                        dx + 18.,
                        y,
                        194.,
                        36.,
                        &name.chars().take(28).collect::<String>(),
                        false,
                    ) {
                        e.activate_custom_template(e.user_templates[*index].mol.clone());
                        e.template_query.clear();
                    }
                    if button(dx + 218., y, 60., 36., "Insert", false) {
                        e.activate_custom_template(e.user_templates[*index].mol.clone());
                        e.template_query.clear();
                    }
                    if button(dx + 284., y, 76., 36., "Rename", false) {
                        e.text_buffer = name.clone();
                        e.text_target = Some(TextTarget::RenameUserTemplate(*index));
                    }
                    let deleting = e.pending_delete_template == Some(*index);
                    if button(
                        dx + 366.,
                        y,
                        116.,
                        36.,
                        if deleting { "Confirm delete" } else { "Delete" },
                        deleting,
                    ) {
                        if deleting {
                            match e.delete_user_template(*index) {
                                Ok(()) => e.pending_delete_template = None,
                                Err(error) => e.status = error,
                            }
                        } else {
                            e.pending_delete_template = Some(*index);
                        }
                    }
                }
                if matching.is_empty() {
                    draw_text("No custom templates yet", dx + 20., dy + 160., 15., muted());
                }
                draw_text(
                    &format!("Page {} / {page_count}", e.template_page + 1),
                    dx + 300.,
                    dy + 438.,
                    13.,
                    muted(),
                );
                if e.template_page > 0 && button(dx + 270., dy + 416., 72., 31., "Previous", false)
                {
                    e.template_page -= 1;
                }
                if e.template_page + 1 < page_count
                    && button(dx + 350., dy + 416., 72., 31., "Next", false)
                {
                    e.template_page += 1;
                }
            }
            if button(
                dx + 18.,
                dy + 416.,
                235.,
                31.,
                "Save selected fragment",
                false,
            ) {
                if e.selected.is_none() && e.selected_atoms.is_empty() {
                    e.status = "Select a connected fragment before saving a template".into();
                } else {
                    e.text_buffer.clear();
                    e.text_target = Some(TextTarget::TemplateName);
                }
            }
            if is_key_pressed(KeyCode::Escape) {
                e.templates_open = false;
                e.template_query.clear();
                e.pending_delete_template = None;
            }
        } else if e.pending_export.is_some() {
            draw_rectangle(0., 0., w, h, Color::from_rgba(0, 0, 0, 90));
            let dx = (w - 560.) / 2.;
            let dy = (h - 190.) / 2.;
            draw_rectangle(dx, dy, 560., 190., WHITE);
            draw_text(
                "Export may lose information",
                dx + 20.,
                dy + 35.,
                24.,
                DARKGRAY,
            );
            draw_text(
                &e.status.chars().take(65).collect::<String>(),
                dx + 20.,
                dy + 75.,
                17.,
                RED,
            );
            draw_text(
                "Save as .mols to preserve the full document.",
                dx + 20.,
                dy + 105.,
                17.,
                DARKGRAY,
            );
            if button(dx + 250., dy + 135., 135., 34., "Cancel", false)
                || is_key_pressed(KeyCode::Escape)
            {
                e.pending_export = None;
            }
            if button(dx + 397., dy + 135., 140., 34., "Export", true)
                || is_key_pressed(KeyCode::Enter)
            {
                if let Some(path) = e.pending_export.take() {
                    e.save_confirmed(&path);
                }
            }
        } else {
            draw_rectangle(0., 0., w, h, Color::from_rgba(0, 0, 0, 90));
            let dx = (w - 760.) / 2.;
            let dy = (h - 190.) / 2.;
            draw_rectangle(dx, dy, 760., 190., WHITE);
            draw_rectangle_lines(dx, dy, 760., 190., 1., GRAY);
            let saving = e.file_dialog == Some(true);
            draw_text(
                if saving {
                    "Save molecule"
                } else {
                    "Open MOL, SDF, RXN, SMILES, SMARTS or MOLS"
                },
                dx + 22.,
                dy + 38.,
                25.,
                DARKGRAY,
            );
            draw_rectangle_lines(dx + 20., dy + 58., 720., 39., 1., accent());
            draw_text(&e.path_input, dx + 29., dy + 84., 20., DARKGRAY);
            while let Some(ch) = get_char_pressed() {
                if !ch.is_control() {
                    e.path_input.push(ch);
                }
            }
            if is_key_pressed(KeyCode::Backspace) {
                e.path_input.pop();
            }
            if is_key_pressed(KeyCode::Escape) {
                e.file_dialog = None;
                e.close_after_save = None;
                e.quit_when_clean = false;
            }
            if is_key_pressed(KeyCode::Enter) {
                let path = e.path_input.clone();
                if saving {
                    e.request_save(path);
                } else {
                    e.load(&path);
                }
                e.file_dialog = None;
            }
            if button(dx + 530., dy + 125., 94., 36., "Cancel", false) {
                e.file_dialog = None;
                e.close_after_save = None;
                e.quit_when_clean = false;
            }
            if button(
                dx + 636.,
                dy + 125.,
                104.,
                36.,
                if saving { "Save" } else { "Open" },
                true,
            ) {
                let path = e.path_input.clone();
                if saving {
                    e.request_save(path);
                } else {
                    e.load(&path);
                }
                e.file_dialog = None;
            }
            draw_text(
                "Use an absolute path or a path relative to the launch folder.",
                dx + 20.,
                dy + 117.,
                15.,
                GRAY,
            );
        }
        if e.file_dialog.is_none()
            && e.pending_export.is_none()
            && e.pending_close.is_none()
            && e.pending_recovery.is_none()
            && !e.pending_quit
            && !e.periodic_open
            && !e.templates_open
            && e.text_target.is_none()
            && e.info_overlay.is_none()
            && e.menu_open.is_none()
        {
            if let Some((i, (_, t))) = tools.iter().enumerate().find(|(i, _)| {
                let y = TOP + 30. + *i as f32 * 35.;
                p.x >= 18. && p.x < 70. && p.y >= y && p.y < y + 33.
            }) {
                let text = tool_help(*t);
                let ty = (TOP + 30. + i as f32 * 35.).min(h - BOTTOM - 40.);
                round_rect(LEFT + 8., ty, 370., 34., 9., navy());
                draw_text(text, LEFT + 18., ty + 23., 14., WHITE);
            }
        }
        e.persist_zoom_preference();
        e.sync_active_document();
        next_frame().await;
    }
}
fn tool_help(t: Tool) -> &'static str {
    match t {
        Tool::Select => "Select and drag atoms or graphic objects (Esc)",
        Tool::Erase => "Delete one atom, bond or object (E)",
        Tool::Bond(1) => "Draw a single bond (1)",
        Tool::Bond(2) => "Draw or change a bond to a double bond (2)",
        Tool::Bond(3) => "Draw a triple bond (3)",
        Tool::StyledBond(BondStyle::Aromatic) => "Draw an aromatic bond (4)",
        Tool::StyledBond(BondStyle::Wedge) => "Draw a filled wedge (5)",
        Tool::StyledBond(BondStyle::Hash) => "Draw a hashed wedge (6)",
        Tool::StyledBond(BondStyle::Wavy) => "Draw an unspecified stereo bond (7)",
        Tool::Atom => "Place or replace an atom (element symbol)",
        Tool::Chain => "Draw a zigzag carbon chain (drag between endpoints)",
        Tool::Ring(6) => "Place a six-membered ring (Shift+R)",
        Tool::Ring(5) => "Place a five-membered ring",
        Tool::Charge(1) => "Increase atom charge (+)",
        Tool::Charge(-1) => "Decrease atom charge (-)",
        Tool::Arrow => "Drag to draw a reaction arrow",
        Tool::Electron(1) => "Drag to draw a single-electron fishhook arrow",
        Tool::Electron(_) => "Drag to draw a two-electron curved arrow",
        Tool::Plus => "Place a reaction plus sign",
        Tool::Text => "Place a text label",
        Tool::Shape(ShapeKind::Line) => "Drag to draw a line",
        Tool::Shape(ShapeKind::Rectangle) => "Drag to draw a rectangle",
        Tool::Shape(ShapeKind::Ellipse) => "Drag to draw an ellipse",
        Tool::Template(_) => "R rotates by 60°; A attaches to one selected atom; Esc cancels",
        _ => "Draw a bond",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chemistry_a_workbench_routes_inputs_to_real_pka_and_ph_calculations() {
        let mut editor = Editor::new();
        editor.menu_action("Chemistry A acid/base tools");
        assert!(editor.acid_base_open);

        editor.start_acid_base_input(TextTarget::BufferPh);
        editor.text_buffer = "4,76;0,10;0,20".into();
        editor.commit_text();
        assert!(editor
            .info_overlay
            .as_deref()
            .is_some_and(|report| report.contains("Buffer pH = 5.0610")));

        editor.start_acid_base_input(TextTarget::PkaFromBuffer);
        editor.text_buffer = "5,06;0,10;0,20".into();
        editor.commit_text();
        assert!(editor
            .info_overlay
            .as_deref()
            .is_some_and(|report| report.contains("pKa = 4.7590")));

        editor.start_acid_base_input(TextTarget::BjerrumDiagram);
        editor.text_buffer = "2,15|7,20|12,35;0;14;2".into();
        editor.commit_text();
        assert!(editor.info_overlay.as_deref().is_some_and(|report| {
            report.contains("Bjerrum distribution table") && report.contains("14.00")
        }));
        let chart = editor.acid_base_chart.as_ref().unwrap();
        assert_eq!(chart.series.len(), 4);
        assert_eq!(chart.points.len(), 8);
        assert!(chart
            .points
            .iter()
            .all(|(_, fractions)| { (fractions.iter().sum::<f64>() - 100.0).abs() < 1e-9 }));

        editor.start_acid_base_input(TextTarget::PolyproticTitrationPh);
        editor.text_buffer = "4,76;25;0,10;0,10;12,5".into();
        editor.commit_text();
        assert!(editor
            .info_overlay
            .as_deref()
            .is_some_and(|report| report.contains("Titration pH =")));
    }

    #[test]
    fn chemistry_a_workbench_localizes_calculation_results_and_keeps_plottable_data() {
        let mut editor = Editor::new();
        editor.language = Language::Danish;
        editor.calculate_acid_base(TextTarget::BufferPh, "4.76;0.1;0.2");
        assert!(editor
            .info_overlay
            .as_deref()
            .is_some_and(|report| report.contains("Buffer-pH = 5.0610")));

        editor.calculate_acid_base(TextTarget::PolyproticTitrationCurve, "4.76;25;0.1;0.1;1.5");
        let curve = editor.acid_base_chart.as_ref().unwrap();
        assert_eq!(curve.points.len(), 41);
        assert_eq!(curve.series, ["pH"]);
        assert!(curve
            .points
            .windows(2)
            .all(|points| { points[0].0 < points[1].0 && points[0].1[0] <= points[1].1[0] }));
        let csv = chemistry_chart_csv(curve);
        assert!(csv.starts_with("Tilsat stærk base (mL),pH\n"));
        assert_eq!(csv.lines().count(), 42);

        editor.calculate_acid_base(TextTarget::WeakAcidPh, "4.76;0.1");
        assert!(editor
            .info_overlay
            .as_deref()
            .is_some_and(|report| report.contains("Metode: ideel vandig monoprotisk ligevægt")));
    }

    #[test]
    fn chemistry_a_calculator_refuses_missing_and_nonphysical_inputs_without_fake_values() {
        let mut editor = Editor::new();
        editor.calculate_acid_base(TextTarget::WeakAcidPh, "4.76;0");
        assert!(editor.info_overlay.as_deref().is_some_and(|report| {
            report.contains("calculation unavailable") && report.contains("positive number")
        }));
        editor.calculate_acid_base(TextTarget::PkaFromBuffer, "7;0.1;0");
        assert!(editor.info_overlay.as_deref().is_some_and(|report| {
            report.contains("calculation unavailable") && report.contains("positive number")
        }));
    }

    #[test]
    fn three_d_atom_symbol_visibility_persists_without_clobbering_other_settings() {
        let path = std::env::temp_dir().join(format!(
            "molecula-3d-settings-{}-{}.conf",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        fs::write(&path, "language=da\ntheme=dark\n").unwrap();
        assert!(load_3d_atom_symbols_visibility_from(&path));
        save_3d_atom_symbols_visibility_to(&path, false).unwrap();
        assert!(!load_3d_atom_symbols_visibility_from(&path));
        let settings = fs::read_to_string(&path).unwrap();
        assert!(settings.contains("language=da"));
        assert!(settings.contains("theme=dark"));
        save_3d_atom_symbols_visibility_to(&path, true).unwrap();
        assert!(load_3d_atom_symbols_visibility_from(&path));
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn three_d_render_modes_provide_cpk_spacefill_and_monochrome_behavior() {
        assert_eq!(THREE_D_STYLE_COUNT, 5);
        assert_eq!(three_d_style_name(3, Language::English), "Spacefill");
        assert_eq!(three_d_style_name(3, Language::Danish), "Rumfyldning");
        assert_eq!(three_d_atom_color(2, "O"), element_color("O"));
        assert_ne!(element_color("O"), element_color("N"));
        assert_eq!(three_d_atom_color(4, "O"), three_d_atom_color(4, "C"));
        assert!(three_d_atom_radius(3, "Br", 70.) > three_d_atom_radius(3, "H", 70.));
        assert_eq!(three_d_atom_radius(0, "C", 70.), 2.5);
    }

    #[test]
    fn three_d_camera_rotation_preserves_coordinate_distance() {
        let point = vec3(1.25, -3.5, 4.75);
        let rotated = rotate_3d(point, vec2(0.83, -0.41));
        assert!((rotated.length() - point.length()).abs() < 0.0001);
        assert_eq!(
            project_3d(point, point, vec2(1.0, 1.0), 80., vec2(5., 9.)),
            vec2(5., 9.)
        );
    }

    #[test]
    fn background_calculation_errors_keep_pubchem_and_rdkit_distinct() {
        let (pubchem_notice, pubchem_status) =
            background_calculation_error("PUBCHEM_LOOKUP\ncurl: could not resolve host");
        assert!(pubchem_notice.contains("PubChem lookup failed"));
        assert!(pubchem_notice.contains("The drawing was not changed"));
        assert!(pubchem_notice.contains("could not resolve host"));
        assert_eq!(pubchem_status, "PubChem lookup failed");

        let (rdkit_notice, rdkit_status) = background_calculation_error("No module named rdkit");
        assert!(rdkit_notice.contains("RDKit unavailable"));
        assert!(rdkit_notice.contains("MOLECULA_PYTHON"));
        assert_eq!(rdkit_status, "Calculation failed");
    }

    #[test]
    fn molecule_data_dialog_updates_removes_and_undoes_sdf_fields() {
        let mut editor = Editor::new();
        editor.menu_action("Edit molecule data");
        assert!(matches!(
            editor.text_target,
            Some(TextTarget::MoleculeProperty)
        ));
        editor.text_buffer = "sample_id=RX-17".into();
        editor.commit_text();
        assert_eq!(
            editor.mol.properties,
            vec![("sample_id".into(), "RX-17".into())]
        );
        editor.text_buffer = "NAME=Ethyl acetate".into();
        editor.text_target = Some(TextTarget::MoleculeProperty);
        editor.commit_text();
        assert_eq!(editor.mol.name, "Ethyl acetate");
        assert!(editor
            .mol
            .properties
            .iter()
            .any(|(key, value)| key == "NAME" && value == "Ethyl acetate"));

        let mut sdf = molfile_for(&editor.mol);
        for (key, value) in &editor.mol.properties {
            sdf += &format!(">  <{key}>\n{value}\n\n");
        }
        sdf += "$$$$\n";
        let loaded = parse_sdf(&sdf).unwrap().remove(0);
        assert_eq!(loaded.properties, editor.mol.properties);
        let snapshot = Snapshot {
            mol: editor.mol.clone(),
            graphics: vec![],
        };
        assert_eq!(
            native::read(&native::write(&snapshot)).unwrap().mol,
            snapshot.mol
        );

        editor.text_buffer = "-sample_id".into();
        editor.text_target = Some(TextTarget::MoleculeProperty);
        editor.commit_text();
        assert!(!editor
            .mol
            .properties
            .iter()
            .any(|(key, _)| key == "sample_id"));
        editor.undo();
        assert!(editor
            .mol
            .properties
            .iter()
            .any(|(key, value)| key == "sample_id" && value == "RX-17"));

        let before = editor.mol.properties.clone();
        editor.text_buffer = "bad key=value".into();
        editor.text_target = Some(TextTarget::MoleculeProperty);
        editor.commit_text();
        assert_eq!(editor.mol.properties, before);
    }

    #[test]
    fn rdkit_cip_labels_apply_to_atoms_and_undo_without_touching_drawn_stereo() {
        let mut editor = Editor::new();
        for _ in 0..3 {
            editor.add_atom(vec2(0., 0.), "C");
        }
        editor.mol.atoms[0].stereo = Some("R".into());
        editor.mol.atoms[1].stereo = Some("@@".into());
        editor.mol.atoms[2].stereo = Some("S".into());
        editor.apply_stereo_labels(
            "ENGINE\tRDKit 2026\nMETHOD\tRDKit stereochemistry assignment\nATOM\t2\tS\n",
        );
        assert_eq!(editor.mol.atoms[0].stereo, None);
        assert_eq!(editor.mol.atoms[1].stereo.as_deref(), Some("S"));
        assert_eq!(editor.mol.atoms[2].stereo, None);
        editor.undo();
        assert_eq!(editor.mol.atoms[0].stereo.as_deref(), Some("R"));
        assert_eq!(editor.mol.atoms[1].stereo.as_deref(), Some("@@"));
        assert_eq!(editor.mol.atoms[2].stereo.as_deref(), Some("S"));

        let before = editor.mol.clone();
        editor.apply_stereo_labels("ATOM\t0\tR\n");
        assert_eq!(editor.mol, before);
        assert!(editor
            .info_overlay
            .as_deref()
            .is_some_and(|text| text.contains("invalid")));
        editor.apply_stereo_labels("ATOM\t1\tR\nATOM\t1\tS\n");
        assert_eq!(editor.mol, before);
    }

    #[test]
    fn unique_atom_map_assignment_respects_selection_and_is_undoable() {
        let mut editor = Editor::new();
        for _ in 0..4 {
            editor.add_atom(vec2(0., 0.), "C");
        }
        editor.mol.atoms[0].atom_map = Some(3);
        editor.selected_atoms = vec![2, 1];
        editor.assign_unique_atom_maps();
        assert_eq!(editor.mol.atoms[0].atom_map, Some(3));
        assert_eq!(editor.mol.atoms[1].atom_map, Some(1));
        assert_eq!(editor.mol.atoms[2].atom_map, Some(2));
        assert_eq!(editor.mol.atoms[3].atom_map, None);
        editor.undo();
        assert_eq!(editor.mol.atoms[1].atom_map, None);
        assert_eq!(editor.mol.atoms[2].atom_map, None);

        editor.clear_selection();
        editor.assign_unique_atom_maps();
        assert_eq!(
            editor
                .mol
                .atoms
                .iter()
                .filter_map(|atom| atom.atom_map)
                .collect::<Vec<_>>(),
            vec![1, 2, 3, 4]
        );
    }

    #[test]
    fn copy_smiles_command_refuses_unrepresentable_query_data_and_has_platform_clipboard() {
        assert!(all_commands().contains(&"Copy as SMILES"));
        assert!(all_commands().contains(&"Copy as SMARTS"));
        assert!(all_commands().contains(&"Copy drawing as image"));
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "C");
        editor.mol.atoms[0].alias = Some("Me".into());
        editor.menu_action("Copy as SMILES");
        assert!(editor.calc_rx.is_none());
        assert!(editor
            .info_overlay
            .as_deref()
            .is_some_and(|text| text.contains("cannot be represented safely")));
        #[cfg(target_os = "macos")]
        assert_eq!(clipboard_candidates()[0].0, "pbcopy");
        #[cfg(target_os = "windows")]
        assert_eq!(clipboard_candidates()[0].0, "clip.exe");
        #[cfg(all(unix, not(target_os = "macos")))]
        assert!(clipboard_candidates()
            .iter()
            .any(|(program, _)| *program == "wl-copy"));
    }

    #[test]
    fn copy_smarts_refuses_unrepresentable_aliases_before_touching_clipboard() {
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "C");
        editor.mol.atoms[0].alias = Some("Ph".into());
        editor.menu_action("Copy as SMARTS");
        assert!(editor.calc_rx.is_none());
        assert!(editor
            .info_overlay
            .as_deref()
            .is_some_and(|text| text.contains("cannot preserve aliases")));
    }

    #[test]
    fn smarts_file_export_writes_query_text_and_warns_about_document_loss() {
        let mut editor = Editor::new();
        let carbon = editor.add_atom(vec2(0., 0.), "C");
        let oxygen = editor.add_atom(vec2(40., 0.), "O");
        editor.add_bond(carbon, oxygen, 1);
        let path =
            std::env::temp_dir().join(format!("molecula-query-{}.smarts", std::process::id()));
        let path_text = path.to_string_lossy().to_string();
        editor.save(&path_text);
        assert_eq!(editor.pending_export.as_deref(), Some(path_text.as_str()));
        editor.save_confirmed(&path_text);
        assert_eq!(fs::read_to_string(&path).unwrap(), "C-O\n");
        assert!(editor.status.starts_with("Exported "));
        assert!(editor
            .export_losses(&path_text)
            .contains(&"document metadata, reaction annotation and drawing objects"));
        let mut imported = Editor::new();
        imported.load(&path_text);
        assert_eq!(imported.mol.atoms.len(), 2);
        assert_eq!(imported.mol.bonds.len(), 1);
        assert!(imported.status.starts_with("Opened SMARTS pattern"));
        let _ = fs::remove_file(path);
    }

    #[test]
    fn local_smiles_file_import_and_export_work_without_rdkit() {
        let path =
            std::env::temp_dir().join(format!("molecula-local-smiles-{}.smi", std::process::id()));
        let mut editor = Editor::new();
        editor.mol = smiles::parse("CC(=O)O").unwrap();
        editor.save(path.to_str().unwrap());
        assert!(editor.pending_export.is_some());
        editor.save_confirmed(path.to_str().unwrap());
        let text = fs::read_to_string(&path).unwrap();
        assert_eq!(text.trim(), "CC(=O)O");
        fs::write(&path, "\n# generated structure\nCC(=O)O acetic-acid\n").unwrap();

        let mut restored = Editor::new();
        restored.load(path.to_str().unwrap());
        assert_eq!(restored.mol.atoms.len(), 4);
        assert_eq!(restored.mol.bonds.len(), 3);
        assert_eq!(restored.status, format!("Opened SMILES {}", path.display()));
        let _ = fs::remove_file(path);
    }

    #[test]
    fn valence_override_control_cycles_back_to_automatic() {
        assert_eq!(cycle_valence_override(None), Some(1));
        assert_eq!(cycle_valence_override(Some(1)), Some(2));
        assert_eq!(cycle_valence_override(Some(14)), Some(0));
        assert_eq!(cycle_valence_override(Some(0)), None);
        assert_eq!(cycle_valence_override(Some(255)), Some(0));
    }

    #[test]
    fn v3000_roundtrip_handles_large_graphs_and_supported_atom_properties() {
        let mut mol = Molecule::default();
        for i in 0..1001 {
            let mut atom = Atom::new(vec2(i as f32, 0.), "C");
            if i == 0 {
                atom.charge = -1;
                atom.isotope = Some(13);
                atom.radical = 2;
                atom.atom_map = Some(777);
                atom.stereo = Some("either".into());
            }
            mol.atoms.push(atom);
        }
        for i in 0..1000 {
            mol.bonds.push(Bond {
                a: i,
                b: i + 1,
                order: 1,
                style: BondStyle::Normal,
                stereo: None,
                query: None,
            });
        }
        mol.bonds[0].style = BondStyle::Wedge;
        mol.bonds[1].order = 4;
        mol.bonds[1].style = BondStyle::Aromatic;
        let text = molfile_for(&mol);
        assert!(text.contains("V3000"));
        let got = parse_mol_at(&text, vec2(0., 0.)).unwrap();
        assert_eq!(got.atoms.len(), 1001);
        assert_eq!(got.bonds.len(), 1000);
        assert_eq!(got.atoms[0].charge, -1);
        assert_eq!(got.atoms[0].isotope, Some(13));
        assert_eq!(got.atoms[0].radical, 2);
        assert_eq!(got.atoms[0].atom_map, Some(777));
        assert_eq!(got.atoms[0].stereo.as_deref(), Some("either"));
        assert!(got.atoms[1].aromatic && got.atoms[2].aromatic);
        assert_eq!(got.bonds[0].style, BondStyle::Wedge);
    }

    #[test]
    fn mol_v2000_roundtrip_keeps_unspecified_atom_stereo() {
        let mut atom = Atom::new(vec2(0., 0.), "C");
        atom.stereo = Some("either".into());
        let mol = Molecule {
            atoms: vec![atom],
            ..Molecule::default()
        };
        let text = molfile_for(&mol);
        let got = parse_mol_at(&text, vec2(0., 0.)).unwrap();
        assert_eq!(got.atoms[0].stereo.as_deref(), Some("either"));
    }

    #[test]
    fn mol_v2000_roundtrip_preserves_nondefault_and_zero_valence() {
        for override_value in [0, 1, 5, 14] {
            let mut atom = Atom::new(vec2(250., 250.), "P");
            atom.valence_override = Some(override_value);
            let mol = Molecule {
                atoms: vec![atom],
                ..Molecule::default()
            };
            let text = molfile_for(&mol);
            let line = text.lines().nth(4).unwrap();
            let encoded = line
                .split_whitespace()
                .nth(10)
                .unwrap()
                .parse::<u8>()
                .unwrap();
            assert_eq!(
                encoded,
                if override_value == 0 {
                    15
                } else {
                    override_value
                }
            );
            let parsed = parse_mol_at(&text, vec2(0., 0.)).unwrap();
            assert_eq!(parsed.atoms[0].valence_override, Some(override_value));
        }
    }

    #[test]
    fn v3000_roundtrip_preserves_zero_valence_override() {
        let mut mol = Molecule::default();
        for i in 0..1000 {
            let mut atom = Atom::new(vec2(i as f32, 0.), "C");
            if i == 0 {
                atom.valence_override = Some(0);
            }
            mol.atoms.push(atom);
        }
        let text = molfile_for(&mol);
        assert!(text.contains("VAL=-1"));
        let parsed = parse_mol_at(&text, vec2(0., 0.)).unwrap();
        assert_eq!(parsed.atoms[0].valence_override, Some(0));
    }

    #[test]
    fn v3000_roundtrip_preserves_hydrogen_count_queries() {
        let mut mol = Molecule::default();
        for i in 0..1000 {
            let mut atom = Atom::new(vec2(i as f32, 0.), "C");
            if i == 0 {
                atom.implicit_h_override = Some(0);
            } else if i == 1 {
                atom.implicit_h_override = Some(5);
            }
            mol.atoms.push(atom);
        }
        let text = molfile_for(&mol);
        assert!(text.contains("HCOUNT=-1"));
        assert!(text.contains("HCOUNT=5"));
        let parsed = parse_mol_at(&text, vec2(0., 0.)).unwrap();
        assert_eq!(parsed.atoms[0].implicit_h_override, Some(0));
        assert_eq!(parsed.atoms[1].implicit_h_override, Some(5));

        let mut editor = Editor::new();
        editor.mol = mol;
        assert!(!editor
            .export_losses("query.mol")
            .contains(&"advanced atom properties"));
        editor.mol.atoms[0].implicit_h_override = Some(6);
        assert!(editor
            .export_losses("query.mol")
            .contains(&"advanced atom properties"));
    }

    #[test]
    fn v3000_hydrogen_count_parser_handles_unspecified_and_invalid_values() {
        let v3000 = |attribute: &str| {
            format!("H\n  Molecula\n\n  0  0  0     0  0            V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 {attribute}\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n")
        };
        assert_eq!(
            parse_mol_at(&v3000("HCOUNT=0"), vec2(0., 0.))
                .unwrap()
                .atoms[0]
                .implicit_h_override,
            None
        );
        assert!(parse_mol_at(&v3000("HCOUNT=6"), vec2(0., 0.)).is_err());
        assert!(parse_mol_at(&v3000("HCOUNT=oops"), vec2(0., 0.)).is_err());
    }

    #[test]
    fn v3000_atom_queries_roundtrip_substitution_unsaturation_and_ring_bonds() {
        let source = "Q\n  Molecula\n\n  0  0  0     0  0            V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 SUBST=3 UNSAT=1 RBCNT=-1\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
        let parsed = parse_mol_at(source, vec2(0., 0.)).unwrap();
        assert_eq!(parsed.atoms[0].query.as_deref(), Some(".s3;.u;.rb0"));

        // A supported V3000 query promotes a small graph automatically, so
        // exporting the structure does not silently downgrade it to V2000.
        let exported = molfile_for(&parsed);
        assert!(exported.contains("V3000"));
        assert!(exported.contains("SUBST=3"));
        assert!(exported.contains("UNSAT=1"));
        assert!(exported.contains("RBCNT=-1"));
        let roundtrip = parse_mol_at(&exported, vec2(0., 0.)).unwrap();
        assert_eq!(roundtrip.atoms[0].query.as_deref(), Some(".rb0;.s3;.u"));

        let mut editor = Editor::new();
        editor.mol = parsed;
        assert!(!editor
            .export_losses("query.mol")
            .contains(&"advanced atom properties"));
        editor.mol.atoms[0].query = Some(".s*".into());
        assert!(editor
            .export_losses("query.mol")
            .contains(&"advanced atom properties"));
    }

    #[test]
    fn v3000_query_import_handles_unspecified_and_rejects_invalid_values() {
        let file = |attributes: &str| {
            format!("Q\n  Molecula\n\n  0  0  0     0  0            V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 {attributes}\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n")
        };
        let unspecified = parse_mol_at(&file("SUBST=0 UNSAT=0 RBCNT=0"), vec2(0., 0.)).unwrap();
        assert_eq!(unspecified.atoms[0].query, None);
        for attributes in ["SUBST=-2", "SUBST=256", "UNSAT=2", "RBCNT=-2", "RBCNT=256"] {
            assert!(
                parse_mol_at(&file(attributes), vec2(0., 0.)).is_err(),
                "{attributes}"
            );
        }
    }

    #[test]
    fn v3000_bond_queries_and_ring_topology_roundtrip() {
        for (record, attr, query) in [
            ("5", "", Some("single-or-double")),
            ("6", "", Some("single-or-aromatic")),
            ("7", "", Some("double-or-aromatic")),
            ("8", "", Some("any")),
            ("1", "TOPO=1", Some("topology-ring")),
            ("1", "TOPO=2", Some("topology-chain")),
        ] {
            let source = format!("Q\n  Molecula\n\n  0  0  0     0  0            V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 0\nM  V30 2 C 1 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 {record} 1 2 {attr}\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n");
            let parsed = parse_mol_at(&source, vec2(0., 0.)).unwrap();
            assert_eq!(parsed.bonds[0].query.as_deref(), query);
            let exported = molfile_for(&parsed);
            assert!(exported.contains("V3000"));
            if !attr.is_empty() {
                assert!(exported.contains(attr));
            } else {
                assert!(exported.contains(&format!("M  V30 1 {record} ")));
            }
            let roundtrip = parse_mol_at(&exported, vec2(0., 0.)).unwrap();
            assert_eq!(roundtrip.bonds[0].query.as_deref(), query);
            let mut editor = Editor::new();
            editor.mol = parsed;
            assert!(!editor
                .export_losses("query.mol")
                .contains(&"bond styles or query data"));
        }
        let invalid = "Q\n  Molecula\n\n  0  0  0     0  0            V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 2 1 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0 0 0 0 0\nM  V30 2 C 1 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2 TOPO=3\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
        assert!(parse_mol_at(invalid, vec2(0., 0.)).is_err());
        let conflicting = invalid.replace("TOPO=3", "TOPO=1 TOPO=2");
        assert!(parse_mol_at(&conflicting, vec2(0., 0.))
            .unwrap_err()
            .contains("Conflicting V3000 TOPO"));
    }

    #[test]
    fn v3000_atom_lists_and_reserved_query_atoms_roundtrip() {
        for (type_text, expected_query) in [
            ("[C,N]", "C,N"),
            ("NOT [O,S]", "!O,S"),
            ("A", "any atom"),
            ("Q", "heteroatom"),
            ("*", "any atom"),
        ] {
            let source = format!("Q\n  Molecula\n\n  0  0  0     0  0            V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 {type_text} 0 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n");
            let parsed = parse_mol_at(&source, vec2(0., 0.)).unwrap();
            assert_eq!(parsed.atoms[0].query.as_deref(), Some(expected_query));
            let exported = molfile_for(&parsed);
            assert!(exported.contains("V3000"));
            let roundtrip = parse_mol_at(&exported, vec2(0., 0.)).unwrap();
            assert_eq!(roundtrip.atoms[0].query.as_deref(), Some(expected_query));
            let mut editor = Editor::new();
            editor.mol = parsed;
            assert!(!editor
                .export_losses("query.mol")
                .contains(&"advanced atom properties"));
        }
        for invalid_type in ["[C,NotAnElement]", "[C,C]", "NOT C", "[]"] {
            let source = format!("Q\n  Molecula\n\n  0  0  0     0  0            V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 {invalid_type} 0 0 0 0 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n");
            assert!(
                parse_mol_at(&source, vec2(0., 0.)).is_err(),
                "{invalid_type}"
            );
        }
    }

    #[test]
    fn mol_v2000_roundtrip_keeps_local_chiral_parity_and_stereo_bond_direction() {
        let mut mol = Molecule::default();
        let mut a = Atom::new(vec2(0., 0.), "C");
        a.stereo = Some("@@".into());
        a.atom_map = Some(42);
        mol.atoms.push(a);
        mol.atoms.push(Atom::new(vec2(40., 0.), "F"));
        mol.bonds.push(Bond {
            a: 0,
            b: 1,
            order: 1,
            style: BondStyle::Hash,
            stereo: None,
            query: None,
        });
        let got = parse_mol_at(&molfile_for(&mol), vec2(0., 0.)).unwrap();
        assert_eq!(got.atoms[0].stereo.as_deref(), Some("@@"));
        assert_eq!(got.atoms[0].atom_map, Some(42));
        assert_eq!(got.bonds[0].style, BondStyle::Hash);
    }

    #[test]
    fn native_v3_roundtrip_keeps_separate_3d_coordinates() {
        let mut mol = Molecule::default();
        let mut atom = Atom::new(vec2(2., 3.), "O");
        atom.coord3d = Some(vec3(1.25, -2.5, 3.75));
        mol.atoms.push(atom);
        let snapshot = Snapshot {
            mol,
            graphics: vec![],
        };
        let got = native::read(&native::write(&snapshot)).unwrap();
        assert_eq!(got.mol.atoms[0].coord3d, Some(vec3(1.25, -2.5, 3.75)));
        assert_eq!(got.mol.atoms[0].pos, vec2(2., 3.));
    }

    #[test]
    fn rxn_roundtrip_preserves_component_roles_and_graphs() {
        let mut mol = Molecule::default();
        for (x, role) in [(0., "reactant"), (100., "product"), (50., "agent")] {
            let mut atom = Atom::new(vec2(x, 0.), "O");
            atom.reaction_role = Some(role.into());
            mol.atoms.push(atom);
        }
        let text = rxn_for(&mol).unwrap();
        let got = parse_rxn(&text, vec2(200., 200.)).unwrap();
        assert_eq!(got.atoms.len(), 3);
        assert_eq!(
            got.atoms
                .iter()
                .map(|a| a.reaction_role.as_deref().unwrap())
                .collect::<Vec<_>>(),
            vec!["reactant", "product", "agent"]
        );
    }

    #[test]
    fn functional_group_and_amino_acid_templates_are_connected_graphs() {
        for (name, atoms, bonds) in [
            ("Ester", 4, 3),
            ("Amide", 4, 3),
            ("Glycine", 5, 4),
            ("Alanine", 6, 5),
            ("Nitrile", 2, 1),
        ] {
            let mut e = Editor::new();
            e.insert_template(name, vec2(0., 0.));
            assert_eq!(e.mol.atoms.len(), atoms, "{name}");
            assert_eq!(e.mol.bonds.len(), bonds, "{name}");
            assert!(
                e.mol.bonds.iter().all(|b| b.a < atoms && b.b < atoms),
                "{name}"
            );
        }
    }

    #[test]
    fn fused_and_bridged_ring_templates_are_real_undoable_chemical_graphs() {
        for (name, formula, atoms, bonds, rings) in [
            ("Naphthalene", "C10H8", 10, 11, 2),
            ("Decalin", "C10H18", 10, 11, 2),
            ("Indole", "C8H7N", 9, 10, 2),
            ("Quinoline", "C9H7N", 10, 11, 2),
            ("Anthracene", "C14H10", 14, 16, 3),
            ("Phenanthrene", "C14H10", 14, 16, 3),
            ("Norbornane", "C7H12", 7, 8, 2),
        ] {
            assert!(builtin_template_matches_query(name, ""));
            assert!(builtin_template_matches_query(name, "ring"));
            let mut editor = Editor::new();
            editor.insert_template(name, vec2(100., 100.));
            assert_eq!(editor.mol.atoms.len(), atoms, "{name} atom count");
            assert_eq!(editor.mol.bonds.len(), bonds, "{name} bond count");
            let analysis = chemistry::analyze(&editor.mol);
            assert_eq!(analysis.formula, formula, "{name}");
            assert_eq!(analysis.connected_components, 1, "{name}");
            assert_eq!(analysis.cycle_rank, rings, "{name}");
            assert!(
                analysis.warnings.is_empty(),
                "{name}: {:?}",
                analysis.warnings
            );
            if matches!(
                name,
                "Naphthalene" | "Indole" | "Quinoline" | "Anthracene" | "Phenanthrene"
            ) {
                assert!(editor.mol.atoms.iter().all(|atom| atom.aromatic));
                assert!(editor.mol.bonds.iter().all(|bond| bond.order == 4));
            }
            editor.undo();
            assert!(
                editor.mol.atoms.is_empty(),
                "{name} insertion is one undoable action"
            );
        }
    }

    #[test]
    fn proteinogenic_amino_acid_templates_have_expected_formulas_and_graphs() {
        for (name, formula) in [
            ("Glycine", "C2H5NO2"),
            ("Alanine", "C3H7NO2"),
            ("Valine", "C5H11NO2"),
            ("Leucine", "C6H13NO2"),
            ("Isoleucine", "C6H13NO2"),
            ("Serine", "C3H7NO3"),
            ("Threonine", "C4H9NO3"),
            ("Cysteine", "C3H7NO2S"),
            ("Methionine", "C5H11NO2S"),
            ("Aspartic acid", "C4H7NO4"),
            ("Asparagine", "C4H8N2O3"),
            ("Glutamic acid", "C5H9NO4"),
            ("Glutamine", "C5H10N2O3"),
            ("Lysine", "C6H14N2O2"),
            ("Arginine", "C6H14N4O2"),
            ("Histidine", "C6H9N3O2"),
            ("Phenylalanine", "C9H11NO2"),
            ("Tyrosine", "C9H11NO3"),
            ("Tryptophan", "C11H12N2O2"),
            ("Proline", "C5H9NO2"),
        ] {
            assert!(
                builtin_template_matches_query(name, "amino acid"),
                "{name} search alias"
            );
            let mut editor = Editor::new();
            editor.insert_template(name, vec2(0., 0.));
            let analysis = chemistry::analyze(&editor.mol);
            assert_eq!(analysis.formula, formula, "{name}");
            assert_eq!(analysis.connected_components, 1, "{name}");
            assert!(!editor.mol.atoms.is_empty(), "{name}");
            editor.undo();
            assert!(editor.mol.atoms.is_empty(), "{name}");
        }
        assert!(!builtin_template_matches_query("Benzene", "proteinogenic"));
    }

    #[test]
    fn built_in_and_custom_templates_rotate_before_placement() {
        let center = vec2(250., 180.);
        let mut normal = Editor::new();
        normal.insert_template("Ethanol", center);
        let mut rotated = Editor::new();
        rotated.template_rotation = 60;
        rotated.insert_template("Ethanol", center);
        assert_eq!(rotated.template_rotation, 0);
        assert_eq!(rotated.mol.bonds, normal.mol.bonds);
        for (got, original) in rotated.mol.atoms.iter().zip(&normal.mol.atoms) {
            assert!(
                got.pos
                    .distance(rotate_template_point(original.pos, center, 60))
                    < 0.01
            );
        }
        for bond in &rotated.mol.bonds {
            assert!(
                (rotated.mol.atoms[bond.a]
                    .pos
                    .distance(rotated.mol.atoms[bond.b].pos)
                    - normal.mol.atoms[bond.a]
                        .pos
                        .distance(normal.mol.atoms[bond.b].pos))
                .abs()
                    < 0.01
            );
        }

        let mut custom_graph = Molecule::default();
        custom_graph.atoms.push(Atom::new(vec2(-20., 0.), "C"));
        custom_graph.atoms.push(Atom::new(vec2(20., 0.), "O"));
        custom_graph.bonds.push(Bond {
            a: 0,
            b: 1,
            order: 1,
            style: BondStyle::Normal,
            stereo: None,
            query: None,
        });
        let mut custom = Editor::new();
        custom.template_rotation = 90;
        custom.pending_custom_template = Some(custom_graph);
        custom.insert_template("User", center);
        assert!((custom.mol.atoms[0].pos.distance(custom.mol.atoms[1].pos) - 40.0).abs() < 0.01);
        assert!((custom.mol.atoms[0].pos.x - center.x).abs() < 0.01);
        assert!((custom.mol.atoms[1].pos.x - center.x).abs() < 0.01);
        assert_eq!(custom.template_rotation, 0);
    }

    #[test]
    fn template_activation_prepares_ghost_graphs_for_builtin_and_custom_templates() {
        let mut editor = Editor::new();
        editor.activate_builtin_template("Benzene");
        let preview = editor.template_preview.as_ref().unwrap();
        assert_eq!(preview.atoms.len(), 6);
        assert_eq!(preview.bonds.len(), 6);
        assert!(matches!(editor.tool, Tool::Template("Benzene")));
        assert!(editor.pending_custom_template.is_none());

        let mut custom = Molecule {
            name: "Offset pair".into(),
            ..Molecule::default()
        };
        custom.atoms.push(Atom::new(vec2(20., -10.), "C"));
        custom.atoms.push(Atom::new(vec2(60., -10.), "O"));
        custom.bonds.push(Bond {
            a: 0,
            b: 1,
            order: 1,
            style: BondStyle::Normal,
            stereo: None,
            query: None,
        });
        editor.activate_custom_template(custom);
        let preview = editor.template_preview.as_ref().unwrap();
        assert_eq!(preview.atoms.len(), 2);
        assert_eq!(preview.bonds.len(), 1);
        assert!((preview.atoms[0].pos.x + 20.).abs() < 0.001);
        assert!((preview.atoms[1].pos.x - 20.).abs() < 0.001);
        assert!(editor.pending_custom_template.is_some());
        assert!(matches!(editor.tool, Tool::Template("User")));
        editor.template_rotation = 120;
        editor.cancel_template_placement();
        assert_eq!(editor.template_preview, None);
        assert_eq!(editor.pending_custom_template, None);
        assert_eq!(editor.template_rotation, 0);
        assert!(matches!(editor.tool, Tool::Select));
    }

    #[test]
    fn attached_template_preview_aligns_anchor_and_insertion_is_one_undoable_graph_edit() {
        let mut editor = Editor::new();
        editor.mol.atoms.push(Atom::new(vec2(200., 180.), "C"));
        editor.selected = Some(0);
        editor.selected_atoms = vec![0];
        editor.activate_builtin_template("Benzene");
        editor.toggle_template_attachment();
        assert_eq!(editor.template_attachment_target, Some(0));
        editor.template_rotation = 60;
        let cursor = vec2(200., 320.);
        let center = editor.template_placement_center(cursor);
        let (anchor_index, anchor) = editor.template_anchor().unwrap();
        let anchor_world = center + rotate_template_point(anchor, vec2(0., 0.), 60);
        assert!((anchor_world.distance(editor.mol.atoms[0].pos) - BOND).abs() < 0.01);

        editor.insert_template("Benzene", cursor);
        assert_eq!(editor.mol.atoms.len(), 7);
        assert_eq!(editor.mol.bonds.len(), 7);
        assert_eq!(molecule_components(&editor.mol).len(), 1);
        let attachment = editor
            .mol
            .bonds
            .iter()
            .find(|bond| bond.a == 0 || bond.b == 0)
            .unwrap();
        let inserted_anchor = if attachment.a == 0 {
            attachment.b
        } else {
            attachment.a
        };
        assert_eq!(inserted_anchor, anchor_index + 1);
        assert!(
            (editor.mol.atoms[inserted_anchor]
                .pos
                .distance(editor.mol.atoms[0].pos)
                - BOND)
                .abs()
                < 0.01
        );

        editor.undo();
        assert_eq!(editor.mol.atoms.len(), 1);
        assert!(editor.mol.bonds.is_empty());

        let mut custom = Editor::new();
        custom.mol.atoms.push(Atom::new(vec2(80., 40.), "C"));
        custom.selected = Some(0);
        custom.selected_atoms = vec![0];
        let mut pair = Molecule {
            name: "C–O pair".into(),
            ..Molecule::default()
        };
        pair.atoms.push(Atom::new(vec2(-20., 0.), "C"));
        pair.atoms.push(Atom::new(vec2(20., 0.), "O"));
        pair.bonds.push(Bond {
            a: 0,
            b: 1,
            order: 1,
            style: BondStyle::Normal,
            stereo: None,
            query: None,
        });
        custom.activate_custom_template(pair);
        custom.toggle_template_attachment();
        custom.template_rotation = 90;
        custom.insert_template("User", vec2(80., 160.));
        assert_eq!(custom.mol.atoms.len(), 3);
        assert_eq!(custom.mol.bonds.len(), 2);
        assert_eq!(molecule_components(&custom.mol).len(), 1);
        custom.undo();
        assert_eq!(custom.mol.atoms.len(), 1);
        assert!(custom.mol.bonds.is_empty());
    }

    #[test]
    fn custom_template_storage_roundtrips_graph_and_rejects_duplicate_names() {
        let directory = std::env::temp_dir().join(format!(
            "molecula-user-template-storage-{}",
            std::process::id()
        ));
        let _ = fs::remove_dir_all(&directory);
        let mut molecule = Molecule::default();
        let mut carbon = Atom::new(vec2(-20., 0.), "C");
        carbon.isotope = Some(13);
        carbon.alias = Some("Me".into());
        molecule.atoms.push(carbon);
        molecule.atoms.push(Atom::new(vec2(20., 0.), "O"));
        molecule.bonds.push(Bond {
            a: 0,
            b: 1,
            order: 1,
            style: BondStyle::Normal,
            stereo: None,
            query: None,
        });
        let saved = save_user_template_to(&directory, "Methyl ether", &molecule).unwrap();
        assert_eq!(saved.name, "Methyl ether");
        let loaded = load_user_templates_from(&directory).unwrap();
        assert_eq!(loaded.len(), 1);
        assert_eq!(loaded[0].name, "Methyl ether");
        assert_eq!(loaded[0].mol.atoms.len(), 2);
        assert_eq!(loaded[0].mol.bonds.len(), 1);
        assert_eq!(loaded[0].mol.atoms[0].isotope, Some(13));
        assert_eq!(loaded[0].mol.atoms[0].alias.as_deref(), Some("Me"));
        assert!(save_user_template_to(&directory, "Methyl ether", &molecule).is_err());
        let _existing = save_user_template_to(&directory, "Existing", &molecule).unwrap();
        assert!(save_user_template_to(&directory, "../unsafe", &molecule).is_err());
        let mut invalid = molecule.clone();
        invalid.bonds[0].b = 99;
        assert!(save_user_template_to(&directory, "Invalid", &invalid).is_err());
        let mut editor = Editor::new();
        editor.user_templates = load_user_templates_from(&directory).unwrap();
        let original_index = editor
            .user_templates
            .iter()
            .position(|template| template.name == "Methyl ether")
            .unwrap();
        assert!(editor
            .rename_user_template_in(original_index, "Existing", &directory)
            .is_err());
        assert!(
            saved.path.exists(),
            "duplicate rename must preserve original"
        );
        assert!(editor
            .rename_user_template_in(original_index, "Ether group", &directory)
            .is_ok());
        assert!(!saved.path.exists());
        let renamed = load_user_templates_from(&directory).unwrap();
        assert_eq!(renamed.len(), 2);
        let ether = renamed
            .iter()
            .find(|template| template.name == "Ether group")
            .unwrap();
        assert_eq!(ether.mol.atoms[0].isotope, Some(13));
        assert_eq!(ether.mol.atoms[0].alias.as_deref(), Some("Me"));
        assert_eq!(ether.mol.bonds.len(), 1);
        let renamed_index = editor
            .user_templates
            .iter()
            .position(|template| template.name == "Ether group")
            .unwrap();
        editor.delete_user_template(renamed_index).unwrap();
        let existing_index = editor
            .user_templates
            .iter()
            .position(|template| template.name == "Existing")
            .unwrap();
        editor.delete_user_template(existing_index).unwrap();
        assert!(!saved.path.exists());
        fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn selected_custom_template_extracts_connected_graph_and_inserts_undoably() {
        let mut editor = Editor::new();
        let c = editor.add_atom(vec2(0., 0.), "C");
        let o = editor.add_atom(vec2(43., 0.), "O");
        let n = editor.add_atom(vec2(220., 0.), "N");
        editor.add_bond(c, o, 2);
        editor.selected_atoms = vec![c, o];
        editor.selected = Some(c);
        editor.mol.atoms[c].reaction_role = Some("reactant".into());
        editor.mol.atoms[c].atom_map = Some(7);

        let fragment = editor.selected_fragment_for_template().unwrap();
        assert_eq!(fragment.atoms.len(), 2);
        assert_eq!(fragment.bonds.len(), 1);
        assert!(fragment
            .atoms
            .iter()
            .all(|atom| atom.reaction_role.is_none()));
        assert!(fragment.atoms.iter().all(|atom| atom.atom_map.is_none()));
        assert!((fragment.atoms[0].pos + fragment.atoms[1].pos).length() < 0.001);
        editor.selected_atoms = vec![c, n];
        assert!(editor.selected_fragment_for_template().is_err());

        editor.insert_custom_template(fragment, vec2(400., 300.));
        assert_eq!(editor.mol.atoms.len(), 5);
        assert_eq!(editor.mol.bonds.len(), 2);
        assert!(editor
            .mol
            .bonds
            .iter()
            .any(|bond| bond.a == 3 && bond.b == 4));
        assert_eq!(editor.mol.atoms[n].element, "N");
        editor.undo();
        assert_eq!(editor.mol.atoms.len(), 3);
        assert_eq!(editor.mol.bonds.len(), 1);
    }

    #[test]
    fn pubchem_query_encoding_is_safe_and_stable() {
        assert_eq!(url_encode("salicylic acid"), "salicylic%20acid");
        assert_eq!(url_encode("C/C=C\\O"), "C%2FC%3DC%5CO");
    }

    #[test]
    fn png_and_pdf_exports_are_real_documents_with_structure_content() {
        let mut mol = Molecule::default();
        mol.atoms.push(Atom::new(vec2(100., 100.), "C"));
        mol.atoms.push(Atom::new(vec2(170., 100.), "O"));
        mol.bonds.push(Bond {
            a: 0,
            b: 1,
            order: 2,
            style: BondStyle::Normal,
            stereo: None,
            query: None,
        });
        let png = render_png(&mol, &[], 640, 360).unwrap();
        assert_eq!(&png[..8], b"\x89PNG\r\n\x1a\n");
        assert!(image::load_from_memory(&png).is_ok());
        let pdf = render_pdf(&mol, &[], 842., 595.);
        let text = String::from_utf8(pdf).unwrap();
        assert!(text.starts_with("%PDF-1.4"));
        assert!(text.contains("Tj ET"));
        assert!(text.ends_with("%%EOF\n"));
    }

    #[test]
    fn svg_png_and_pdf_exports_preserve_bond_display_styles() {
        let styles = [
            (BondStyle::Normal, 2, 0),
            (BondStyle::Aromatic, 9, 0),
            (BondStyle::Wedge, 0, 1),
            (BondStyle::Hash, 11, 0),
            (BondStyle::Wavy, 24, 0),
            (BondStyle::Coordinate, 3, 0),
            (BondStyle::Any, 8, 0),
        ];
        for (style, expected_lines, expected_triangles) in styles {
            let mut editor = Editor::new();
            editor.mol.atoms = vec![
                Atom::new(vec2(100., 100.), "C"),
                Atom::new(vec2(170., 100.), "C"),
            ];
            editor.mol.bonds.push(Bond {
                a: 0,
                b: 1,
                order: if style == BondStyle::Normal { 2 } else { 1 },
                style,
                stereo: None,
                query: None,
            });
            let svg = editor.svg_with_size(640., 360.);
            assert_eq!(svg.matches("<line ").count(), expected_lines, "{style:?}");
            assert_eq!(
                svg.matches("<polygon ").count(),
                expected_triangles,
                "{style:?}"
            );

            let pdf = String::from_utf8(render_pdf(&editor.mol, &[], 842., 595.)).unwrap();
            assert_eq!(
                pdf.matches(" m ").count(),
                expected_lines + expected_triangles,
                "{style:?}"
            );
            assert_eq!(
                pdf.matches(" h f\n").count(),
                expected_triangles,
                "{style:?}"
            );

            let png = render_png(&editor.mol, &[], 640, 360).unwrap();
            let image = image::load_from_memory(&png).unwrap().to_rgba8();
            assert!(image.pixels().any(|pixel| pixel[0] < 100), "{style:?}");
        }
    }

    #[test]
    fn drawing_export_keeps_the_editable_document_dirty_and_untitled() {
        let path = std::env::temp_dir().join(format!("molecula-export-{}.png", std::process::id()));
        let mut e = Editor::new();
        e.mol.atoms.push(Atom::new(vec2(80., 80.), "N"));
        e.dirty = true;
        e.save(path.to_str().unwrap());
        assert!(e.pending_export.is_some());
        e.save_confirmed(path.to_str().unwrap());
        assert!(e.dirty);
        assert!(e.documents[0].path.is_none());
        assert!(path.exists());
        let _ = fs::remove_file(path);
    }

    #[test]
    fn naming_lesson_reveals_real_steps_progressively() {
        let mut e = Editor::new();
        let a = e.add_atom(vec2(0., 0.), "C");
        let b = e.add_atom(vec2(BOND, 0.), "C");
        let o = e.add_atom(vec2(BOND * 2., 0.), "O");
        e.add_bond(a, b, 1);
        e.add_bond(b, o, 1);
        e.start_naming_lesson();
        assert!(e.lesson_open);
        let lesson = e.lesson.as_ref().unwrap();
        assert_eq!(lesson.step, 0);
        assert!(
            matches!(&lesson.content,LessonContent::Naming{name,steps} if !name.is_empty()&&steps.len()==6)
        );
    }

    #[test]
    fn reaction_lesson_grades_balanced_graphs_and_never_grades_missing_data() {
        let mut e = Editor::new();
        for (x, role) in [(0., "reactant"), (100., "product")] {
            let a = e.add_atom(vec2(x, 0.), "C");
            let b = e.add_atom(vec2(x + BOND, 0.), "O");
            e.add_bond(a, b, 1);
            for i in [a, b] {
                e.mol.atoms[i].reaction_role = Some(role.into());
            }
        }
        e.graphics.push(Graphic::Arrow {
            start: vec2(40., 100.),
            end: vec2(140., 100.),
            kind: "reaction".into(),
            above: String::new(),
            below: String::new(),
        });
        e.start_reaction_lesson();
        e.lesson.as_mut().unwrap().step = 3;
        e.check_reaction_lesson();
        assert!(
            matches!(&e.lesson.as_ref().unwrap().content,LessonContent::Reaction{feedback:Some(s)} if s.starts_with("Correct:"))
        );
        e.mol.atoms[0].query = Some("*".into());
        e.check_reaction_lesson();
        assert!(
            matches!(&e.lesson.as_ref().unwrap().content,LessonContent::Reaction{feedback:Some(s)} if s.starts_with("Not evaluated:"))
        );
    }

    #[test]
    fn language_setting_roundtrips_and_localizes_menu_labels() {
        let path = std::env::temp_dir().join(format!(
            "molecula-language-test-{}.conf",
            std::process::id()
        ));
        save_language_to(&path, Language::Danish).unwrap();
        save_setting(&path, "theme", "dark").unwrap();
        save_implicit_hydrogen_visibility_to(&path, true).unwrap();
        save_zoom_preference_to(&path, 1.75).unwrap();
        assert_eq!(load_language_from(&path), Language::Danish);
        assert!(load_implicit_hydrogen_visibility_from(&path));
        assert_eq!(load_zoom_preference_from(&path), 1.75);
        save_language_to(&path, Language::English).unwrap();
        let settings = fs::read_to_string(&path).unwrap();
        assert!(settings.contains("theme=dark"));
        assert!(settings.contains("show_implicit_hydrogens=true"));
        assert!(settings.contains("zoom=1.7500"));
        assert_eq!(localized(Language::Danish, "Calculations"), "Beregninger");
        assert_eq!(localized(Language::English, "Calculations"), "Calculations");
        assert_eq!(
            localized(Language::Danish, "Reaction coefficient: custom"),
            "Reaktionskoefficient: brugerdefineret"
        );
        assert_eq!(
            localized(Language::Danish, "Gasteiger partial charges"),
            "Gasteiger-partialladninger"
        );
        let _ = fs::remove_file(path);
    }

    #[test]
    fn quick_start_guide_localizes_and_explains_in_place_bond_editing() {
        let english = quick_start_guide(Language::English);
        assert!(
            english.contains("To change an existing bond, choose Double (=) and click its middle.")
        );
        assert!(english.contains("UNDO Cmd/Ctrl+Z"));

        let danish = quick_start_guide(Language::Danish);
        assert!(danish.contains("MOLECULA · HURTIG VEJLEDNING"));
        assert!(danish.contains("vælg Dobbeltbinding (=), og klik midt på bindingen"));
        assert!(danish.contains("Når et atom er valgt, skriv @ efterfulgt af et massetal"));
        assert!(!danish.contains("With an atom selected"));
        assert!(danish.contains("FORTRYD Cmd/Ctrl+Z"));
    }

    #[test]
    fn zoom_preference_clamps_values_and_rejects_non_finite_writes() {
        let path =
            std::env::temp_dir().join(format!("molecula-zoom-test-{}.conf", std::process::id()));
        save_zoom_preference_to(&path, 9.0).unwrap();
        assert_eq!(load_zoom_preference_from(&path), 4.0);
        assert!(save_zoom_preference_to(&path, f32::NAN).is_err());
        fs::write(&path, "zoom=not-a-number\n").unwrap();
        assert_eq!(load_zoom_preference_from(&path), 1.0);
        fs::write(&path, "zoom=0.01\n").unwrap();
        assert_eq!(load_zoom_preference_from(&path), 0.2);
        let _ = fs::remove_file(path);
    }

    #[test]
    fn zoom_changes_become_the_default_for_new_documents() {
        let mut editor = Editor::new();
        editor.zoom = 1.0;
        editor.preferred_zoom = 1.0;
        editor.adjust_zoom(1.5, vec2(640., 360.), 1280., 720.);
        assert_eq!(editor.preferred_zoom, 1.5);
        assert!(editor.zoom_preference_dirty);

        editor.new_document();
        assert_eq!(editor.zoom, 1.5);
    }

    #[test]
    fn changed_zoom_is_written_once_without_clobbering_other_preferences() {
        let path =
            std::env::temp_dir().join(format!("molecula-zoom-persist-{}.conf", std::process::id()));
        save_setting(&path, "theme", "dark").unwrap();
        let mut editor = Editor::new();
        editor.preferred_zoom = 2.25;
        editor.zoom_preference_dirty = true;

        editor.persist_zoom_preference_to(&path);

        assert!(!editor.zoom_preference_dirty);
        assert_eq!(load_zoom_preference_from(&path), 2.25);
        assert!(fs::read_to_string(&path).unwrap().contains("theme=dark"));
        let _ = fs::remove_file(path);
    }

    #[test]
    fn implicit_hydrogen_labels_are_optional_and_respect_atom_aliases() {
        let carbon = Atom::new(vec2(0., 0.), "C");
        assert_eq!(atom_display_label(&carbon, Some(3), false), "C");
        assert_eq!(atom_display_label(&carbon, Some(3), true), "CH3");
        let mut oxygen = Atom::new(vec2(0., 0.), "O");
        oxygen.charge = -1;
        assert_eq!(atom_display_label(&oxygen, Some(0), true), "O-");
        let mut alias = Atom::new(vec2(0., 0.), "C");
        alias.alias = Some("Me".into());
        assert_eq!(atom_display_label(&alias, Some(3), true), "Me");
    }

    #[test]
    fn dirty_tab_close_requires_an_explicit_decision_and_keeps_other_tabs() {
        let mut e = Editor::new();
        e.new_document();
        e.switch_document(0);
        e.checkpoint();
        e.add_atom(vec2(10., 10.), "C");
        e.dirty = true;
        e.sync_active_document();
        e.request_close_document(0);
        assert_eq!(e.pending_close, Some(0));
        assert_eq!(e.documents.len(), 2);
        e.pending_close = None;
        e.request_close_document(0);
        e.close_document(0);
        assert_eq!(e.documents.len(), 1);
        assert_eq!(e.active_document, 0);
        assert!(e.mol.atoms.is_empty());
    }

    #[test]
    fn saving_a_dirty_tab_completes_close_only_after_write_succeeds() {
        let path =
            std::env::temp_dir().join(format!("molecula-close-save-{}.mols", std::process::id()));
        let mut e = Editor::new();
        e.add_atom(vec2(10., 10.), "O");
        e.dirty = true;
        e.sync_active_document();
        e.documents[0].path = Some(path.to_string_lossy().into_owned());
        e.request_close_document(0);
        e.close_pending_with_save();
        assert!(path.exists());
        assert_eq!(e.documents.len(), 1);
        assert!(e.mol.atoms.is_empty());
        let _ = fs::remove_file(path);
    }

    #[test]
    fn autosave_recovery_waits_for_user_choice() {
        let mut e = Editor::new();
        let mut snapshot = Snapshot {
            mol: Molecule::default(),
            graphics: vec![],
        };
        snapshot.mol.atoms.push(Atom::new(vec2(4., 5.), "N"));
        let docs = vec![DocumentTab {
            title: "Recovered".into(),
            path: None,
            snapshot,
            dirty: true,
            undo: vec![],
            redo: vec![],
        }];
        e.stage_recovery(&native::write_workspace(&docs, 0))
            .unwrap();
        assert!(e.mol.atoms.is_empty());
        assert!(e.pending_recovery.is_some());
        e.apply_recovery(true);
        assert_eq!(e.mol.atoms.len(), 1);
        assert!(e.dirty);
    }

    #[test]
    fn two_click_arrows_are_real_undoable_objects() {
        let mut e = Editor::new();
        e.tool = Tool::Arrow;
        e.drag_start = Some(vec2(100., 100.));
        e.finish_arrow(vec2(100., 100.));
        assert!(e.graphics.is_empty() && e.pending_arrow.is_some());
        e.drag_start = Some(vec2(300., 100.));
        e.finish_arrow(vec2(300., 100.));
        assert!(matches!(e.graphics.first(), Some(Graphic::Arrow { .. })));
        e.undo();
        assert!(e.graphics.is_empty());
        e.tool = Tool::Electron(2);
        e.drag_start = Some(vec2(100., 100.));
        e.finish_arrow(vec2(100., 100.));
        e.drag_start = Some(vec2(200., 100.));
        e.finish_arrow(vec2(200., 100.));
        assert!(matches!(
            e.graphics.first(),
            Some(Graphic::ElectronArrow { electrons: 2, .. })
        ));
    }

    #[test]
    fn removing_atom_reindexes_remaining_bonds() {
        let mut editor = Editor::new();
        for x in [0., 40., 80., 120.] {
            editor.add_atom(vec2(x, 0.), "C");
        }
        editor.add_bond(0, 2, 1);
        editor.add_bond(2, 3, 2);
        editor.remove_atom(1);
        assert_eq!(editor.mol.atoms.len(), 3);
        assert_eq!((editor.mol.bonds[0].a, editor.mol.bonds[0].b), (0, 1));
        assert_eq!((editor.mol.bonds[1].a, editor.mol.bonds[1].b), (1, 2));
    }

    #[test]
    fn editing_existing_bond_does_not_duplicate_it() {
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "C");
        editor.add_atom(vec2(40., 0.), "O");
        editor.add_bond(0, 1, 1);
        editor.add_bond(1, 0, 2);
        assert_eq!(editor.mol.bonds.len(), 1);
        assert_eq!(editor.mol.bonds[0].order, 2);
    }

    #[test]
    fn bond_inspector_edits_order_style_and_ez_with_undo_redo() {
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "C");
        editor.add_atom(vec2(40., 0.), "C");
        editor.add_bond(0, 1, 1);
        editor.selected_bond = Some(0);
        editor.selected_bonds = vec![0];

        editor.set_bond_order(0, 2);
        assert_eq!(editor.mol.bonds[0].order, 2);
        editor.cycle_double_bond_stereo(0);
        assert_eq!(editor.mol.bonds[0].stereo.as_deref(), Some("E"));
        editor.cycle_double_bond_stereo(0);
        assert_eq!(editor.mol.bonds[0].stereo.as_deref(), Some("Z"));
        editor.set_bond_order(0, 4);
        assert_eq!(editor.mol.bonds[0].style, BondStyle::Aromatic);
        assert_eq!(editor.mol.bonds[0].stereo, None);
        assert!(editor.mol.atoms.iter().all(|atom| atom.aromatic));
        editor.cycle_bond_style(0);
        assert_eq!(editor.mol.bonds[0].order, 1);
        assert_eq!(editor.mol.bonds[0].style, BondStyle::Wedge);
        assert!(editor.mol.atoms.iter().all(|atom| !atom.aromatic));

        editor.undo();
        assert_eq!(editor.mol.bonds[0].order, 4);
        assert_eq!(editor.mol.bonds[0].style, BondStyle::Aromatic);
        editor.redo();
        assert_eq!(editor.mol.bonds[0].style, BondStyle::Wedge);
        assert_eq!(editor.selected_bond, None);
    }

    #[test]
    fn clicking_existing_bond_with_double_tool_changes_order_in_one_undoable_edit() {
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "C");
        editor.add_atom(vec2(40., 0.), "C");
        editor.add_bond(0, 1, 1);
        editor.tool = Tool::Bond(2);

        assert!(editor.apply_bond_tool_to_hit_bond(vec2(20., 0.)));
        assert_eq!(editor.mol.bonds.len(), 1);
        assert_eq!(editor.mol.bonds[0].order, 2);
        assert_eq!(editor.selected_bond, Some(0));

        editor.undo();
        assert_eq!(editor.mol.bonds[0].order, 1);
        editor.redo();
        assert_eq!(editor.mol.bonds[0].order, 2);
        assert_eq!(editor.mol.bonds.len(), 1);
    }

    #[test]
    fn setting_a_concrete_bond_order_clears_only_order_alternatives_and_is_undoable() {
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "C");
        editor.add_atom(vec2(40., 0.), "C");
        editor.mol.bonds.push(Bond {
            a: 0,
            b: 1,
            order: 1,
            style: BondStyle::Any,
            stereo: Some("E".into()),
            query: Some("single-or-double;topology-ring".into()),
        });
        editor.tool = Tool::Bond(2);

        assert!(editor.apply_bond_tool_to_hit_bond(vec2(20., 0.)));
        let bond = &editor.mol.bonds[0];
        assert_eq!(bond.order, 2);
        assert_eq!(bond.style, BondStyle::Normal);
        assert_eq!(bond.stereo, None);
        assert_eq!(bond.query.as_deref(), Some("topology-ring"));

        editor.undo();
        assert_eq!(editor.mol.bonds[0].order, 1);
        assert_eq!(editor.mol.bonds[0].style, BondStyle::Any);
        assert_eq!(
            editor.mol.bonds[0].query.as_deref(),
            Some("single-or-double;topology-ring")
        );
        editor.redo();
        assert_eq!(editor.mol.bonds[0].order, 2);
        assert_eq!(editor.mol.bonds[0].query.as_deref(), Some("topology-ring"));
    }

    #[test]
    fn query_bond_tool_updates_an_existing_bond_instead_of_adding_a_parallel_edge() {
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "C");
        editor.add_atom(vec2(40., 0.), "C");
        editor.add_bond(0, 1, 1);
        editor.tool = Tool::StyledBond(BondStyle::Any);
        editor.bond_query = Some("double-or-aromatic".into());

        assert!(editor.apply_bond_tool_to_hit_bond(vec2(20., 0.)));
        assert_eq!(editor.mol.bonds.len(), 1);
        assert_eq!(editor.mol.bonds[0].style, BondStyle::Any);
        assert_eq!(
            editor.mol.bonds[0].query.as_deref(),
            Some("double-or-aromatic")
        );
        editor.undo();
        assert_eq!(editor.mol.bonds[0].style, BondStyle::Normal);
        assert_eq!(editor.mol.bonds[0].query, None);
    }

    #[test]
    fn bond_tool_click_on_atom_keeps_draw_from_atom_behavior() {
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "C");
        editor.add_atom(vec2(40., 0.), "C");
        editor.add_bond(0, 1, 1);
        editor.tool = Tool::Bond(2);

        assert!(!editor.apply_bond_tool_to_hit_bond(vec2(0., 0.)));
        assert_eq!(editor.mol.bonds[0].order, 1);
        assert_eq!(editor.selected_bond, None);
    }

    #[test]
    fn changing_double_bond_to_wedge_clears_ez_and_is_undoable() {
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "C");
        editor.add_atom(vec2(40., 0.), "C");
        editor.add_bond(0, 1, 2);
        editor.mol.bonds[0].stereo = Some("E".into());

        editor.cycle_bond_style(0);
        assert_eq!(editor.mol.bonds[0].order, 1);
        assert_eq!(editor.mol.bonds[0].style, BondStyle::Wedge);
        assert_eq!(editor.mol.bonds[0].stereo, None);

        editor.undo();
        assert_eq!(editor.mol.bonds[0].order, 2);
        assert_eq!(editor.mol.bonds[0].stereo.as_deref(), Some("E"));
        editor.redo();
        assert_eq!(editor.mol.bonds[0].style, BondStyle::Wedge);
        assert_eq!(editor.mol.bonds[0].stereo, None);
    }

    #[test]
    fn bond_inspector_cycles_query_topology_and_undoes_it() {
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "C");
        editor.add_atom(vec2(40., 0.), "C");
        editor.add_bond(0, 1, 1);
        editor.mol.bonds[0].query = Some("single-or-double".into());

        editor.cycle_bond_topology(0);
        assert_eq!(
            editor.mol.bonds[0].query.as_deref(),
            Some("single-or-double;topology-ring")
        );
        editor.cycle_bond_topology(0);
        assert_eq!(
            editor.mol.bonds[0].query.as_deref(),
            Some("single-or-double;topology-chain")
        );
        editor.undo();
        assert_eq!(
            editor.mol.bonds[0].query.as_deref(),
            Some("single-or-double;topology-ring")
        );
        editor.cycle_bond_topology(0);
        editor.cycle_bond_topology(0);
        assert_eq!(
            editor.mol.bonds[0].query.as_deref(),
            Some("single-or-double")
        );
    }

    #[test]
    fn selected_atom_accumulates_query_predicates_and_conflicts_are_non_destructive() {
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "C");
        editor.selected = Some(0);
        editor.apply_quick(shortcuts::QuickAction::Query(".a".into()));
        editor.apply_quick(shortcuts::QuickAction::Query(".H0".into()));
        assert_eq!(editor.mol.atoms[0].query.as_deref(), Some(".a;.H0"));
        let undo_len = editor.undo.len();
        editor.apply_quick(shortcuts::QuickAction::Query(".A".into()));
        assert_eq!(editor.mol.atoms[0].query.as_deref(), Some(".a;.H0"));
        assert_eq!(editor.undo.len(), undo_len);
        assert!(editor.status.contains("Conflicting query predicates"));
        editor.undo();
        assert_eq!(editor.mol.atoms[0].query.as_deref(), Some(".a"));
    }

    #[test]
    fn atom_query_inspector_actions_add_and_clear_predicates_with_undo() {
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "C");
        editor.selected = Some(0);
        editor.apply_quick(shortcuts::QuickAction::Query(".a".into()));
        editor.apply_quick(shortcuts::QuickAction::Query(".H1".into()));
        assert_eq!(editor.mol.atoms[0].query.as_deref(), Some(".a;.H1"));
        editor.clear_atom_query(0);
        assert_eq!(editor.mol.atoms[0].query, None);
        editor.undo();
        assert_eq!(editor.mol.atoms[0].query.as_deref(), Some(".a;.H1"));
    }

    #[test]
    fn atom_inspector_cycles_supported_radical_states_with_undo_redo() {
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "C");

        for expected in [1, 2, 3, 0] {
            editor.cycle_atom_radical(0);
            assert_eq!(editor.mol.atoms[0].radical, expected);
        }
        editor.undo();
        assert_eq!(editor.mol.atoms[0].radical, 3);
        editor.redo();
        assert_eq!(editor.mol.atoms[0].radical, 0);
    }

    #[test]
    fn atom_inspector_cycles_implicit_h_override_and_undoes_it() {
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "C");
        for expected in [Some(0), Some(1), Some(2), Some(3), Some(4), Some(5), None] {
            editor.cycle_implicit_h_override(0);
            assert_eq!(editor.mol.atoms[0].implicit_h_override, expected);
        }
        editor.undo();
        assert_eq!(editor.mol.atoms[0].implicit_h_override, Some(5));
        editor.redo();
        assert_eq!(editor.mol.atoms[0].implicit_h_override, None);

        let mut v3000 = Editor::new();
        v3000.add_atom(vec2(0., 0.), "C");
        v3000.mol.atoms[0].implicit_h_override = Some(5);
        let decoded = parse_mol_at(&v3000.molfile(), vec2(0., 0.)).unwrap();
        assert_eq!(decoded.atoms[0].implicit_h_override, Some(5));
    }

    #[test]
    fn isotope_quick_input_edits_selected_atom_and_is_undoable() {
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "N");
        editor.selected = Some(0);
        editor.apply_quick(shortcuts::QuickAction::Isotope(Some(15)));
        assert_eq!(editor.mol.atoms[0].isotope, Some(15));
        editor.undo();
        assert_eq!(editor.mol.atoms[0].isotope, None);
        editor.redo();
        assert_eq!(editor.mol.atoms[0].isotope, Some(15));

        let exported = editor.molfile();
        let imported = parse_mol_at(&exported, vec2(0., 0.)).unwrap();
        assert_eq!(imported.atoms[0].element, "N");
        assert_eq!(imported.atoms[0].isotope, Some(15));

        editor.selected = Some(0);
        editor.apply_quick(shortcuts::QuickAction::Isotope(None));
        assert_eq!(editor.mol.atoms[0].isotope, None);
        editor.undo();
        assert_eq!(editor.mol.atoms[0].isotope, Some(15));
    }

    #[test]
    fn alias_quick_input_displays_on_carbon_and_is_undoable() {
        let mut editor = Editor::new();
        let carbon = editor.add_atom(vec2(0., 0.), "C");
        let oxygen = editor.add_atom(vec2(40., 0.), "O");
        editor.add_bond(carbon, oxygen, 1);
        editor.selected = Some(carbon);
        editor.apply_quick(shortcuts::QuickAction::Alias(Some("Ph".into())));
        assert_eq!(editor.mol.atoms[carbon].alias.as_deref(), Some("Ph"));
        assert!(show_label(&editor.mol, carbon));
        editor.undo();
        assert_eq!(editor.mol.atoms[carbon].alias, None);
        assert!(!show_label(&editor.mol, carbon));
        editor.redo();
        assert_eq!(editor.mol.atoms[carbon].alias.as_deref(), Some("Ph"));

        editor.selected = Some(carbon);
        editor.apply_quick(shortcuts::QuickAction::Alias(None));
        assert_eq!(editor.mol.atoms[carbon].alias, None);
    }

    #[test]
    fn ez_cycle_is_not_applied_to_non_double_bonds() {
        let mut editor = Editor::new();
        editor.add_atom(vec2(0., 0.), "C");
        editor.add_atom(vec2(40., 0.), "C");
        editor.add_bond(0, 1, 1);
        editor.cycle_double_bond_stereo(0);
        assert_eq!(editor.mol.bonds[0].stereo, None);
        assert!(editor.undo.is_empty());
    }

    #[test]
    fn mol_roundtrip_keeps_connectivity_charge_and_isotope() {
        let mut editor = Editor::new();
        let c = editor.add_atom(vec2(200., 200.), "C");
        let o = editor.add_atom(vec2(243., 200.), "O");
        editor.add_bond(c, o, 2);
        editor.mol.atoms[o].charge = -1;
        editor.mol.atoms[c].isotope = Some(13);
        editor.mol.atoms[c].radical = 2;
        let read = parse_mol_at(&editor.molfile(), vec2(300., 300.)).unwrap();
        assert_eq!(read.atoms.len(), 2);
        assert_eq!(read.bonds[0].order, 2);
        assert_eq!(read.atoms[1].charge, -1);
        assert_eq!(read.atoms[0].isotope, Some(13));
        assert_eq!(read.atoms[0].radical, 2);
    }

    #[test]
    fn mol_v2000_roundtrip_preserves_supported_hydrogen_count_overrides() {
        for hydrogen_count in 0..=4 {
            let mut editor = Editor::new();
            let mut carbon = Atom::new(vec2(250., 250.), "C");
            carbon.implicit_h_override = Some(hydrogen_count);
            editor.mol.atoms.push(carbon);

            let molfile = editor.molfile();
            let atom_line = molfile.lines().nth(4).unwrap();
            assert_eq!(
                atom_line
                    .split_whitespace()
                    .nth(7)
                    .unwrap()
                    .parse::<u8>()
                    .unwrap(),
                hydrogen_count + 1
            );
            let imported = parse_mol_at(&molfile, vec2(0., 0.)).unwrap();
            assert_eq!(
                imported.atoms[0].implicit_h_override,
                Some(hydrogen_count),
                "H{hydrogen_count}"
            );
            assert!(!editor
                .export_losses("test.mol")
                .contains(&"advanced atom properties"));
        }
    }

    #[test]
    fn mol_export_warns_for_hydrogen_overrides_outside_v2000_support() {
        let mut editor = Editor::new();
        let mut carbon = Atom::new(vec2(0., 0.), "C");
        carbon.implicit_h_override = Some(6);
        editor.mol.atoms.push(carbon);
        assert!(editor
            .export_losses("test.mol")
            .contains(&"advanced atom properties"));

        let mut large = Editor::new();
        large
            .mol
            .atoms
            .extend((0..1000).map(|_| Atom::new(vec2(0., 0.), "C")));
        large.mol.atoms[0].implicit_h_override = Some(2);
        assert!(!large
            .export_losses("test.mol")
            .contains(&"advanced atom properties"));
        large.mol.atoms[0].implicit_h_override = Some(6);
        assert!(large
            .export_losses("test.mol")
            .contains(&"advanced atom properties"));
    }

    #[test]
    fn mol_radical_records_support_multiple_atoms_and_reject_invalid_types() {
        let mol = "MOL\n  Rust\n\n  2  0  0  0  0  0            999 V2000\n       0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n       1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\nM  RAD  2   1   2   2   3\nM  END\n";
        let read = parse_mol_at(mol, vec2(0., 0.)).unwrap();
        assert_eq!(read.atoms[0].radical, 2);
        assert_eq!(read.atoms[1].radical, 3);

        let invalid = mol.replace("   2   3", "   2   4");
        assert!(parse_mol_at(&invalid, vec2(0., 0.)).is_err());
    }

    #[test]
    fn mol_import_rejects_non_finite_coordinates_and_invalid_bonds() {
        let bad_coordinate = "MOL\n  Rust\n\n  1  0  0  0  0  0            999 V2000\n       NaN    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\nM  END\n";
        assert!(parse_mol_at(bad_coordinate, vec2(0., 0.)).is_err());

        let invalid_order = "MOL\n  Rust\n\n  2  1  0  0  0  0            999 V2000\n       0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n       1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  0  0  0  0  0\nM  END\n";
        assert!(parse_mol_at(invalid_order, vec2(0., 0.)).is_err());
    }

    #[test]
    fn undo_redo_restores_graphics_and_molecule() {
        let mut editor = Editor::new();
        editor.checkpoint();
        editor.add_atom(vec2(0., 0.), "N");
        editor.graphics.push(Graphic::Plus { pos: vec2(20., 0.) });
        editor.undo();
        assert!(editor.mol.atoms.is_empty() && editor.graphics.is_empty());
        editor.redo();
        assert_eq!(editor.mol.atoms.len(), 1);
        assert_eq!(editor.graphics.len(), 1);
    }

    #[test]
    fn lossy_export_is_flagged() {
        let mut editor = Editor::new();
        editor.graphics.push(Graphic::Plus { pos: vec2(0., 0.) });
        assert!(editor.export_losses("test.mols").is_empty());
        assert!(editor
            .export_losses("test.mol")
            .contains(&"graphics, reaction arrows and text"));
        editor.mol.atoms.push(Atom::new(vec2(0., 0.), "C"));
        editor.mol.atoms[0].valence_override = Some(15);
        assert!(!editor
            .export_losses("test.mol")
            .contains(&"advanced atom properties"));
        assert!(molfile_for(&editor.mol).contains("V3000"));
        assert!(!editor
            .export_losses("test.mols")
            .contains(&"advanced atom properties"));
    }

    #[test]
    fn lossy_save_requires_confirmation_before_writing() {
        let path =
            std::env::temp_dir().join(format!("molecula-loss-check-{}.mol", std::process::id()));
        let mut editor = Editor::new();
        editor.add_atom(vec2(20., 20.), "C");
        editor.graphics.push(Graphic::Plus {
            pos: vec2(80., 20.),
        });
        let path_text = path.to_string_lossy().to_string();
        editor.save(&path_text);
        assert!(
            !path.exists(),
            "lossy file must not be written before confirmation"
        );
        assert_eq!(editor.pending_export.as_deref(), Some(path_text.as_str()));
        assert!(editor.status.contains("graphics, reaction arrows and text"));

        editor.save_confirmed(&path_text);
        assert!(path.exists());
        assert!(editor.status.starts_with("Saved "));
        let _ = fs::remove_file(path);
    }

    #[test]
    fn smiles_export_warns_about_document_data_that_the_notation_cannot_store() {
        let mut editor = Editor::new();
        editor.mol = smiles::parse("CCO").unwrap();
        let losses = editor.export_losses("ethanol.smi");
        assert!(losses.contains(&"2D/3D coordinates and drawing layout"));

        editor.mol.name = "Ethanol drawing".into();
        editor
            .mol
            .properties
            .push(("source".into(), "class handout".into()));
        let losses = editor.export_losses("ethanol.smiles");
        assert!(losses.contains(&"molecule name and data fields"));

        editor.mol.atoms[0].alias = Some("ethyl-end".into());
        let losses = editor.export_losses("ethanol.smi");
        assert!(losses.contains(&"query, alias, R-group, radical, valence or stereo annotations"));
    }

    #[test]
    fn mol_export_warning_distinguishes_preserved_stereo_from_lost_atom_data() {
        let mut editor = Editor::new();
        let mut carbon = Atom::new(vec2(0., 0.), "C");
        carbon.stereo = Some("@@".into());
        carbon.radical = 2;
        carbon.isotope = Some(13);
        editor.mol.atoms.push(carbon);
        editor.mol.atoms.push(Atom::new(vec2(40., 0.), "F"));
        editor.mol.bonds.push(Bond {
            a: 0,
            b: 1,
            order: 1,
            style: BondStyle::Wedge,
            stereo: None,
            query: None,
        });
        assert!(!editor
            .export_losses("test.mol")
            .contains(&"advanced atom properties"));
        assert!(!editor
            .export_losses("test.mol")
            .contains(&"bond styles or query data"));

        editor.mol.atoms[0].alias = Some("Me".into());
        assert!(editor
            .export_losses("test.sdf")
            .contains(&"advanced atom properties"));
    }

    #[test]
    fn aromatic_mol_roundtrip_restores_atom_aromaticity_without_loss_warning() {
        let mut editor = Editor::new();
        editor.insert_template("Benzene", vec2(300., 300.));
        let imported = parse_mol_at(&editor.molfile(), vec2(0., 0.)).unwrap();
        assert!(imported.atoms.iter().all(|atom| atom.aromatic));
        assert!(imported
            .bonds
            .iter()
            .all(|bond| bond.order == 4 && bond.style == BondStyle::Aromatic));
        assert!(!editor
            .export_losses("benzene.mol")
            .contains(&"advanced atom properties"));
        assert!(!editor
            .export_losses("benzene.mol")
            .contains(&"bond styles or query data"));
    }

    #[test]
    fn benzene_and_pyridine_templates_have_expected_graphs() {
        let mut benzene = Editor::new();
        benzene.insert_template("Benzene", vec2(300., 300.));
        assert_eq!((benzene.mol.atoms.len(), benzene.mol.bonds.len()), (6, 6));
        assert_eq!(chemistry::analyze(&benzene.mol).formula, "C6H6");
        let mut pyridine = Editor::new();
        pyridine.insert_template("Pyridine", vec2(300., 300.));
        assert_eq!(chemistry::analyze(&pyridine.mol).formula, "C5H5N");
    }

    #[test]
    fn explicit_hydrogen_add_remove_preserves_formula() {
        let mut e = Editor::new();
        let c1 = e.add_atom(vec2(0., 0.), "C");
        let c2 = e.add_atom(vec2(BOND, 0.), "C");
        let o = e.add_atom(vec2(BOND * 2., 0.), "O");
        e.add_bond(c1, c2, 1);
        e.add_bond(c2, o, 1);
        let formula = chemistry::analyze(&e.mol).formula;
        assert_eq!(formula, "C2H6O");
        e.add_explicit_hydrogens();
        assert_eq!(chemistry::analyze(&e.mol).formula, formula);
        assert_eq!(e.mol.atoms.iter().filter(|a| a.element == "H").count(), 6);
        e.remove_explicit_hydrogens();
        assert_eq!(chemistry::analyze(&e.mol).formula, formula);
    }

    #[test]
    fn mirror_inverts_stored_stereo_and_rotation_preserves_it() {
        let mut e = Editor::new();
        let i = e.add_atom(vec2(10., 20.), "C");
        e.mol.atoms[i].stereo = Some("R".into());
        e.transform_document(std::f32::consts::FRAC_PI_2, false, false);
        assert_eq!(e.mol.atoms[i].stereo.as_deref(), Some("R"));
        e.transform_document(0., true, false);
        assert_eq!(e.mol.atoms[i].stereo.as_deref(), Some("S"));
        e.mol.atoms[i].stereo = Some("@".into());
        e.transform_document(0., true, false);
        assert_eq!(e.mol.atoms[i].stereo.as_deref(), Some("@@"));
    }

    #[test]
    fn heterocycle_templates_are_real_aromatic_graphs() {
        let mut e = Editor::new();
        e.insert_template("Imidazole", vec2(300., 300.));
        assert_eq!((e.mol.atoms.len(), e.mol.bonds.len()), (5, 5));
        assert_eq!(e.mol.atoms.iter().filter(|a| a.element == "N").count(), 2);
        assert!(e.mol.atoms.iter().all(|a| a.aromatic));
        assert!(e.mol.bonds.iter().all(|b| b.style == BondStyle::Aromatic));
    }

    #[test]
    fn aromatic_heterocycle_templates_have_correct_graph_formulas() {
        for (template, expected_formula, hetero_element, expected_h) in [
            ("Imidazole", "C3H4N2", "N", 2usize),
            ("Pyrrole", "C4H5N", "N", 1),
            ("Furan", "C4H4O", "O", 1),
            ("Thiophene", "C4H4S", "S", 1),
        ] {
            let mut editor = Editor::new();
            editor.insert_template(template, vec2(300., 300.));
            assert_eq!(
                chemistry::analyze(&editor.mol).formula,
                expected_formula,
                "{template}"
            );
            assert_eq!(
                editor
                    .mol
                    .atoms
                    .iter()
                    .filter(|atom| atom.element == hetero_element)
                    .count(),
                expected_h,
                "{template} heteroatom count"
            );
            assert!(
                editor.mol.atoms.iter().all(|atom| atom.aromatic),
                "{template}"
            );
            assert!(
                editor
                    .mol
                    .bonds
                    .iter()
                    .all(|bond| bond.style == BondStyle::Aromatic),
                "{template}"
            );
            assert!(
                chemistry::analyze(&editor.mol).warnings.is_empty(),
                "{template}: {:?}",
                chemistry::analyze(&editor.mol).warnings
            );
        }
    }

    #[test]
    fn ion_charge_does_not_modify_unrelated_atoms() {
        let mut e = Editor::new();
        e.add_atom(vec2(0., 0.), "O");
        e.add_atom(vec2(40., 0.), "C");
        e.selected = Some(0);
        e.apply_quick(shortcuts::QuickAction::Charge(-1));
        assert_eq!(e.mol.atoms[0].charge, -1);
        assert_eq!(e.mol.atoms[1].charge, 0);
        e.undo();
        assert_eq!(e.mol.atoms[0].charge, 0);
    }

    #[test]
    fn every_top_level_menu_has_working_commands() {
        for (i, name) in MENU_NAMES.iter().enumerate() {
            assert!(!menu_items(i).is_empty(), "{} menu has no commands", name);
        }
        let calculation_commands = menu_items(7);
        assert!(calculation_commands.contains(&"Formula & mass"));
        assert!(calculation_commands.contains(&"Molecular properties"));
        assert!(calculation_commands.contains(&"Local topological descriptors"));
        assert!(calculation_commands.contains(&"Reaction analysis"));
        assert!(calculation_commands.contains(&"Reaction coefficient: custom"));
        assert!(calculation_commands.contains(&"Gasteiger partial charges"));
        assert!(calculation_commands.contains(&"View next reaction step"));
        assert!(calculation_commands.contains(&"View previous reaction step"));
        assert!(calculation_commands.contains(&"RDKit descriptors"));
        assert!(calculation_commands.contains(&"Enumerate tautomers"));
        assert!(calculation_commands.contains(&"Enumerate stereoisomers"));
        assert!(menu_items(6).contains(&"Find substructure matches"));
        assert!(menu_items(6).contains(&"Auto-map reaction (Complete)"));
        assert!(menu_items(6).contains(&"Auto-map reaction (Matching)"));
        assert!(menu_items(6).contains(&"Auto-map reaction (Changing)"));
        let mut editor = Editor::new();
        editor.menu_action("Recent files");
        assert!(editor.recent_open);
    }

    #[test]
    fn automatic_atom_mapping_modes_apply_atomically_and_undo() {
        let mut e = Editor::new();
        for (i, element) in ["C", "O", "C", "O", "N"].into_iter().enumerate() {
            e.mol
                .atoms
                .push(Atom::new(vec2(i as f32 * 30., 0.), element));
        }
        e.mol.bonds.push(Bond {
            a: 0,
            b: 1,
            order: 1,
            style: BondStyle::Normal,
            stereo: None,
            query: None,
        });
        e.mol.bonds.push(Bond {
            a: 2,
            b: 3,
            order: 2,
            style: BondStyle::Normal,
            stereo: None,
            query: None,
        });
        e.mol.atoms[0].atom_map = Some(7);
        let before = e.mol.clone();
        let report = "STATUS\tcomplete\nOPTIMALITY\toptimal_given_seed\nENGINE\tAtomMapper test\nPAIR\t0\t0\t0\nPAIR\t1\t1\t1\n";
        e.apply_auto_atom_mapping(report, &[0, 1], &[2, 3], "Matching");
        assert_eq!(e.mol.atoms[0].atom_map, Some(1));
        assert_eq!(e.mol.atoms[2].atom_map, Some(1));
        assert_eq!(e.mol.atoms[1].atom_map, Some(2));
        assert_eq!(e.mol.atoms[3].atom_map, Some(2));
        assert_eq!(e.mol.atoms[4].atom_map, None);
        e.undo();
        assert_eq!(e.mol, before);

        e.apply_auto_atom_mapping(report, &[0, 1], &[2, 3], "Changing");
        assert_eq!(e.mol.atoms[0].atom_map, None);
        assert_eq!(e.mol.atoms[2].atom_map, None);
        assert_eq!(e.mol.atoms[1].atom_map, Some(1));
        assert_eq!(e.mol.atoms[3].atom_map, Some(1));
        e.undo();
        assert_eq!(e.mol, before);
    }

    #[test]
    fn incomplete_or_out_of_range_auto_mapping_never_changes_graph() {
        let mut e = Editor::new();
        for i in 0..2 {
            e.mol.atoms.push(Atom::new(vec2(i as f32 * 30., 0.), "C"));
        }
        e.mol.atoms[0].atom_map = Some(12);
        let before = e.mol.clone();
        e.apply_auto_atom_mapping("STATUS\ttimed_out\n", &[0], &[1], "Complete");
        assert_eq!(e.mol, before);
        e.apply_auto_atom_mapping("STATUS\tcomplete\nPAIR\t9\t0\t0\n", &[0], &[1], "Complete");
        assert_eq!(e.mol, before);
        e.apply_auto_atom_mapping(
            "STATUS\tcomplete\nPAIR\t0\t0\t0\nPAIR\tbad\t1\t0\n",
            &[0],
            &[1],
            "Complete",
        );
        assert_eq!(e.mol, before);
        e.apply_auto_atom_mapping("STATUS\tcomplete\nPAIR\t0\t0\t2\n", &[0], &[1], "Complete");
        assert_eq!(e.mol, before);
        e.apply_auto_atom_mapping("STATUS\tcomplete\nPAIR\t0\t0\t0\n", &[0], &[1], "Changing");
        assert_eq!(e.mol, before);
        assert!(e.undo.is_empty());
        e.apply_auto_atom_mapping(
            "STATUS\tcomplete\nPAIR\t0\t0\t1\nPAIR\t0\t0\t1\n",
            &[0],
            &[1],
            "Changing",
        );
        assert_eq!(e.mol, before);
        assert!(e.undo.is_empty());
    }

    #[test]
    fn overlapping_background_chemistry_jobs_do_not_replace_mapping_results() {
        let mut e = Editor::new();
        let (tx, rx) = mpsc::channel();
        e.calc_rx = Some((e.mol.clone(), rx));
        e.menu_action("Auto-map reaction (Complete)");
        assert!(e.calc_rx.is_some());
        assert!(e.auto_map_pending.is_none());
        assert!(e.status.contains("Wait for the current background"));
        let _ = tx.send(Ok("original worker result".into()));
    }

    #[test]
    fn selected_multicomponent_reaction_runs_optional_mapper_and_maps_undoably() {
        let python = std::env::var("MOLECULA_PYTHON").unwrap_or_else(|_| "python3".into());
        let available = Command::new(python)
            .args(["-c", "import rdkit, atommapper"])
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|status| status.success());
        if !available {
            return;
        }

        let mut e = Editor::new();
        for (i, element) in ["C", "C", "O", "C", "C", "O"].into_iter().enumerate() {
            let mut atom = Atom::new(vec2(i as f32 * 28., 0.), element);
            atom.reaction_role = Some(if i < 3 { "reactant" } else { "product" }.into());
            atom.reaction_step = 1;
            e.mol.atoms.push(atom);
        }
        // Each reaction side has an ethane component plus a disconnected oxygen.
        for (a, b, order) in [(0, 1, 1), (3, 4, 1)] {
            e.mol.bonds.push(Bond {
                a,
                b,
                order,
                style: BondStyle::Normal,
                stereo: None,
                query: None,
            });
        }
        e.selected_atoms = (0..6).collect();
        e.start_auto_atom_mapping("Complete");
        let (captured, rx) = e.calc_rx.take().expect("mapping job should be queued");
        let report = rx
            .recv_timeout(std::time::Duration::from_secs(40))
            .expect("mapping worker should finish")
            .expect("optional mapper should complete");
        assert_eq!(captured, e.mol);
        e.apply_auto_atom_mapping(&report, &[0, 1, 2], &[3, 4, 5], "Complete");
        assert!(e.mol.atoms[..3].iter().all(|atom| atom.atom_map.is_some()));
        assert!(e.mol.atoms[3..].iter().all(|atom| atom.atom_map.is_some()));
        assert_eq!(e.mol.atoms[0].atom_map, e.mol.atoms[3].atom_map);
        assert_eq!(e.mol.atoms[1].atom_map, e.mol.atoms[4].atom_map);
        assert_eq!(e.mol.atoms[2].atom_map, e.mol.atoms[5].atom_map);
        let mapped = e.mol.clone();
        e.undo();
        assert!(e.mol.atoms.iter().all(|atom| atom.atom_map.is_none()));
        e.redo();
        assert_eq!(e.mol, mapped);
    }

    #[test]
    fn automatic_mapping_refuses_selected_agents_and_cross_step_components() {
        let mut editor = Editor::new();
        for (role, step) in [("reactant", 1), ("product", 1), ("agent", 1)] {
            let mut atom = Atom::new(vec2(editor.mol.atoms.len() as f32 * 40., 0.), "C");
            atom.reaction_role = Some(role.into());
            atom.reaction_step = step;
            editor.mol.atoms.push(atom);
        }
        editor.selected_atoms = vec![0, 1, 2];
        let before = editor.mol.clone();
        editor.start_auto_atom_mapping("Complete");
        assert_eq!(editor.mol, before);
        assert!(editor.calc_rx.is_none());
        assert!(editor.auto_map_pending.is_none());
        assert!(editor
            .info_overlay
            .as_deref()
            .is_some_and(|message| message.contains("reactant or product")));

        editor.mol.atoms[2].reaction_role = Some("product".into());
        editor.mol.atoms[2].reaction_step = 2;
        editor.selected_atoms = vec![0, 1, 2];
        let before = editor.mol.clone();
        editor.start_auto_atom_mapping("Complete");
        assert_eq!(editor.mol, before);
        assert!(editor.calc_rx.is_none());
        assert!(editor.auto_map_pending.is_none());
        assert!(editor
            .info_overlay
            .as_deref()
            .is_some_and(|message| message.contains("same reaction step")));
    }

    #[test]
    fn reaction_step_navigation_is_not_limited_to_three_steps() {
        let mut editor = Editor::new();
        for _ in 0..5 {
            editor.graphics.push(Graphic::Arrow {
                start: vec2(200., 100.),
                end: vec2(300., 100.),
                kind: "reaction".into(),
                above: String::new(),
                below: String::new(),
            });
        }
        editor.menu_action("View next reaction step");
        assert_eq!(editor.reaction_step_view, Some(2));
        for _ in 0..5 {
            editor.menu_action("View next reaction step");
        }
        assert_eq!(editor.reaction_step_view, Some(5));
        editor.menu_action("View previous reaction step");
        assert_eq!(editor.reaction_step_view, Some(4));
        editor.menu_action("View all reaction steps");
        assert_eq!(editor.reaction_step_view, None);
    }

    #[test]
    fn topological_descriptor_command_displays_computed_method_and_values() {
        let mut editor = Editor::new();
        let carbon = editor.add_atom(vec2(0., 0.), "C");
        let oxygen = editor.add_atom(vec2(40., 0.), "O");
        editor.add_bond(carbon, oxygen, 1);
        editor.menu_action("Local topological descriptors");
        let report = editor.info_overlay.unwrap();
        assert!(report.contains("Wiener index: 1"));
        assert!(report.contains("First Zagreb index (M₁): 2"));
        assert!(report.contains("Randić connectivity index"));
        assert!(report.contains("explicit hydrogens excluded"));
    }

    #[test]
    fn substructure_command_reports_and_selects_actual_target_match() {
        let mut editor = Editor::new();
        for element in ["C", "O", "C", "C", "O"] {
            editor.add_atom(vec2(0., 0.), element);
        }
        for (a, b) in [(0, 1), (2, 3), (3, 4)] {
            editor.mol.bonds.push(Bond {
                a,
                b,
                order: 1,
                style: BondStyle::Normal,
                stereo: None,
                query: None,
            });
        }
        editor.selected_atoms = vec![0, 1];
        editor.selected = Some(0);
        editor.menu_action("Find substructure matches");
        assert_eq!(editor.selected_atoms, vec![3, 4]);
        assert!(editor
            .info_overlay
            .as_deref()
            .unwrap()
            .contains("1 match(es)"));
        assert!(editor
            .info_overlay
            .as_deref()
            .unwrap()
            .contains("target atoms 4, 5"));
    }

    #[test]
    fn long_chemical_identifiers_wrap_without_losing_characters() {
        let smiles = "C".repeat(171);
        let wrapped = wrap_info_lines(&format!("CANDIDATE 1 {smiles}"), 32);
        assert!(wrapped.len() > 1);
        assert_eq!(wrapped.concat(), format!("CANDIDATE 1 {smiles}"));
    }

    #[test]
    fn calculation_menu_commands_return_real_results_or_clear_unavailable_states() {
        let mut e = Editor::new();
        e.add_atom(vec2(0., 0.), "C");
        e.menu_action("Formula & mass");
        assert!(e.info_overlay.as_deref().unwrap().contains("CH4"));
        e.menu_action("Molecular properties");
        assert!(e.info_overlay.as_deref().unwrap().contains("Atoms: 1"));
        e.menu_action("Reaction analysis");
        assert!(e
            .info_overlay
            .as_deref()
            .unwrap()
            .contains("Draw a reaction"));
        e.menu_action("RDKit descriptors");
        assert!(e
            .info_overlay
            .as_deref()
            .unwrap()
            .contains("Calculating RDKit"));
        assert!(e.calc_rx.is_some());
    }

    #[test]
    fn rdkit_candidate_actions_run_as_optional_background_jobs() {
        let mut e = Editor::new();
        e.add_atom(vec2(0., 0.), "C");
        e.menu_action("Enumerate tautomers");
        assert!(e.calc_rx.is_some());
        assert!(e
            .info_overlay
            .as_deref()
            .unwrap()
            .contains("tautomer-report"));
        e.menu_action("Enumerate stereoisomers");
        assert!(e.calc_rx.is_some());
        assert!(e.status.contains("Wait for the current background"));
        assert!(e
            .info_overlay
            .as_deref()
            .unwrap()
            .contains("tautomer-report"));
    }

    #[test]
    fn gasteiger_command_runs_in_background_without_mutating_the_molecule() {
        let python = std::env::var("MOLECULA_PYTHON").unwrap_or_else(|_| "python3".into());
        let available = Command::new(python)
            .args(["-c", "import rdkit"])
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|status| status.success());
        if !available {
            return;
        }
        let mut e = Editor::new();
        for (i, element) in ["C", "C", "O"].into_iter().enumerate() {
            e.add_atom(vec2(i as f32 * BOND, 0.), element);
        }
        e.add_bond(0, 1, 1);
        e.add_bond(1, 2, 1);
        let before = e.mol.clone();
        e.menu_action("Gasteiger partial charges");
        let (captured, rx) = e
            .calc_rx
            .take()
            .expect("calculation should be queued in background");
        let report = rx
            .recv_timeout(std::time::Duration::from_secs(30))
            .unwrap()
            .unwrap();
        assert_eq!(captured, before);
        assert_eq!(e.mol, before);
        assert!(report.contains("METHOD\tGasteiger-Marsili"));
        assert!(report.contains("SUMMARY\tNet model charge"));
        assert_eq!(
            report
                .lines()
                .filter(|line| line.starts_with("ATOM\t"))
                .count(),
            3
        );
    }

    #[test]
    fn drag_selection_finds_atoms_in_a_rectangle_in_either_direction() {
        let mut e = Editor::new();
        e.add_atom(vec2(20., 30.), "C");
        e.add_atom(vec2(80., 90.), "O");
        e.add_atom(vec2(200., 200.), "N");
        assert_eq!(
            e.atoms_in_rect(vec2(10., 10.), vec2(100., 100.)),
            vec![0, 1]
        );
        e.select_box(vec2(100., 100.), vec2(10., 10.));
        assert_eq!(e.selected_atoms, vec![0, 1]);
        assert_eq!(e.selected, None);
        assert!(e.status.contains("2 atom(s)"));
    }

    #[test]
    fn click_selects_a_bond_for_the_inspector_and_delete_is_undoable() {
        let mut e = Editor::new();
        e.add_atom(vec2(100., 100.), "C");
        e.add_atom(vec2(143., 100.), "O");
        e.add_bond(0, 1, 2);
        e.select_at(vec2(121.5, 100.), false);
        assert_eq!(e.selected_bond, Some(0));
        assert!(e.selected.is_none());

        e.menu_action("Delete selection");
        assert!(e.mol.bonds.is_empty());
        e.undo();
        assert_eq!(e.mol.bonds.len(), 1);
    }

    #[test]
    fn zoom_keeps_cursor_anchor_and_pan_inverse_are_consistent() {
        let mut e = Editor::new();
        let anchor = vec2(640., 360.);
        let before = e.to_world(anchor, 1280., 720.);
        e.adjust_zoom(1.5, anchor, 1280., 720.);
        let after = e.to_world(anchor, 1280., 720.);
        assert!(before.distance(after) < 0.001);
        e.pan += vec2(23., -11.);
        let world = e.to_world(vec2(700., 400.), 1280., 720.);
        let screen = e.to_screen(world, 1280., 720.);
        assert!(screen.distance(vec2(700., 400.)) < 0.001);
    }

    #[test]
    fn selection_can_move_duplicate_rotate_and_delete_graph_and_graphics() {
        let mut e = Editor::new();
        let a = e.add_atom(vec2(100., 100.), "C");
        let b = e.add_atom(vec2(143., 100.), "O");
        e.add_bond(a, b, 1);
        e.graphics.push(Graphic::Plus {
            pos: vec2(100., 200.),
        });
        e.select_box(vec2(80., 80.), vec2(160., 120.));
        assert_eq!(e.selected_atoms, vec![0, 1]);
        assert_eq!(e.selected_bonds, vec![0]);
        e.move_selection(vec2(10., 20.));
        assert_eq!(e.mol.atoms[0].pos, vec2(110., 120.));
        assert_eq!(e.mol.atoms[1].pos, vec2(153., 120.));
        assert_eq!(
            e.graphics[0],
            Graphic::Plus {
                pos: vec2(100., 200.)
            }
        );
        e.select_box(vec2(80., 80.), vec2(170., 150.));
        e.duplicate_selection();
        assert_eq!(e.mol.atoms.len(), 4);
        assert_eq!(e.mol.bonds.len(), 2);
        assert!(e.mol.atoms[2].atom_map.is_none());
        e.rotate_selection(std::f32::consts::FRAC_PI_2);
        e.delete_selection();
        assert_eq!(e.mol.atoms.len(), 2);
        assert_eq!(e.mol.bonds.len(), 1);
        e.undo();
        assert_eq!(e.mol.atoms.len(), 4);
    }

    #[test]
    fn reaction_component_commands_are_undoable_and_apply_to_whole_fragment() {
        let mut e = Editor::new();
        let a = e.add_atom(vec2(0., 0.), "C");
        let b = e.add_atom(vec2(40., 0.), "O");
        e.add_bond(a, b, 1);
        e.selected = Some(a);
        e.selected_atoms = vec![a];
        e.menu_action("Reaction component: agent");
        assert_eq!(e.mol.atoms[0].reaction_role.as_deref(), Some("agent"));
        assert_eq!(e.mol.atoms[1].reaction_role.as_deref(), Some("agent"));
        e.menu_action("Reaction coefficient: 3");
        assert_eq!(
            e.mol
                .atoms
                .iter()
                .map(|x| x.coefficient)
                .collect::<Vec<_>>(),
            vec![3, 3]
        );
        e.undo();
        assert_eq!(e.mol.atoms[0].coefficient, 1);
    }

    #[test]
    fn custom_reaction_coefficient_edits_component_and_validates_range() {
        let mut e = Editor::new();
        let a = e.add_atom(vec2(0., 0.), "C");
        let b = e.add_atom(vec2(40., 0.), "O");
        let c = e.add_atom(vec2(200., 0.), "C");
        let d = e.add_atom(vec2(240., 0.), "O");
        e.add_bond(a, b, 1);
        e.add_bond(c, d, 1);
        e.selected_atoms = vec![a];
        e.menu_action("Reaction coefficient: custom");
        assert!(matches!(
            e.text_target,
            Some(TextTarget::ReactionCoefficient)
        ));
        e.text_buffer = "42".into();
        e.commit_text();
        assert_eq!(
            e.mol
                .atoms
                .iter()
                .map(|atom| atom.coefficient)
                .collect::<Vec<_>>(),
            vec![42, 42, 1, 1]
        );
        e.undo();
        assert!(e.mol.atoms.iter().all(|atom| atom.coefficient == 1));

        e.selected_atoms = vec![a];
        e.text_target = Some(TextTarget::ReactionCoefficient);
        e.text_buffer = "100".into();
        e.commit_text();
        assert!(e.mol.atoms.iter().all(|atom| atom.coefficient == 1));
        assert!(e.status.contains("from 1 to 99"));
    }

    #[test]
    fn automatic_reaction_balance_updates_coefficients_and_undo_restores_them() {
        let mut e = Editor::new();
        for (element, x, role) in [
            ("H", 80., "reactant"),
            ("H", 90., "reactant"),
            ("O", 120., "reactant"),
            ("O", 130., "reactant"),
            ("O", 420., "product"),
            ("H", 410., "product"),
            ("H", 430., "product"),
        ] {
            let mut atom = Atom::new(vec2(x, 100.), element);
            atom.reaction_role = Some(role.into());
            e.mol.atoms.push(atom);
        }
        for (a, b, order) in [(0, 1, 1), (2, 3, 2), (4, 5, 1), (4, 6, 1)] {
            e.mol.bonds.push(Bond {
                a,
                b,
                order,
                style: BondStyle::Normal,
                stereo: None,
                query: None,
            });
        }
        e.graphics.push(Graphic::Arrow {
            start: vec2(220., 100.),
            end: vec2(300., 100.),
            kind: "reaction".into(),
            above: String::new(),
            below: String::new(),
        });
        e.menu_action("Auto-balance reaction");
        assert_eq!(e.mol.atoms[0].coefficient, 2);
        assert_eq!(e.mol.atoms[2].coefficient, 1);
        assert_eq!(e.mol.atoms[4].coefficient, 2);
        assert_eq!(
            reaction::analyze_step(&e.mol, &e.graphics, 1)
                .unwrap()
                .balanced,
            Some(true)
        );
        e.undo();
        assert!(e.mol.atoms.iter().all(|atom| atom.coefficient == 1));
    }

    #[test]
    fn automatic_balance_all_steps_is_atomic_and_undoable() {
        let mut e = Editor::new();
        for step in 1..=2 {
            let y = step as f32 * 220.0;
            let offset = e.mol.atoms.len();
            for (element, x, role) in [
                ("H", 80., "reactant"),
                ("H", 90., "reactant"),
                ("O", 120., "reactant"),
                ("O", 130., "reactant"),
                ("O", 420., "product"),
                ("H", 410., "product"),
                ("H", 430., "product"),
            ] {
                let mut atom = Atom::new(vec2(x, y), element);
                atom.reaction_role = Some(role.into());
                atom.reaction_step = step;
                e.mol.atoms.push(atom);
            }
            for (a, b, order) in [(0, 1, 1), (2, 3, 2), (4, 5, 1), (4, 6, 1)] {
                e.mol.bonds.push(Bond {
                    a: offset + a,
                    b: offset + b,
                    order,
                    style: BondStyle::Normal,
                    stereo: None,
                    query: None,
                });
            }
            e.graphics.push(Graphic::Arrow {
                start: vec2(220., y),
                end: vec2(300., y),
                kind: "reaction".into(),
                above: String::new(),
                below: String::new(),
            });
        }
        e.menu_action("Auto-balance all reaction steps");
        assert_eq!(e.mol.atoms[0].coefficient, 2);
        assert_eq!(e.mol.atoms[7].coefficient, 2);
        assert!(e
            .info_overlay
            .as_deref()
            .unwrap()
            .contains("all 2 reaction steps"));
        e.undo();
        assert!(e.mol.atoms.iter().all(|atom| atom.coefficient == 1));

        // If a later step cannot be solved, an earlier valid step stays untouched.
        e.graphics.push(Graphic::Arrow {
            start: vec2(220., 700.),
            end: vec2(300., 700.),
            kind: "reaction".into(),
            above: String::new(),
            below: String::new(),
        });
        e.menu_action("Auto-balance all reaction steps");
        assert!(e.mol.atoms.iter().all(|atom| atom.coefficient == 1));
        assert!(e
            .info_overlay
            .as_deref()
            .unwrap()
            .contains("No coefficients were changed"));
    }

    #[test]
    fn save_open_and_export_workflow_uses_real_files() {
        let base = std::env::temp_dir().join(format!(
            "molecula-e2e-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let mols = base.with_extension("mols");
        let mol = base.with_extension("mol");
        let cml_path = base.with_extension("cml");
        let svg = base.with_extension("svg");

        let mut original = Editor::new();
        let c = original.add_atom(vec2(300., 300.), "C");
        let o = original.add_atom(vec2(343., 300.), "O");
        original.add_bond(c, o, 2);
        original.mol.atoms[c].isotope = Some(13);
        original.mol.atoms[o].charge = -1;
        original.graphics.push(Graphic::Arrow {
            start: vec2(300., 250.),
            end: vec2(350., 250.),
            kind: "reaction".into(),
            above: "heat".into(),
            below: "".into(),
        });

        original.save(mols.to_str().unwrap());
        assert!(original.status.starts_with("Saved "), "{}", original.status);
        let mut reopened = Editor::new();
        reopened.load(mols.to_str().unwrap());
        assert_eq!(reopened.mol, original.mol);
        assert_eq!(reopened.graphics, original.graphics);

        reopened.save(mol.to_str().unwrap());
        assert!(reopened.pending_export.is_some());
        reopened.save_confirmed(mol.to_str().unwrap());
        assert!(reopened.status.starts_with("Saved "), "{}", reopened.status);
        let mol_data = fs::read_to_string(&mol).unwrap();
        let imported = parse_mol_at(&mol_data, vec2(300., 300.)).unwrap();
        assert_eq!(imported.atoms.len(), 2);
        assert_eq!(imported.bonds.len(), 1);
        assert_eq!(imported.atoms[0].isotope, Some(13));
        assert_eq!(imported.atoms[1].charge, -1);

        let mut cml_editor = Editor::new();
        cml_editor.mol = reopened.mol.clone();
        cml_editor.graphics.clear();
        cml_editor.save(cml_path.to_str().unwrap());
        assert!(
            cml_editor.status.starts_with("Exported "),
            "{}",
            cml_editor.status
        );
        let cml_data = fs::read_to_string(&cml_path).unwrap();
        let cml_reopened = cml::parse(&cml_data, vec2(300., 300.)).unwrap();
        assert_eq!(cml_reopened.atoms.len(), 2);
        assert_eq!(cml_reopened.bonds.len(), 1);
        assert_eq!(cml_reopened.atoms[0].isotope, Some(13));
        assert_eq!(cml_reopened.atoms[1].charge, -1);

        reopened
            .write_svg(svg.to_str().unwrap(), 1920., 1080.)
            .unwrap();
        let svg_data = fs::read_to_string(&svg).unwrap();
        assert!(svg_data.starts_with("<svg "));
        assert!(svg_data.contains("<line "));
        assert!(svg_data.contains(">heat</text>"));

        for path in [mols, mol, cml_path, svg] {
            let _ = fs::remove_file(path);
        }
    }

    #[test]
    fn sdf_records_and_properties_roundtrip() {
        let mut mol = Molecule::default();
        mol.atoms.push(Atom::new(vec2(188., 218.), "O"));
        mol.properties
            .push(("SOURCE".into(), "teaching sample".into()));
        let data = format!(
            "{}>  <SOURCE>\nteaching sample\n\n$$$$\n",
            molfile_for(&mol)
        );
        let loaded = parse_sdf(&data).unwrap();
        assert_eq!(loaded.len(), 1);
        assert_eq!(loaded[0].atoms.len(), 1);
        assert_eq!(loaded[0].atoms[0].element, "O");
        assert_eq!(
            loaded[0].properties,
            vec![("SOURCE".into(), "teaching sample".into())]
        );
    }

    #[test]
    fn undo_history_is_isolated_per_document_tab() {
        let mut e = Editor::new();
        e.checkpoint();
        e.add_atom(vec2(10., 10.), "C");
        e.new_document();
        e.checkpoint();
        e.add_atom(vec2(20., 20.), "O");
        e.switch_document(0);
        assert_eq!(e.mol.atoms.len(), 1);
        e.undo();
        assert!(e.mol.atoms.is_empty());
        e.switch_document(1);
        assert_eq!(e.mol.atoms.len(), 1);
        assert_eq!(e.mol.atoms[0].element, "O");
        e.undo();
        assert!(e.mol.atoms.is_empty());
    }

    #[test]
    fn undo_redo_remains_consistent_at_history_limit() {
        let mut e = Editor::new();
        for i in 0..125 {
            e.checkpoint();
            e.add_atom(vec2(i as f32 * 2., 0.), "C");
        }
        assert_eq!(e.mol.atoms.len(), 125);
        for _ in 0..100 {
            e.undo();
        }
        assert_eq!(e.mol.atoms.len(), 25);
        for _ in 0..100 {
            e.redo();
        }
        assert_eq!(e.mol.atoms.len(), 125);
        assert_eq!(chemistry::analyze(&e.mol).atom_count, 125);
    }
}
