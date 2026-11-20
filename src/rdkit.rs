use std::io::Write;
use std::process::{Command, Stdio};

const SCRIPT: &str = include_str!("../scripts/rdkit_descriptors.py");
const CONVERT_SCRIPT: &str = include_str!("../scripts/rdkit_convert.py");
const MAPPING_SCRIPT: &str = include_str!("../scripts/atom_mapping.py");
const CHARGES_SCRIPT: &str = include_str!("../scripts/rdkit_charges.py");

fn execute(script: &str, operation: &str, input_text: &str) -> Result<String, String> {
    let candidates = if let Ok(path) = std::env::var("MOLECULA_PYTHON") {
        vec![path]
    } else {
        vec!["python3.12".into(), "python3".into()]
    };
    let mut failure =
        "No Python interpreter found. Set MOLECULA_PYTHON to a Python with RDKit installed."
            .to_string();
    for python in candidates {
        let child = Command::new(&python)
            .arg("-c")
            .arg(script)
            .arg(operation)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn();
        let Ok(mut child) = child else { continue };
        if let Some(mut input) = child.stdin.take() {
            input
                .write_all(input_text.as_bytes())
                .map_err(|e| e.to_string())?;
        }
        let output = child.wait_with_output().map_err(|e| e.to_string())?;
        let body = String::from_utf8(output.stdout).map_err(|_| "RDKit returned invalid UTF-8")?;
        if output.status.success() {
            return Ok(body);
        }
        failure = body
            .strip_prefix("ERROR\t")
            .unwrap_or("RDKit conversion failed")
            .trim()
            .to_string();
    }
    Err(failure)
}

pub fn from_smiles_with_stereo(smiles: &str) -> Result<(String, Vec<(usize, String)>), String> {
    let output = execute(CONVERT_SCRIPT, "from-smiles", smiles)?;
    let mut molblock = String::new();
    let mut stereo = Vec::new();
    for line in output.lines() {
        if let Some(record) = line.strip_prefix("BONDSTEREO\t") {
            let mut fields = record.split('\t');
            let index = fields
                .next()
                .and_then(|value| value.parse::<usize>().ok())
                .filter(|value| *value > 0)
                .ok_or("RDKit returned an invalid stereo bond index")?;
            let value = fields
                .next()
                .filter(|value| matches!(*value, "E" | "Z"))
                .ok_or("RDKit returned an invalid alkene stereo value")?;
            if fields.next().is_some() {
                return Err("RDKit returned a malformed stereo record".into());
            }
            stereo.push((index - 1, value.to_string()));
        } else {
            molblock.push_str(line);
            molblock.push('\n');
        }
    }
    if molblock.trim().is_empty() {
        return Err("RDKit returned no MOL structure".into());
    }
    Ok((molblock, stereo))
}

pub fn to_smiles(molblock: &str) -> Result<String, String> {
    execute(CONVERT_SCRIPT, "to-smiles", molblock).map(|s| s.trim().to_string())
}

pub fn to_smiles_with_stereo(molblock: &str, molecule: &crate::Molecule) -> Result<String, String> {
    let mut input = String::new();
    let mut has_stereo = false;
    for (index, bond) in molecule.bonds.iter().enumerate() {
        if let Some(stereo) = bond.stereo.as_deref() {
            match stereo {
                "E" | "Z" => {
                    has_stereo = true;
                    input.push_str(&format!("BONDSTEREO\t{}\t{}\n", index + 1, stereo));
                }
                _ => return Err(format!("Unsupported bond stereo annotation: {stereo}")),
            }
        }
    }
    if !has_stereo {
        return to_smiles(molblock);
    }
    input.push_str(molblock);
    execute(CONVERT_SCRIPT, "to-smiles", &input).map(|s| s.trim().to_string())
}

pub fn generate_3d(molblock: &str) -> Result<String, String> {
    execute(CONVERT_SCRIPT, "generate-3d", molblock)
}
pub fn aromatize(molblock: &str) -> Result<String, String> {
    execute(CONVERT_SCRIPT, "aromatize", molblock)
}
pub fn kekulize(molblock: &str) -> Result<String, String> {
    execute(CONVERT_SCRIPT, "kekulize", molblock)
}
pub fn stereo_report(molblock: &str) -> Result<String, String> {
    execute(CONVERT_SCRIPT, "stereo-report", molblock)
}

/// Map all selected reactant-side components against all product-side
/// components with optional AtomMapper. Disconnected structures are supported;
/// result indices address the input side MOL graphs in their atom order.
pub fn map_reaction_pair(reactant: &str, product: &str) -> Result<String, String> {
    let input = format!("{reactant}\n$$$$\n{product}");
    execute(MAPPING_SCRIPT, "map-pair", &input).map_err(|error| {
        if error.contains("No module named") {
            "Automatic atom mapping requires optional Python packages RDKit and AtomMapper (MIT). Install them in the Python environment selected by MOLECULA_PYTHON. Manual atom mapping remains available.".into()
        } else {
            error
        }
    })
}

/// Calculate RDKit Gasteiger-Marsili partial charges for each explicitly drawn
/// atom, with the implicit hydrogen charge total shown separately per atom.
pub fn partial_charge_report(molblock: &str) -> Result<String, String> {
    execute(CHARGES_SCRIPT, "gasteiger", molblock)
}

pub fn import_external(operation: &str, input: &str) -> Result<String, String> {
    execute(CONVERT_SCRIPT, operation, input)
}

pub fn export_external(operation: &str, molblock: &str) -> Result<String, String> {
    let output = execute(CONVERT_SCRIPT, operation, molblock)?;
    let (_, body) = output
        .split_once("RESULT\n")
        .ok_or("RDKit returned no export result")?;
    Ok(body.to_string())
}

pub fn descriptors(molblock: &str) -> Result<String, String> {
    let candidates = if let Ok(path) = std::env::var("MOLECULA_PYTHON") {
        vec![path]
    } else {
        vec!["python3.12".into(), "python3".into()]
    };
    let mut failure =
        "No Python interpreter found. Set MOLECULA_PYTHON to a Python with RDKit installed."
            .to_string();
    for python in candidates {
        let child = Command::new(&python)
            .arg("-c")
            .arg(SCRIPT)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::null())
            .spawn();
        let Ok(mut child) = child else { continue };
        if let Some(mut input) = child.stdin.take() {
            input
                .write_all(molblock.as_bytes())
                .map_err(|e| e.to_string())?;
        }
        let output = child.wait_with_output().map_err(|e| e.to_string())?;
        let body = String::from_utf8(output.stdout).map_err(|_| "RDKit returned invalid UTF-8")?;
        if output.status.success() {
            let mut lines = Vec::new();
            for line in body.lines() {
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() == 2 && fields[0] == "ENGINE" {
                    lines.push(fields[1].to_string());
                } else if fields.len() == 5 && fields[0] == "VALUE" {
                    lines.push(format!(
                        "{}: {} {} ({})",
                        fields[1], fields[2], fields[3], fields[4]
                    ));
                } else if fields.len() == 3 && fields[0] == "UNAVAILABLE" {
                    lines.push(format!("{}: unavailable ({})", fields[1], fields[2]));
                }
            }
            if lines.len() < 2 {
                return Err("RDKit returned an incomplete descriptor report".into());
            }
            return Ok(lines.join("\n"));
        }
        failure = body
            .strip_prefix("ERROR\t")
            .unwrap_or("RDKit descriptor worker failed")
            .trim()
            .to_string();
    }
    Err(failure)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn has_rdkit(python: &str) -> bool {
        Command::new(python)
            .args(["-c", "import rdkit"])
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|status| status.success())
    }

    fn has_rdkit_freesasa(python: &str) -> bool {
        Command::new(python)
            .args(["-c", "from rdkit.Chem import rdFreeSASA"])
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|status| status.success())
    }

    #[test]
    fn installed_rdkit_workers_roundtrip_real_structures() {
        let explicit = std::env::var("MOLECULA_PYTHON").ok();
        let python = explicit.clone().or_else(|| {
            ["python3.12", "python3"]
                .into_iter()
                .find(|candidate| has_rdkit(candidate))
                .map(str::to_owned)
        });
        let Some(python) = python else { return };
        assert!(
            has_rdkit(&python),
            "MOLECULA_PYTHON must point to a Python environment with RDKit installed"
        );

        let mut molecule = crate::Molecule::default();
        for (index, element) in ["C", "C", "O"].into_iter().enumerate() {
            molecule.atoms.push(crate::Atom::new(
                macroquad::prelude::vec2(index as f32 * 40., 0.),
                element,
            ));
        }
        for (a, b) in [(0, 1), (1, 2)] {
            molecule.bonds.push(crate::Bond {
                a,
                b,
                order: 1,
                style: crate::BondStyle::Normal,
                stereo: None,
                query: None,
            });
        }
        let molblock = crate::molfile_for(&molecule);

        let report = descriptors(&molblock).expect("RDKit descriptor worker should run");
        assert!(report.contains("RDKit "));
        assert!(report.contains("Molar mass:"));
        assert!(report.contains("logP:"));
        assert!(report.contains("Solvent-accessible surface area:"));
        if has_rdkit_freesasa(&python) {
            assert!(report.contains("Shrake-Rupley"));
            assert!(report.contains("seed 20261120"));
        } else {
            assert!(report.contains("unavailable (This RDKit build does not include rdFreeSASA)"));
        }
        assert_eq!(to_smiles(&molblock).unwrap(), "CCO");
        molecule.atoms[0].atom_map = Some(5);
        let mapped_smiles = to_smiles(&crate::molfile_for(&molecule)).unwrap();
        assert!(
            mapped_smiles.contains(":5]"),
            "unexpected mapped SMILES: {mapped_smiles}"
        );
        assert!(
            !mapped_smiles.contains(":0]"),
            "unexpected map zero: {mapped_smiles}"
        );
        let mut v3000_molecule = crate::Molecule::default();
        for index in 0..1000 {
            v3000_molecule.atoms.push(crate::Atom::new(
                macroquad::prelude::vec2(index as f32, 0.),
                "C",
            ));
        }
        let v3000 = crate::molfile_for(&v3000_molecule);
        assert!(v3000.contains("V3000"));
        let v3000_smiles = to_smiles(&v3000).unwrap();
        assert!(!v3000_smiles.contains(":0]"));
        v3000_molecule.atoms[0].atom_map = Some(7);
        let mapped_v3000 = to_smiles(&crate::molfile_for(&v3000_molecule)).unwrap();
        assert!(
            mapped_v3000.contains(":7]"),
            "unexpected V3000 SMILES: {mapped_v3000}"
        );
        assert!(from_smiles_with_stereo("CCO").unwrap().0.contains("M  END"));
        let generated = generate_3d(&molblock).unwrap();
        assert!(generated.contains("ENGINE\tRDKit "));
        assert!(generated.contains("METHOD\tMMFF94") || generated.contains("METHOD\tUFF"));
        let inchi = export_external("to-inchi", &molblock).unwrap();
        assert!(inchi.starts_with("InChI=1S/C2H6O"));
    }

    #[test]
    fn unavailable_engine_is_reported() {
        if std::env::var("MOLECULA_PYTHON").is_ok() {
            return;
        }
        let result = descriptors("invalid");
        assert!(result.is_err());
    }

    #[test]
    fn malformed_optional_format_conversions_never_succeed_silently() {
        assert!(import_external("from-inchi", "not an InChI").is_err());
        assert!(export_external("to-pdb", "not a MOL structure").is_err());
    }

    #[test]
    fn atommapper_pair_worker_returns_verified_mapping_records_when_installed() {
        let python = std::env::var("MOLECULA_PYTHON").unwrap_or_else(|_| "python3".to_string());
        let available = Command::new(python)
            .args(["-c", "import rdkit, atommapper"])
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|status| status.success());
        if !available {
            return;
        }
        let reactant = from_smiles_with_stereo("CCO").unwrap().0;
        let product = from_smiles_with_stereo("CC=O").unwrap().0;
        let report = map_reaction_pair(&reactant, &product).unwrap();
        assert!(report.contains("STATUS\tcomplete"));
        assert!(report.contains("ENGINE\tAtomMapper"));
        assert!(
            report
                .lines()
                .filter(|line| line.starts_with("PAIR\t"))
                .count()
                >= 3
        );
        assert!(report.lines().any(|line| line.ends_with("\t1")));
    }

    #[test]
    fn atommapper_worker_maps_disconnected_reaction_sides_when_installed() {
        let python = std::env::var("MOLECULA_PYTHON").unwrap_or_else(|_| "python3".to_string());
        let available = Command::new(python)
            .args(["-c", "import rdkit, atommapper"])
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|status| status.success());
        if !available {
            return;
        }
        let reactants = from_smiles_with_stereo("CC.O").unwrap().0;
        let products = from_smiles_with_stereo("CC.O").unwrap().0;
        let report = map_reaction_pair(&reactants, &products).unwrap();
        let pairs: Vec<_> = report
            .lines()
            .filter(|line| line.starts_with("PAIR\t"))
            .collect();
        assert!(report.contains("STATUS\tcomplete"));
        assert_eq!(pairs.len(), 3, "unexpected disconnected mapping: {report}");
        assert!(pairs.iter().any(|line| line.starts_with("PAIR\t2\t2\t")));
    }

    #[test]
    fn gasteiger_worker_reports_real_charges_with_model_and_units() {
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
        let block = from_smiles_with_stereo("CCO").unwrap().0;
        let report = partial_charge_report(&block).unwrap();
        assert!(report.contains("ENGINE\tRDKit "));
        assert!(report.contains("METHOD\tGasteiger-Marsili"));
        assert!(report.contains("UNIT\te (elementary charge)"));
        assert_eq!(
            report
                .lines()
                .filter(|line| line.starts_with("ATOM\t"))
                .count(),
            3
        );
        assert!(report.contains("\tC\t-0.041838\t0.076120"));
        assert!(report.contains("SUMMARY\tNet model charge\t0.000000 e\tFormal charge\t0"));
        assert!(partial_charge_report("invalid MOL block").is_err());
        let ammonium =
            partial_charge_report(&from_smiles_with_stereo("[NH4+]").unwrap().0).unwrap();
        assert!(ammonium.contains("SUMMARY\tNet model charge\t1.000000 e\tFormal charge\t1"));
    }

    #[test]
    fn rdkit_smiles_stereo_survives_mol_import_and_editor_graph_roundtrip() {
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
        for input in ["N[C@@H](C)C(=O)O", "N[C@H](C)C(=O)O", "F/C=C/F", "F/C=C\\F"] {
            let (rdkit_mol, stereo) = from_smiles_with_stereo(input).unwrap();
            let expected = to_smiles(&rdkit_mol).unwrap();
            let mut editor_graph =
                crate::parse_mol_at(&rdkit_mol, macroquad::math::Vec2::new(0., 0.)).unwrap();
            for (index, value) in stereo {
                editor_graph.bonds[index].stereo = Some(value);
            }
            let saved = crate::native::write(&crate::Snapshot {
                mol: editor_graph.clone(),
                graphics: vec![],
            });
            editor_graph = crate::native::read(&saved).unwrap().mol;
            let roundtrip = crate::smiles::export(&editor_graph)
                .or_else(|_| {
                    to_smiles_with_stereo(&crate::molfile_for(&editor_graph), &editor_graph)
                })
                .unwrap();
            assert_eq!(roundtrip, expected, "stereo lost while importing {input}");
        }
    }

    #[test]
    fn explicit_ez_annotation_overrides_conflicting_2d_layout_on_export() {
        let python = std::env::var("MOLECULA_PYTHON").unwrap_or_else(|_| "python3".into());
        if !Command::new(python)
            .args(["-c", "import rdkit"])
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .status()
            .is_ok_and(|status| status.success())
        {
            return;
        }
        let (trans_block, stereo) = from_smiles_with_stereo("F/C=C/F").unwrap();
        let (_, cis_block_stereo) = from_smiles_with_stereo("F/C=C\\F").unwrap();
        let mut trans =
            crate::parse_mol_at(&trans_block, macroquad::math::Vec2::new(0., 0.)).unwrap();
        let cis = crate::parse_mol_at(
            &from_smiles_with_stereo("F/C=C\\F").unwrap().0,
            macroquad::math::Vec2::new(0., 0.),
        )
        .unwrap();
        assert_eq!(stereo.len(), 1);
        assert_eq!(cis_block_stereo.len(), 1);
        for (atom, cis_atom) in trans.atoms.iter_mut().zip(cis.atoms.iter()) {
            atom.pos = cis_atom.pos;
        }
        for (index, value) in stereo {
            trans.bonds[index].stereo = Some(value);
        }
        let exported = to_smiles_with_stereo(&crate::molfile_for(&trans), &trans).unwrap();
        assert_eq!(exported, to_smiles(&trans_block).unwrap());
    }
}
