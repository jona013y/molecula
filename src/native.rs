use crate::{Atom, Bond, BondStyle, DocumentTab, Graphic, Molecule, ShapeKind, Snapshot};
use macroquad::prelude::{vec2, Vec2};

fn hex(s: &str) -> String {
    s.as_bytes().iter().map(|b| format!("{b:02x}")).collect()
}
fn unhex(s: &str) -> Result<String, String> {
    if !s.len().is_multiple_of(2) {
        return Err("Invalid string encoding".into());
    }
    let mut bytes = Vec::new();
    if !s.is_ascii() {
        return Err("Invalid string encoding".into());
    }
    for i in (0..s.len()).step_by(2) {
        bytes.push(u8::from_str_radix(&s[i..i + 2], 16).map_err(|_| "Invalid string encoding")?);
    }
    String::from_utf8(bytes).map_err(|_| "Invalid UTF-8".into())
}
fn optstr(x: &Option<String>) -> String {
    x.as_ref().map(|s| hex(s)).unwrap_or_else(|| "-".into())
}
fn deopt(s: &str) -> Result<Option<String>, String> {
    if s == "-" {
        Ok(None)
    } else {
        Ok(Some(unhex(s)?))
    }
}
fn num<T: std::str::FromStr>(s: &str) -> Result<T, String> {
    s.parse().map_err(|_| format!("Invalid numeric field: {s}"))
}
fn optnum<T: std::str::FromStr>(s: &str) -> Result<Option<T>, String> {
    if s == "-" {
        Ok(None)
    } else {
        Ok(Some(num(s)?))
    }
}
fn on<T: std::fmt::Display>(x: Option<T>) -> String {
    x.map(|v| v.to_string()).unwrap_or_else(|| "-".into())
}
fn point(p: Vec2) -> String {
    format!("{}\t{}", p.x, p.y)
}

pub fn write(s: &Snapshot) -> String {
    let mut out = String::from("MOLECULAR-SKETCH\t4\n");
    out += &format!("NAME\t{}\n", hex(&s.mol.name));
    for (k, v) in &s.mol.properties {
        out += &format!("PROP\t{}\t{}\n", hex(k), hex(v));
    }
    for a in &s.mol.atoms {
        out += &format!(
            "ATOM\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            a.pos.x,
            a.pos.y,
            on(a.z),
            hex(&a.element),
            a.charge,
            on(a.isotope),
            a.radical,
            a.explicit_h,
            on(a.implicit_h_override),
            optstr(&a.alias),
            on(a.atom_map),
            a.aromatic as u8,
            optstr(&a.stereo),
            on(a.r_group),
            optstr(&a.query),
            optstr(&a.reaction_role),
            a.reaction_step,
            a.coefficient,
            a.coord3d.map(|p|p.x.to_string()).unwrap_or_else(||"-".into()),
            a.coord3d.map(|p|p.y.to_string()).unwrap_or_else(||"-".into()),
            a.coord3d.map(|p|p.z.to_string()).unwrap_or_else(||"-".into()),
            on(a.valence_override)
        );
    }
    for b in &s.mol.bonds {
        out += &format!(
            "BOND\t{}\t{}\t{}\t{}\t{}\t{}\n",
            b.a,
            b.b,
            b.order,
            b.style.code(),
            optstr(&b.stereo),
            optstr(&b.query)
        );
    }
    for g in &s.graphics {
        match g {
            Graphic::Arrow {
                start,
                end,
                kind,
                above,
                below,
            } => {
                out += &format!(
                    "ARROW\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    start.x,
                    start.y,
                    end.x,
                    end.y,
                    hex(kind),
                    hex(above),
                    hex(below)
                )
            }
            Graphic::Text {
                pos,
                text,
                size,
                bold,
                italic,
            } => {
                out += &format!(
                    "TEXT\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    pos.x,
                    pos.y,
                    size,
                    *bold as u8,
                    *italic as u8,
                    hex(text)
                )
            }
            Graphic::Plus { pos } => out += &format!("PLUS\t{}\n", point(*pos)),
            Graphic::ElectronArrow {
                start,
                control,
                end,
                electrons,
            } => {
                out += &format!(
                    "ELECTRON\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    start.x, start.y, control.x, control.y, end.x, end.y, electrons
                )
            }
            Graphic::Shape { kind, start, end } => {
                out += &format!(
                    "SHAPE\t{}\t{}\t{}\t{}\t{}\n",
                    kind.code(),
                    start.x,
                    start.y,
                    end.x,
                    end.y
                )
            }
        }
    }
    out += "END\n";
    out
}

pub fn write_workspace(tabs: &[DocumentTab], active: usize) -> String {
    let mut out = format!(
        "MOLECULA-WORKSPACE\t1\nACTIVE\t{}\n",
        active.min(tabs.len().saturating_sub(1))
    );
    for tab in tabs {
        out.push_str(&format!(
            "TAB\t{}\t{}\t{}\t{}\n",
            hex(&tab.title),
            tab.path.as_deref().map(hex).unwrap_or_else(|| "-".into()),
            tab.dirty as u8,
            hex(&write(&tab.snapshot))
        ));
    }
    out.push_str("END\n");
    out
}

pub fn read_workspace(input: &str) -> Result<(Vec<DocumentTab>, usize), String> {
    if input.len() > 128 * 1024 * 1024 {
        return Err("Workspace exceeds the 128 MiB safety limit".into());
    }
    let mut lines = input.lines();
    if lines.next() != Some("MOLECULA-WORKSPACE\t1") {
        return Err("Unsupported workspace version".into());
    }
    let active_line = lines.next().ok_or("Workspace is truncated")?;
    let active = active_line
        .strip_prefix("ACTIVE\t")
        .ok_or("Missing active tab")?
        .parse::<usize>()
        .map_err(|_| "Invalid active tab")?;
    let mut tabs = Vec::new();
    let mut ended = false;
    for line in lines {
        if line == "END" {
            ended = true;
            break;
        }
        let f: Vec<_> = line.split('\t').collect();
        if f.len() != 5 || f[0] != "TAB" {
            return Err("Malformed workspace tab".into());
        }
        let title = unhex(f[1])?;
        let path = deopt(f[2])?;
        let dirty = match f[3] {
            "0" => false,
            "1" => true,
            _ => return Err("Invalid dirty flag".into()),
        };
        let snapshot = read(&unhex(f[4])?)?;
        tabs.push(DocumentTab {
            title,
            path,
            snapshot,
            dirty,
            undo: Vec::new(),
            redo: Vec::new(),
        });
        if tabs.len() > 128 {
            return Err("Workspace exceeds 128 tabs".into());
        }
    }
    if !ended || tabs.is_empty() || active >= tabs.len() {
        return Err("Workspace is incomplete or has an invalid active tab".into());
    }
    Ok((tabs, active))
}

pub fn read(input: &str) -> Result<Snapshot, String> {
    // Reject unreasonable files before splitting lines or allocating decoded
    // fields. The record-count limits below still enforce semantic bounds.
    const MAX_DOCUMENT_BYTES: usize = 64 * 1024 * 1024;
    if input.len() > MAX_DOCUMENT_BYTES {
        return Err("Document exceeds the 64 MiB safety limit".into());
    }
    let mut lines = input.lines();
    let version = lines.next();
    if !matches!(
        version,
        Some(
            "MOLECULAR-SKETCH\t1"
                | "MOLECULAR-SKETCH\t2"
                | "MOLECULAR-SKETCH\t3"
                | "MOLECULAR-SKETCH\t4"
        )
    ) {
        return Err("Unsupported Molecular Sketch document version".into());
    }
    let mut mol = Molecule::default();
    let mut graphics = Vec::new();
    let mut ended = false;
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        match f.first().copied().unwrap_or("") {
            "NAME" if f.len() == 2 => mol.name = unhex(f[1])?,
            "PROP" if f.len() == 3 => mol.properties.push((unhex(f[1])?, unhex(f[2])?)),
            "ATOM" if f.len() == 16 || f.len() == 19 || f.len() == 22 || f.len() == 23 => {
                let mut a = Atom::new(vec2(num(f[1])?, num(f[2])?), &unhex(f[4])?);
                a.z = optnum(f[3])?;
                a.charge = num(f[5])?;
                a.isotope = optnum(f[6])?;
                a.radical = num(f[7])?;
                a.explicit_h = num(f[8])?;
                a.implicit_h_override = optnum(f[9])?;
                a.alias = deopt(f[10])?;
                a.atom_map = optnum(f[11])?;
                a.aromatic = f[12] == "1";
                a.stereo = deopt(f[13])?;
                a.r_group = optnum(f[14])?;
                a.query = deopt(f[15])?;
                if f.len() >= 19 {
                    a.reaction_role = deopt(f[16])?;
                    a.reaction_step = num(f[17])?;
                    a.coefficient = num(f[18])?;
                    if a.coefficient == 0 || a.coefficient > 99 || a.reaction_step == 0 {
                        return Err("Invalid reaction component metadata".into());
                    }
                    if a.reaction_role
                        .as_deref()
                        .is_some_and(|role| !matches!(role, "reactant" | "product" | "agent"))
                    {
                        return Err("Unknown reaction role".into());
                    }
                }
                if f.len() >= 22 {
                    if [f[19], f[20], f[21]].iter().all(|x| *x == "-") {
                        a.coord3d = None;
                    } else if [f[19], f[20], f[21]].contains(&"-") {
                        return Err("Partial 3D coordinate triplet".into());
                    } else {
                        let p = macroquad::prelude::vec3(num(f[19])?, num(f[20])?, num(f[21])?);
                        if !p.is_finite() {
                            return Err("Non-finite 3D coordinates".into());
                        }
                        a.coord3d = Some(p);
                    }
                }
                if f.len() == 23 {
                    a.valence_override = optnum(f[22])?;
                }
                mol.atoms.push(a);
            }
            "BOND" if f.len() == 7 => mol.bonds.push(Bond {
                a: num(f[1])?,
                b: num(f[2])?,
                order: num(f[3])?,
                style: BondStyle::parse(f[4]).ok_or("Unknown bond style")?,
                stereo: deopt(f[5])?,
                query: deopt(f[6])?,
            }),
            "ARROW" if f.len() == 8 => graphics.push(Graphic::Arrow {
                start: vec2(num(f[1])?, num(f[2])?),
                end: vec2(num(f[3])?, num(f[4])?),
                kind: unhex(f[5])?,
                above: unhex(f[6])?,
                below: unhex(f[7])?,
            }),
            "TEXT" if f.len() == 7 => graphics.push(Graphic::Text {
                pos: vec2(num(f[1])?, num(f[2])?),
                size: num(f[3])?,
                bold: f[4] == "1",
                italic: f[5] == "1",
                text: unhex(f[6])?,
            }),
            "PLUS" if f.len() == 3 => graphics.push(Graphic::Plus {
                pos: vec2(num(f[1])?, num(f[2])?),
            }),
            "ELECTRON" if f.len() == 8 => {
                let electrons: u8 = num(f[7])?;
                if electrons != 1 && electrons != 2 {
                    return Err("Electron arrow must carry one or two electrons".into());
                }
                graphics.push(Graphic::ElectronArrow {
                    start: vec2(num(f[1])?, num(f[2])?),
                    control: vec2(num(f[3])?, num(f[4])?),
                    end: vec2(num(f[5])?, num(f[6])?),
                    electrons,
                });
            }
            "SHAPE" if f.len() == 6 => graphics.push(Graphic::Shape {
                kind: ShapeKind::parse(f[1]).ok_or("Unknown graphic shape")?,
                start: vec2(num(f[2])?, num(f[3])?),
                end: vec2(num(f[4])?, num(f[5])?),
            }),
            "END" if f.len() == 1 => {
                ended = true;
                break;
            }
            _ => {
                return Err(format!(
                    "Malformed document line: {}",
                    f.first().unwrap_or(&"")
                ))
            }
        }
    }
    if !ended {
        return Err("Document is truncated".into());
    }
    if mol.atoms.len() > 100_000 || mol.bonds.len() > 200_000 {
        return Err("Document is too large".into());
    }
    for b in &mol.bonds {
        if b.a >= mol.atoms.len()
            || b.b >= mol.atoms.len()
            || b.a == b.b
            || !(1..=4).contains(&b.order)
        {
            return Err("Invalid bond references".into());
        }
    }
    if mol.atoms.iter().any(|a| {
        !a.pos.x.is_finite() || !a.pos.y.is_finite() || a.z.is_some_and(|z| !z.is_finite())
    }) || graphics.iter().any(|g| match g {
        Graphic::Arrow { start, end, .. } => !finite(*start) || !finite(*end),
        Graphic::Text { pos, size, .. } => !finite(*pos) || !size.is_finite() || *size <= 0.,
        Graphic::Plus { pos } => !finite(*pos),
        Graphic::ElectronArrow {
            start,
            control,
            end,
            ..
        } => !finite(*start) || !finite(*control) || !finite(*end),
        Graphic::Shape { start, end, .. } => !finite(*start) || !finite(*end),
    }) {
        return Err("Document contains non-finite coordinates or invalid text size".into());
    }
    Ok(Snapshot { mol, graphics })
}

fn finite(p: Vec2) -> bool {
    p.x.is_finite() && p.y.is_finite()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn roundtrip_preserves_chemistry_and_graphics() {
        let mut a = Atom::new(vec2(12., 35.), "N");
        a.isotope = Some(15);
        a.charge = 1;
        a.atom_map = Some(27);
        a.stereo = Some("R".into());
        a.reaction_role = Some("agent".into());
        a.reaction_step = 3;
        a.coefficient = 2;
        a.valence_override = Some(5);
        let mut s = Snapshot {
            mol: Molecule::default(),
            graphics: vec![
                Graphic::Arrow {
                    start: vec2(1., 2.),
                    end: vec2(3., 4.),
                    kind: "equilibrium".into(),
                    above: "heat".into(),
                    below: "".into(),
                },
                Graphic::ElectronArrow {
                    start: vec2(4., 5.),
                    control: vec2(7., 1.),
                    end: vec2(9., 6.),
                    electrons: 2,
                },
                Graphic::Shape {
                    kind: ShapeKind::Ellipse,
                    start: vec2(20., 30.),
                    end: vec2(80., 60.),
                },
            ],
        };
        s.mol.atoms.push(a);
        s.mol.name = "Test ø".into();
        assert_eq!(read(&write(&s)).unwrap().mol, s.mol);
        assert_eq!(read(&write(&s)).unwrap().graphics, s.graphics);
    }
    #[test]
    fn rejects_bad_bond_reference() {
        let s = "MOLECULAR-SKETCH\t1\nBOND\t0\t1\t1\tnormal\t-\t-\nEND\n";
        assert!(read(s).is_err());
    }

    #[test]
    fn reads_legacy_v1_documents_with_default_reaction_metadata() {
        let old = "MOLECULAR-SKETCH\t1\nATOM\t1\t2\t-\t43\t0\t-\t0\t0\t-\t-\t-\t0\t-\t-\t-\nEND\n";
        let loaded = read(old).unwrap();
        assert_eq!(loaded.mol.atoms[0].reaction_role, None);
        assert_eq!(loaded.mol.atoms[0].reaction_step, 1);
        assert_eq!(loaded.mol.atoms[0].coefficient, 1);
    }
    #[test]
    fn reads_legacy_v3_and_roundtrips_v4_valence_override() {
        let legacy = "MOLECULAR-SKETCH\t3\nATOM\t1\t2\t-\t43\t0\t-\t0\t0\t-\t-\t-\t0\t-\t-\t-\t-\t1\t1\t-\t-\t-\nEND\n";
        let loaded = read(legacy).unwrap();
        assert_eq!(loaded.mol.atoms[0].valence_override, None);
        let mut atom = Atom::new(vec2(0., 0.), "P");
        atom.valence_override = Some(5);
        let snapshot = Snapshot {
            mol: Molecule {
                atoms: vec![atom],
                ..Molecule::default()
            },
            graphics: vec![],
        };
        let encoded = write(&snapshot);
        assert!(encoded.starts_with("MOLECULAR-SKETCH\t4\n"));
        assert_eq!(read(&encoded).unwrap().mol, snapshot.mol);
    }
    #[test]
    fn rejects_invalid_electron_count() {
        let s = "MOLECULAR-SKETCH\t1\nELECTRON\t1\t2\t3\t4\t5\t6\t3\nEND\n";
        assert!(read(s).is_err());
    }

    #[test]
    fn rejects_non_finite_geometry_and_invalid_bond_orders() {
        let bad_position =
            "MOLECULAR-SKETCH\t1\nATOM\tNaN\t0\t-\t43\t0\t-\t0\t0\t-\t-\t-\t0\t-\t-\t-\nEND\n";
        assert!(read(bad_position).is_err());

        let bad_order = "MOLECULAR-SKETCH\t1\nATOM\t0\t0\t-\t43\t0\t-\t0\t0\t-\t-\t-\t0\t-\t-\t-\nATOM\t1\t0\t-\t43\t0\t-\t0\t0\t-\t-\t-\t0\t-\t-\t-\nBOND\t0\t1\t0\tnormal\t-\t-\nEND\n";
        assert!(read(bad_order).is_err());
    }

    #[test]
    fn rejects_oversized_documents_before_parsing() {
        let input = format!("MOLECULAR-SKETCH\t1\n{}", " ".repeat(64 * 1024 * 1024));
        assert!(matches!(read(&input), Err(error) if error.contains("64 MiB")));
    }

    #[test]
    fn workspace_roundtrip_preserves_tabs_active_index_and_dirty_state() {
        let snapshot = Snapshot {
            mol: Molecule::default(),
            graphics: vec![Graphic::Plus { pos: vec2(8., 9.) }],
        };
        let tabs = vec![
            DocumentTab {
                title: "One ø".into(),
                path: Some("/tmp/one.mols".into()),
                snapshot: snapshot.clone(),
                dirty: false,
                undo: vec![],
                redo: vec![],
            },
            DocumentTab {
                title: "Untitled".into(),
                path: None,
                snapshot,
                dirty: true,
                undo: vec![],
                redo: vec![],
            },
        ];
        let (loaded, active) = read_workspace(&write_workspace(&tabs, 1)).unwrap();
        assert_eq!(active, 1);
        assert_eq!(loaded.len(), 2);
        assert_eq!(loaded[0].title, "One ø");
        assert_eq!(loaded[0].path.as_deref(), Some("/tmp/one.mols"));
        assert!(!loaded[0].dirty);
        assert!(loaded[1].dirty);
        assert_eq!(loaded[1].snapshot.graphics.len(), 1);
    }
}
