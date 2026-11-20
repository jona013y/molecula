//! A bounded SMILES reader/writer for ordinary chemical graphs.
//! Stereochemical and query SMILES are rejected rather than silently flattened.

use crate::{chemistry, smarts, Atom, Bond, BondStyle, Molecule};
use macroquad::prelude::vec2;
use std::collections::{BTreeSet, HashMap};

fn element_at(input: &str, cursor: &mut usize, aromatic: bool) -> Result<String, String> {
    let bytes = input.as_bytes();
    let start = *cursor;
    if *cursor >= bytes.len() {
        return Err("Missing SMILES atom symbol".into());
    }
    if aromatic {
        let element = match bytes[*cursor] as char {
            'b' => "B",
            'c' => "C",
            'n' => "N",
            'o' => "O",
            'p' => "P",
            's' => "S",
            other => return Err(format!("Unsupported aromatic SMILES atom `{other}`.")),
        };
        *cursor += 1;
        return Ok(element.into());
    }
    if !bytes[*cursor].is_ascii_uppercase() {
        return Err("Expected an element symbol in SMILES".into());
    }
    *cursor += 1;
    if *cursor < bytes.len() && bytes[*cursor].is_ascii_lowercase() {
        *cursor += 1;
    }
    let element = &input[start..*cursor];
    if !chemistry::valid_element(element) {
        return Err(format!("Unknown SMILES element `{element}`."));
    }
    Ok(element.into())
}

fn parse_bracket_atom(content: &str) -> Result<Atom, String> {
    let bytes = content.as_bytes();
    if bytes.is_empty() {
        return Err("Empty SMILES bracket atom.".into());
    }
    let mut cursor = 0usize;
    let isotope_start = cursor;
    while cursor < bytes.len() && bytes[cursor].is_ascii_digit() {
        cursor += 1;
    }
    let isotope = if cursor > isotope_start {
        let mass = content[isotope_start..cursor]
            .parse::<u16>()
            .map_err(|_| "SMILES isotope mass is outside the supported range")?;
        if mass == 0 {
            return Err("SMILES isotope mass must be greater than zero.".into());
        }
        Some(mass)
    } else {
        None
    };
    if cursor >= bytes.len() || bytes[cursor] == b'*' {
        return Err("Wildcard/query atoms are not supported in chemical SMILES import.".into());
    }
    let aromatic = bytes[cursor].is_ascii_lowercase();
    let element = element_at(content, &mut cursor, aromatic)?;
    let mut atom = Atom::new(vec2(0., 0.), &element);
    atom.aromatic = aromatic;
    atom.isotope = isotope;
    // In SMILES bracket atoms have no implicit hydrogens unless H is stated.
    atom.implicit_h_override = Some(0);

    if cursor < bytes.len() && bytes[cursor] == b'@' {
        return Err("SMILES atom chirality (@/@@) is not yet supported by local import.".into());
    }
    if cursor < bytes.len() && bytes[cursor] == b'H' {
        cursor += 1;
        let start = cursor;
        while cursor < bytes.len() && bytes[cursor].is_ascii_digit() {
            cursor += 1;
        }
        let count = if cursor == start {
            1
        } else {
            content[start..cursor]
                .parse::<u8>()
                .map_err(|_| "SMILES hydrogen count is outside the supported range")?
        };
        if count > 8 {
            return Err("SMILES bracket hydrogen counts are limited to H0–H8.".into());
        }
        atom.explicit_h = count;
    }
    if cursor < bytes.len() && matches!(bytes[cursor], b'+' | b'-') {
        let sign_byte = bytes[cursor];
        let sign = if sign_byte == b'+' { 1i16 } else { -1i16 };
        cursor += 1;
        let digit_start = cursor;
        while cursor < bytes.len() && bytes[cursor].is_ascii_digit() {
            cursor += 1;
        }
        let magnitude = if cursor > digit_start {
            content[digit_start..cursor]
                .parse::<i16>()
                .map_err(|_| "SMILES formal charge is outside the supported range")?
        } else {
            let mut repeated = 1i16;
            while cursor < bytes.len() && bytes[cursor] == sign_byte {
                repeated += 1;
                cursor += 1;
            }
            repeated
        };
        let charge = sign * magnitude;
        atom.charge = i8::try_from(charge)
            .map_err(|_| "SMILES formal charge is outside the supported range")?;
    }
    if cursor < bytes.len() && bytes[cursor] == b':' {
        cursor += 1;
        let start = cursor;
        while cursor < bytes.len() && bytes[cursor].is_ascii_digit() {
            cursor += 1;
        }
        if cursor == start {
            return Err("SMILES atom-map labels must be numeric.".into());
        }
        atom.atom_map = Some(
            content[start..cursor]
                .parse::<u16>()
                .map_err(|_| "SMILES atom-map label is outside the supported range")?,
        );
    }
    if cursor != bytes.len() {
        return Err(format!(
            "Unsupported SMILES bracket-atom syntax near `{}`.",
            &content[cursor..]
        ));
    }
    Ok(atom)
}

fn add_bond(molecule: &mut Molecule, a: usize, b: usize, order: Option<u8>) -> Result<(), String> {
    if a == b
        || molecule
            .bonds
            .iter()
            .any(|bond| (bond.a == a && bond.b == b) || (bond.a == b && bond.b == a))
    {
        return Err("SMILES contains a self-closure or duplicate bond.".into());
    }
    let order = order.unwrap_or_else(|| {
        if molecule.atoms[a].aromatic && molecule.atoms[b].aromatic {
            4
        } else {
            1
        }
    });
    if order == 4 && !(molecule.atoms[a].aromatic && molecule.atoms[b].aromatic) {
        return Err("An aromatic SMILES bond must join two aromatic atoms.".into());
    }
    molecule.bonds.push(Bond {
        a,
        b,
        order,
        style: if order == 4 {
            BondStyle::Aromatic
        } else {
            BondStyle::Normal
        },
        stereo: None,
        query: None,
    });
    Ok(())
}

/// Parse a documented subset of chemical SMILES: organic/aromatic atoms,
/// bracket isotope/charge/H/map annotations, branches, ring closures and dots.
/// Atom/bond stereochemistry and wildcard query atoms are deliberately refused.
pub fn parse(input: &str) -> Result<Molecule, String> {
    let input = input.trim();
    if input.is_empty() {
        return Err("The SMILES input is empty.".into());
    }
    if input.len() > 100_000 || !input.is_ascii() {
        return Err("SMILES input is limited to 100,000 ASCII bytes.".into());
    }
    let bytes = input.as_bytes();
    let mut molecule = Molecule::default();
    let mut current = None;
    let mut pending_bond = None;
    let mut branches = Vec::<(usize, usize)>::new();
    let mut rings = HashMap::<String, (usize, Option<u8>)>::new();
    let mut last_was_dot = false;
    let mut cursor = 0usize;
    while cursor < bytes.len() {
        let ch = bytes[cursor] as char;
        if ch.is_ascii_uppercase() || ch.is_ascii_lowercase() || ch == '[' {
            let atom = if ch == '[' {
                let close = input[cursor + 1..]
                    .find(']')
                    .map(|offset| cursor + 1 + offset)
                    .ok_or("Unclosed SMILES bracket atom")?;
                let atom = parse_bracket_atom(&input[cursor + 1..close])?;
                cursor = close + 1;
                atom
            } else {
                let aromatic = ch.is_ascii_lowercase();
                let element = element_at(input, &mut cursor, aromatic)?;
                if !aromatic
                    && !matches!(
                        element.as_str(),
                        "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
                    )
                {
                    return Err(format!(
                        "Element `{element}` must use bracket notation in SMILES."
                    ));
                }
                let mut atom = Atom::new(vec2(0., 0.), &element);
                atom.aromatic = aromatic;
                atom
            };
            let index = molecule.atoms.len();
            molecule.atoms.push(atom);
            if index >= 512 {
                return Err("Local SMILES import is limited to 512 atoms.".into());
            }
            if let Some(previous) = current {
                add_bond(&mut molecule, previous, index, pending_bond.take())?;
            } else if pending_bond.is_some() {
                return Err("A SMILES bond expression cannot begin a component.".into());
            }
            current = Some(index);
            last_was_dot = false;
            continue;
        }
        if matches!(ch, '-' | '=' | '#' | ':') {
            if current.is_none() || pending_bond.is_some() {
                return Err("Misplaced or repeated SMILES bond expression.".into());
            }
            pending_bond = Some(match ch {
                '-' => 1,
                '=' => 2,
                '#' => 3,
                ':' => 4,
                _ => unreachable!(),
            });
            cursor += 1;
            continue;
        }
        if ch.is_ascii_digit() || ch == '%' {
            let label = if ch == '%' {
                if cursor + 2 >= bytes.len()
                    || !bytes[cursor + 1..cursor + 3].iter().all(u8::is_ascii_digit)
                {
                    return Err("Local SMILES percent ring labels must contain two digits.".into());
                }
                let label = input[cursor + 1..cursor + 3].to_string();
                cursor += 3;
                label
            } else {
                cursor += 1;
                ch.to_string()
            };
            let atom = current.ok_or("SMILES ring label has no current atom")?;
            if let Some((other, opened_bond)) = rings.remove(&label) {
                if let (Some(first), Some(second)) = (opened_bond, pending_bond) {
                    if first != second {
                        return Err(
                            "The ends of a SMILES ring closure specify different bonds.".into()
                        );
                    }
                }
                add_bond(
                    &mut molecule,
                    other,
                    atom,
                    pending_bond.take().or(opened_bond),
                )?;
            } else {
                rings.insert(label, (atom, pending_bond.take()));
            }
            continue;
        }
        match ch {
            '(' => {
                if pending_bond.is_some() {
                    return Err("A SMILES branch must begin with an atom.".into());
                }
                let parent = current.ok_or("SMILES branch has no parent atom")?;
                branches.push((parent, molecule.atoms.len()));
                cursor += 1;
            }
            ')' => {
                if pending_bond.is_some() {
                    return Err("A SMILES branch ends with an incomplete bond.".into());
                }
                let (parent, start_len) = branches.pop().ok_or("Unmatched SMILES closing parenthesis")?;
                if molecule.atoms.len() == start_len {
                    return Err("An empty SMILES branch is invalid.".into());
                }
                current = Some(parent);
                cursor += 1;
            }
            '.' => {
                if pending_bond.is_some() || !branches.is_empty() || current.is_none() || last_was_dot {
                    return Err("A SMILES component separator cannot occur inside a branch or bond.".into());
                }
                current = None;
                last_was_dot = true;
                cursor += 1;
            }
            '/' | '\\' => return Err("Directional bond stereochemistry (/ and \\) is not yet supported by local SMILES import.".into()),
            _ => return Err(format!("Unsupported SMILES token `{ch}` at byte {}.", cursor + 1)),
        }
    }
    if !branches.is_empty() {
        return Err("Unclosed SMILES branch.".into());
    }
    if !rings.is_empty() {
        return Err("A SMILES ring label appears only once.".into());
    }
    if pending_bond.is_some() {
        return Err("SMILES input ends with an incomplete bond.".into());
    }
    if last_was_dot {
        return Err("SMILES input ends with an empty component.".into());
    }
    smarts::layout_imported_molecule(&mut molecule);
    Ok(molecule)
}

fn atom_token(atom: &Atom, inferred_implicit_h: Option<u8>) -> Result<String, String> {
    if atom.query.is_some()
        || atom.alias.is_some()
        || atom.r_group.is_some()
        || atom.stereo.is_some()
        || atom.radical != 0
        || atom.valence_override.is_some()
    {
        return Err("SMILES cannot preserve this atom's query, alias, R-group, stereochemical, radical, or editor-only data.".into());
    }
    let symbol = if atom.aromatic {
        match atom.element.as_str() {
            "B" => "b",
            "C" => "c",
            "N" => "n",
            "O" => "o",
            "P" => "p",
            "S" => "s",
            _ => {
                return Err(format!(
                    "Aromatic {} is not supported by local SMILES export.",
                    atom.element
                ))
            }
        }
    } else {
        atom.element.as_str()
    };
    let special = atom.isotope.is_some()
        || atom.charge != 0
        || atom.explicit_h != 0
        || atom.implicit_h_override.is_some()
        || atom.atom_map.is_some();
    if !special
        && !atom.aromatic
        && matches!(
            atom.element.as_str(),
            "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
        )
    {
        return Ok(symbol.into());
    }
    if !special && atom.aromatic {
        return Ok(symbol.into());
    }
    let mut token = String::from("[");
    if let Some(isotope) = atom.isotope {
        token.push_str(&isotope.to_string());
    }
    token.push_str(symbol);
    let implicit_h = atom
        .implicit_h_override
        .or(inferred_implicit_h)
        .ok_or_else(|| {
            format!(
                "Implicit hydrogen count for bracketed {} atom is unresolved.",
                atom.element
            )
        })?;
    let hydrogen_count = atom.explicit_h.saturating_add(implicit_h);
    if hydrogen_count > 0 {
        token.push('H');
        if hydrogen_count > 1 {
            token.push_str(&hydrogen_count.to_string());
        }
    }
    if atom.charge != 0 {
        token.push(if atom.charge > 0 { '+' } else { '-' });
        if atom.charge.unsigned_abs() > 1 {
            token.push_str(&atom.charge.unsigned_abs().to_string());
        }
    }
    if let Some(map) = atom.atom_map {
        token.push(':');
        token.push_str(&map.to_string());
    }
    token.push(']');
    Ok(token)
}

fn bond_token(bond: &Bond, first: &Atom, second: &Atom) -> Result<&'static str, String> {
    if bond.query.is_some() || bond.stereo.is_some() {
        return Err(
            "SMILES cannot preserve query bonds or bond stereochemistry in the local exporter."
                .into(),
        );
    }
    match (bond.style, bond.order) {
        (BondStyle::Normal, 1) => Ok(if first.aromatic && second.aromatic {
            "-"
        } else {
            ""
        }),
        (BondStyle::Normal, 2) => Ok("="),
        (BondStyle::Normal, 3) => Ok("#"),
        (BondStyle::Aromatic, 4) if first.aromatic && second.aromatic => Ok(""),
        (BondStyle::Aromatic, 4) => Ok(":"),
        _ => Err(
            "SMILES export supports only ordinary single, double, triple and aromatic bonds."
                .into(),
        ),
    }
}

fn ring_label(label: usize) -> String {
    if label < 10 {
        label.to_string()
    } else {
        format!("%{label:02}")
    }
}

fn render_tree(
    atom: usize,
    atoms: &[String],
    children: &[Vec<(usize, usize)>],
    annotations: &[Vec<(usize, String)>],
    molecule: &Molecule,
) -> Result<String, String> {
    let mut out = atoms[atom].clone();
    for (label, token) in &annotations[atom] {
        out.push_str(token);
        out.push_str(&ring_label(*label));
    }
    if let Some((&(continuation, edge), branches)) = children[atom].split_last() {
        for &(child, branch_edge) in branches {
            out.push('(');
            out.push_str(bond_token(
                &molecule.bonds[branch_edge],
                &molecule.atoms[atom],
                &molecule.atoms[child],
            )?);
            out.push_str(&render_tree(child, atoms, children, annotations, molecule)?);
            out.push(')');
        }
        out.push_str(bond_token(
            &molecule.bonds[edge],
            &molecule.atoms[atom],
            &molecule.atoms[continuation],
        )?);
        out.push_str(&render_tree(
            continuation,
            atoms,
            children,
            annotations,
            molecule,
        )?);
    }
    Ok(out)
}

/// Serialize ordinary molecules locally. Stereochemistry, query data, aliases,
/// radicals and nonstandard editor annotations are rejected without dropping them.
pub fn export(molecule: &Molecule) -> Result<String, String> {
    if molecule.atoms.is_empty() {
        return Err("Draw a structure before exporting SMILES.".into());
    }
    if molecule.atoms.len() > 512 {
        return Err("Local SMILES export is limited to 512 atoms.".into());
    }
    let analysis = chemistry::analyze(molecule);
    let atoms = molecule
        .atoms
        .iter()
        .enumerate()
        .map(|(index, atom)| atom_token(atom, analysis.implicit_h.get(index).copied().flatten()))
        .collect::<Result<Vec<_>, _>>()?;
    let mut adjacency = vec![Vec::<(usize, usize)>::new(); molecule.atoms.len()];
    let mut seen = BTreeSet::new();
    for (edge, bond) in molecule.bonds.iter().enumerate() {
        if bond.a >= molecule.atoms.len() || bond.b >= molecule.atoms.len() || bond.a == bond.b {
            return Err("The structure contains an invalid bond reference.".into());
        }
        let _ = bond_token(bond, &molecule.atoms[bond.a], &molecule.atoms[bond.b])?;
        if !seen.insert((bond.a.min(bond.b), bond.a.max(bond.b))) {
            return Err("SMILES export refuses duplicate bonds.".into());
        }
        adjacency[bond.a].push((bond.b, edge));
        adjacency[bond.b].push((bond.a, edge));
    }
    for neighbors in &mut adjacency {
        neighbors.sort_unstable_by_key(|(neighbor, _)| *neighbor);
    }
    let mut children = vec![Vec::new(); molecule.atoms.len()];
    let mut parent_edges = BTreeSet::new();
    let mut visited = vec![false; molecule.atoms.len()];
    let mut roots = Vec::new();
    for root in 0..molecule.atoms.len() {
        if visited[root] {
            continue;
        }
        roots.push(root);
        visited[root] = true;
        let mut stack = vec![(root, 0usize)];
        while let Some((current, next)) = stack.pop() {
            if next >= adjacency[current].len() {
                continue;
            }
            stack.push((current, next + 1));
            let (neighbor, edge) = adjacency[current][next];
            if !visited[neighbor] {
                visited[neighbor] = true;
                children[current].push((neighbor, edge));
                parent_edges.insert(edge);
                stack.push((neighbor, 0));
            }
        }
    }
    let ring_edges = (0..molecule.bonds.len())
        .filter(|edge| !parent_edges.contains(edge))
        .collect::<Vec<_>>();
    if ring_edges.len() > 99 {
        return Err("SMILES export is limited to 99 ring-closure labels.".into());
    }
    let mut annotations = vec![Vec::<(usize, String)>::new(); molecule.atoms.len()];
    for (index, &edge) in ring_edges.iter().enumerate() {
        let label = index + 1;
        let bond = &molecule.bonds[edge];
        let token = bond_token(bond, &molecule.atoms[bond.a], &molecule.atoms[bond.b])?.to_string();
        annotations[bond.a].push((label, token.clone()));
        annotations[bond.b].push((label, token));
    }
    for list in &mut annotations {
        list.sort_by_key(|(label, _)| *label);
    }
    roots
        .into_iter()
        .map(|root| render_tree(root, &atoms, &children, &annotations, molecule))
        .collect::<Result<Vec<_>, _>>()
        .map(|parts| parts.join("."))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn imports_linear_branching_aromatic_and_disconnected_smiles() {
        let ethanol = parse("CCO").unwrap();
        assert_eq!(ethanol.atoms.len(), 3);
        assert_eq!(
            ethanol
                .bonds
                .iter()
                .map(|bond| bond.order)
                .collect::<Vec<_>>(),
            [1, 1]
        );
        let branched = parse("CC(C)O").unwrap();
        assert_eq!(branched.bonds.len(), 3);
        let benzene = parse("c1ccccc1").unwrap();
        assert!(benzene.atoms.iter().all(|atom| atom.aromatic));
        assert!(benzene.bonds.iter().all(|bond| bond.order == 4));
        assert_eq!(parse("[Na+].[Cl-]").unwrap().atoms.len(), 2);
    }

    #[test]
    fn bracket_atoms_keep_isotope_charge_hydrogens_and_maps() {
        let molecule = parse("[13CH3+:7]O").unwrap();
        assert_eq!(molecule.atoms[0].isotope, Some(13));
        assert_eq!(molecule.atoms[0].charge, 1);
        assert_eq!(molecule.atoms[0].explicit_h, 3);
        assert_eq!(molecule.atoms[0].atom_map, Some(7));
        assert_eq!(molecule.atoms[1].element, "O");
        assert_eq!(chemistry::analyze(&parse("[C]").unwrap()).formula, "C");
        assert_eq!(chemistry::analyze(&parse("[CH4]").unwrap()).formula, "CH4");
    }

    #[test]
    fn rejects_unsupported_stereo_queries_and_malformed_smiles() {
        for text in ["C[C@H](O)F", "F/C=C/F", "C~N", "C(", "C1CC", "C..O", "[*]O"] {
            assert!(parse(text).is_err(), "{text}");
        }
    }

    #[test]
    fn smiles_roundtrips_supported_graph_features() {
        for text in [
            "CCO",
            "CC(C)O",
            "c1ccccc1",
            "[13CH3+:7]O",
            "[Na+].[Cl-]",
            "C1CCCCC1",
        ] {
            let parsed = parse(text).unwrap();
            let written = export(&parsed).unwrap();
            let restored = parse(&written).unwrap();
            assert_eq!(
                restored
                    .atoms
                    .iter()
                    .map(|atom| (
                        &atom.element,
                        atom.aromatic,
                        atom.isotope,
                        atom.charge,
                        atom.explicit_h,
                        atom.implicit_h_override,
                        atom.atom_map
                    ))
                    .collect::<Vec<_>>(),
                parsed
                    .atoms
                    .iter()
                    .map(|atom| (
                        &atom.element,
                        atom.aromatic,
                        atom.isotope,
                        atom.charge,
                        atom.explicit_h,
                        atom.implicit_h_override,
                        atom.atom_map
                    ))
                    .collect::<Vec<_>>(),
                "{text} => {written}"
            );
            assert_eq!(
                restored
                    .bonds
                    .iter()
                    .map(|bond| (bond.a.min(bond.b), bond.a.max(bond.b), bond.order))
                    .collect::<BTreeSet<_>>(),
                parsed
                    .bonds
                    .iter()
                    .map(|bond| (bond.a.min(bond.b), bond.a.max(bond.b), bond.order))
                    .collect::<BTreeSet<_>>(),
                "{text} => {written}"
            );
        }
        let mut mapped = parse("C").unwrap();
        mapped.atoms[0].atom_map = Some(9);
        let mapped_text = export(&mapped).unwrap();
        assert_eq!(mapped_text, "[CH4:9]");
        let restored = parse(&mapped_text).unwrap();
        assert_eq!(chemistry::analyze(&restored).formula, "CH4");

        let mut isotope = parse("C").unwrap();
        isotope.atoms[0].isotope = Some(13);
        let isotope_text = export(&isotope).unwrap();
        assert_eq!(isotope_text, "[13CH4]");
        let restored = parse(&isotope_text).unwrap();
        assert_eq!(chemistry::analyze(&restored).formula, "CH4");
    }
}
