//! Deliberately bounded systematic naming for simple, acyclic C/H/O/N/halogen graphs.
//! This is not a general IUPAC or preferred-name engine.
use crate::{BondStyle, Molecule};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NameResult {
    pub name: String,
    pub explanation: String,
    pub steps: Vec<String>,
    pub highlighted_atoms: Vec<usize>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Principal {
    CarboxylicAcid(usize),
    Nitrile(usize),
    Aldehyde(usize),
    Ketone(usize),
    Alcohol(usize),
}

fn root(n: usize) -> Option<&'static str> {
    Some(match n {
        1 => "meth",
        2 => "eth",
        3 => "prop",
        4 => "but",
        5 => "pent",
        6 => "hex",
        7 => "hept",
        8 => "oct",
        9 => "non",
        10 => "dec",
        11 => "undec",
        12 => "dodec",
        _ => return None,
    })
}

fn carbon_neighbors(m: &Molecule) -> Vec<Vec<usize>> {
    let mut out = vec![Vec::new(); m.atoms.len()];
    for b in &m.bonds {
        if m.atoms[b.a].element == "C" && m.atoms[b.b].element == "C" {
            out[b.a].push(b.b);
            out[b.b].push(b.a);
        }
    }
    out
}

fn has_cycle(neighbors: &[Vec<usize>], i: usize, parent: Option<usize>, seen: &mut [bool]) -> bool {
    seen[i] = true;
    for &j in &neighbors[i] {
        if !seen[j] {
            if has_cycle(neighbors, j, Some(i), seen) {
                return true;
            }
        } else if Some(j) != parent {
            return true;
        }
    }
    false
}

fn enumerate_paths(
    neighbors: &[Vec<usize>],
    current: usize,
    previous: Option<usize>,
    path: &mut Vec<usize>,
    all: &mut Vec<Vec<usize>>,
) {
    path.push(current);
    let mut extended = false;
    for &next in &neighbors[current] {
        if Some(next) != previous {
            extended = true;
            enumerate_paths(neighbors, next, Some(current), path, all);
        }
    }
    if !extended {
        all.push(path.clone());
    }
    path.pop();
}

fn bond_order(m: &Molecule, a: usize, b: usize) -> Option<u8> {
    m.bonds
        .iter()
        .find(|x| (x.a == a && x.b == b) || (x.a == b && x.b == a))
        .map(|x| x.order)
}

fn prefix_group(groups: &mut BTreeMap<String, Vec<usize>>, prefix: &str, locant: usize) {
    groups.entry(prefix.to_owned()).or_default().push(locant);
}

fn format_prefixes(groups: &mut BTreeMap<String, Vec<usize>>) -> String {
    let mut entries: Vec<(String, String)> = groups
        .iter_mut()
        .map(|(base, locants)| {
            locants.sort_unstable();
            let loc = locants
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",");
            let pref = if base.contains('-') {
                let mult = match locants.len() {
                    1 => "",
                    2 => "bis",
                    3 => "tris",
                    4 => "tetrakis",
                    _ => "",
                };
                if mult.is_empty() {
                    format!("{loc}-({base})")
                } else {
                    format!("{loc}-{mult}({base})")
                }
            } else {
                let mult = match locants.len() {
                    2 => "di",
                    3 => "tri",
                    4 => "tetra",
                    _ => "",
                };
                format!("{loc}-{mult}{base}")
            };
            (base.clone(), pref)
        })
        .collect();
    entries.sort_by(|a, b| a.0.cmp(&b.0));
    entries
        .into_iter()
        .map(|x| x.1)
        .collect::<Vec<_>>()
        .join("-")
}

/// Name a small saturated, carbon-only branch as a systematic alkyl substituent.
/// Only tree branches whose off-chain atoms are terminal methyl groups are accepted.
fn name_alkyl_substituent(m: &Molecule, branch: &[usize], attachment: usize) -> Option<String> {
    if branch.is_empty() || branch.len() > 6 || !branch.contains(&attachment) {
        return None;
    }
    let branch_set: BTreeSet<usize> = branch.iter().copied().collect();
    let mut neighbors = vec![Vec::new(); m.atoms.len()];
    for bond in &m.bonds {
        if branch_set.contains(&bond.a) && branch_set.contains(&bond.b) {
            if bond.order != 1 || bond.style != BondStyle::Normal {
                return None;
            }
            neighbors[bond.a].push(bond.b);
            neighbors[bond.b].push(bond.a);
        }
    }
    for &atom in branch {
        for bond in m
            .bonds
            .iter()
            .filter(|bond| bond.a == atom || bond.b == atom)
        {
            let other = if bond.a == atom { bond.b } else { bond.a };
            if branch_set.contains(&other) {
                continue;
            }
            if atom == attachment && m.atoms[other].element == "C" {
                if bond.order != 1 || bond.style != BondStyle::Normal {
                    return None;
                }
            } else if m.atoms[other].element != "H" || bond.order != 1 {
                return None;
            }
        }
    }
    let mut paths = Vec::new();
    for &atom in branch {
        enumerate_paths(&neighbors, atom, None, &mut Vec::new(), &mut paths);
    }
    paths.retain(|path| path.contains(&attachment));
    let max_len = paths.iter().map(Vec::len).max()?;
    paths.retain(|path| path.len() == max_len);

    let mut choices = Vec::new();
    for path in paths {
        for oriented in [path.clone(), path.into_iter().rev().collect()] {
            let Some(attachment_locant) = oriented.iter().position(|&atom| atom == attachment)
            else {
                continue;
            };
            let chain: BTreeSet<usize> = oriented.iter().copied().collect();
            let mut methyl_locants = Vec::new();
            let mut valid = true;
            for (position, &atom) in oriented.iter().enumerate() {
                for &neighbor in &neighbors[atom] {
                    if !branch_set.contains(&neighbor) || chain.contains(&neighbor) {
                        continue;
                    }
                    if neighbors[neighbor]
                        .iter()
                        .filter(|&&n| branch_set.contains(&n))
                        .count()
                        != 1
                    {
                        valid = false;
                        break;
                    }
                    methyl_locants.push(position + 1);
                }
                if !valid {
                    break;
                }
            }
            if !valid {
                continue;
            }
            methyl_locants.sort_unstable();
            let stem = root(oriented.len())?;
            let locants = methyl_locants
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",");
            let methyl = if methyl_locants.is_empty() {
                String::new()
            } else {
                let multiplier = match methyl_locants.len() {
                    1 => "",
                    2 => "di",
                    3 => "tri",
                    4 => "tetra",
                    _ => return None,
                };
                format!("{locants}-{multiplier}methyl")
            };
            let name = if oriented.len() == 1 {
                "methyl".to_string()
            } else if oriented.len() == 2 && methyl.is_empty() {
                "ethyl".to_string()
            } else {
                format!("{methyl}{stem}an-{}-yl", attachment_locant + 1)
            };
            choices.push((attachment_locant, methyl_locants, name));
        }
    }
    choices.sort_by(|a, b| a.0.cmp(&b.0).then_with(|| a.1.cmp(&b.1)));
    choices.into_iter().next().map(|choice| choice.2)
}

fn name_simple_carbocycle(m: &Molecule, carbons: &[usize]) -> Result<Option<NameResult>, String> {
    if carbons.len() < 3
        || m.atoms.iter().any(|atom| {
            !matches!(
                atom.element.as_str(),
                "C" | "H" | "O" | "F" | "Cl" | "Br" | "I"
            ) || atom.charge != 0
                || atom.isotope.is_some()
                || atom.aromatic
        })
        || m.bonds.iter().any(|bond| bond.style != BondStyle::Normal)
    {
        return Ok(None);
    }
    let cn = carbon_neighbors(m);
    let mut degree: Vec<usize> = cn.iter().map(Vec::len).collect();
    let mut removed = vec![false; m.atoms.len()];
    let mut leaves = std::collections::VecDeque::new();
    for &carbon in carbons {
        if degree[carbon] <= 1 {
            removed[carbon] = true;
            leaves.push_back(carbon);
        }
    }
    while let Some(atom) = leaves.pop_front() {
        for &neighbor in &cn[atom] {
            if !removed[neighbor] {
                degree[neighbor] -= 1;
                if degree[neighbor] <= 1 {
                    removed[neighbor] = true;
                    leaves.push_back(neighbor);
                }
            }
        }
    }
    let ring: BTreeSet<usize> = carbons
        .iter()
        .copied()
        .filter(|&atom| !removed[atom])
        .collect();
    if !(3..=12).contains(&ring.len()) || ring.iter().any(|&atom| degree[atom] != 2) {
        return Ok(None);
    }
    let ring_bonds: Vec<_> = m
        .bonds
        .iter()
        .filter(|bond| ring.contains(&bond.a) && ring.contains(&bond.b))
        .collect();
    if ring_bonds.len() != ring.len()
        || ring_bonds.iter().any(|bond| !matches!(bond.order, 1 | 2))
        || m.bonds.iter().any(|bond| {
            m.atoms[bond.a].element == "C"
                && m.atoms[bond.b].element == "C"
                && !(ring.contains(&bond.a) && ring.contains(&bond.b))
                && bond.order != 1
        })
    {
        return Ok(None);
    }
    let double_bonds = ring_bonds.iter().filter(|bond| bond.order == 2).count();
    if double_bonds > 3 {
        return Err("The built-in naming subset supports at most three ring double bonds.".into());
    }
    if ring.iter().any(|&atom| {
        ring_bonds
            .iter()
            .filter(|bond| bond.order == 2 && (bond.a == atom || bond.b == atom))
            .count()
            > 1
    }) {
        return Ok(None);
    }

    let mut side_groups = BTreeMap::<String, Vec<usize>>::new();
    let mut side_atoms = Vec::new();
    let mut seen_side_carbons = BTreeSet::new();
    for &ring_atom in &ring {
        for &root_atom in cn[ring_atom].iter().filter(|&&atom| !ring.contains(&atom)) {
            let mut branch = Vec::new();
            let mut previous = ring_atom;
            let mut current = root_atom;
            loop {
                if ring.contains(&current) || !seen_side_carbons.insert(current) {
                    return Ok(None);
                }
                let ring_neighbors = cn[current]
                    .iter()
                    .filter(|atom| ring.contains(atom))
                    .count();
                if ring_neighbors != usize::from(ring.contains(&previous)) {
                    return Ok(None);
                }
                branch.push(current);
                let next: Vec<usize> = cn[current]
                    .iter()
                    .copied()
                    .filter(|&atom| atom != previous && !ring.contains(&atom))
                    .collect();
                if next.len() > 1 {
                    return Ok(None);
                }
                let Some(&next_atom) = next.first() else {
                    break;
                };
                previous = current;
                current = next_atom;
            }
            let group = match branch.len() {
                1 => "methyl",
                2 => "ethyl",
                3 => "propyl",
                4 => "butyl",
                5 => "pentyl",
                6 => "hexyl",
                _ => return Ok(None),
            };
            if branch.len() > ring.len() {
                return Ok(None);
            }
            side_atoms.extend(branch);
            side_groups.entry(group.into()).or_default().push(ring_atom);
        }
    }
    if seen_side_carbons.len() + ring.len() != carbons.len() {
        return Ok(None);
    }

    for (index, atom) in m.atoms.iter().enumerate() {
        if matches!(atom.element.as_str(), "F" | "Cl" | "Br" | "I") {
            let bonds: Vec<_> = m
                .bonds
                .iter()
                .filter(|bond| bond.a == index || bond.b == index)
                .collect();
            if bonds.len() != 1 || bonds[0].order != 1 {
                return Ok(None);
            }
            let center = if bonds[0].a == index {
                bonds[0].b
            } else {
                bonds[0].a
            };
            if !ring.contains(&center) {
                return Ok(None);
            }
            let prefix = match atom.element.as_str() {
                "F" => "fluoro",
                "Cl" => "chloro",
                "Br" => "bromo",
                _ => "iodo",
            };
            side_groups.entry(prefix.into()).or_default().push(center);
            side_atoms.push(index);
        }
    }

    let mut hydroxyls = Vec::new();
    let mut carbonyls = Vec::new();
    for (oxygen, _) in m
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.element == "O")
    {
        let neighbors: Vec<_> = m
            .bonds
            .iter()
            .filter(|bond| bond.a == oxygen || bond.b == oxygen)
            .collect();
        if neighbors.len() != 1 {
            return Ok(None);
        }
        let bond = neighbors[0];
        let carbon = if bond.a == oxygen { bond.b } else { bond.a };
        if !ring.contains(&carbon) {
            return Ok(None);
        }
        match bond.order {
            1 => hydroxyls.push(carbon),
            2 => carbonyls.push(carbon),
            _ => return Ok(None),
        }
        side_atoms.push(oxygen);
    }
    if hydroxyls.len() > 1
        || carbonyls.len() > 1
        || (!carbonyls.is_empty() && !hydroxyls.is_empty())
    {
        return Ok(None);
    }
    if carbonyls.iter().any(|carbon| {
        ring_bonds
            .iter()
            .any(|bond| bond.order == 2 && (bond.a == *carbon || bond.b == *carbon))
    }) {
        return Ok(None);
    }
    for bond in &m.bonds {
        let first = m.atoms[bond.a].element.as_str();
        let second = m.atoms[bond.b].element.as_str();
        let valid = match (first, second) {
            ("C", "C") => true,
            ("C", "H") | ("H", "C") => bond.order == 1,
            ("C", "O") | ("O", "C") => bond.order == 1 || bond.order == 2,
            ("C", "F" | "Cl" | "Br" | "I") | ("F" | "Cl" | "Br" | "I", "C") => bond.order == 1,
            _ => false,
        };
        if !valid {
            return Ok(None);
        }
    }

    let mut numbering_candidates = Vec::new();
    for &start in &ring {
        for &first in cn[start].iter().filter(|neighbor| ring.contains(neighbor)) {
            let mut order = vec![start, first];
            while order.len() < ring.len() {
                let current = *order.last().unwrap();
                let previous = order[order.len() - 2];
                let next = *cn[current]
                    .iter()
                    .find(|&&atom| ring.contains(&atom) && atom != previous)
                    .ok_or("Invalid ring graph")?;
                if next == start {
                    break;
                }
                order.push(next);
            }
            if order.len() != ring.len() {
                continue;
            }
            let locants: BTreeMap<usize, usize> = order
                .iter()
                .enumerate()
                .map(|(index, &atom)| (atom, index + 1))
                .collect();
            let double_locants: Vec<usize> = (0..order.len())
                .filter(|&index| {
                    bond_order(m, order[index], order[(index + 1) % order.len()]) == Some(2)
                })
                .map(|index| index + 1)
                .collect();
            let mut substituent_locants: Vec<usize> = side_groups
                .values()
                .flat_map(|atoms| atoms.iter().filter_map(|atom| locants.get(atom).copied()))
                .collect();
            substituent_locants.sort_unstable();
            let principal_locant = carbonyls
                .first()
                .or_else(|| hydroxyls.first())
                .and_then(|atom| locants.get(atom).copied())
                .unwrap_or(0);
            numbering_candidates.push((
                principal_locant,
                double_locants,
                substituent_locants,
                locants,
            ));
        }
    }
    let (principal_locant, double_locants, _, locants) = numbering_candidates
        .into_iter()
        .min_by(|first, second| {
            first
                .0
                .cmp(&second.0)
                .then_with(|| first.1.cmp(&second.1))
                .then_with(|| first.2.cmp(&second.2))
        })
        .ok_or("Could not determine ring numbering")?;
    let stem =
        root(ring.len()).ok_or("This ring size is not supported by the built-in naming subset.")?;
    let substituent_count: usize = side_groups.values().map(Vec::len).sum();
    let mut prefix_groups = BTreeMap::<String, Vec<usize>>::new();
    for (group, atoms) in side_groups {
        for atom in atoms {
            prefix_groups
                .entry(group.clone())
                .or_default()
                .push(locants[&atom]);
        }
    }
    let mut prefix = format_prefixes(&mut prefix_groups);
    if substituent_count == 1 && double_bonds == 0 && carbonyls.is_empty() && hydroxyls.is_empty() {
        prefix = prefix.replacen("1-", "", 1);
    }
    let parent = if !carbonyls.is_empty() || !hydroxyls.is_empty() {
        let suffix = if !carbonyls.is_empty() { "one" } else { "ol" };
        if double_bonds > 0 {
            let locants = double_locants
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",");
            let ending = match double_bonds {
                1 => "en",
                2 => "dien",
                _ => "trien",
            };
            let connecting_a = if double_bonds > 1 { "a" } else { "" };
            format!("cyclo{stem}{connecting_a}-{locants}-{ending}-{principal_locant}-{suffix}")
        } else if !carbonyls.is_empty() {
            if substituent_count == 0 {
                format!("cyclo{stem}anone")
            } else {
                format!("cyclo{stem}an-{principal_locant}-one")
            }
        } else if substituent_count == 0 {
            format!("cyclo{stem}anol")
        } else {
            format!("cyclo{stem}an-{principal_locant}-ol")
        }
    } else if double_bonds == 1 {
        if substituent_count == 0 {
            format!("cyclo{stem}ene")
        } else {
            format!("cyclo{stem}-{}-ene", double_locants[0])
        }
    } else if double_bonds > 1 {
        let locants = double_locants
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(",");
        let ending = if double_bonds == 2 { "diene" } else { "triene" };
        format!("cyclo{stem}a-{locants}-{ending}")
    } else {
        format!("cyclo{stem}ane")
    };
    let name = format!("{prefix}{parent}");
    let principal_group = if !carbonyls.is_empty() {
        "cyclic ketone (-one)"
    } else if !hydroxyls.is_empty() {
        "alcohol (-ol)"
    } else {
        "no supported suffix group"
    };
    let steps = vec![
        format!("1. Principal group priority: {principal_group} is the supported principal group."),
        format!(
            "2. Parent: the {n}-carbon ring is the parent; attached straight alkyl chains are prefixes.",
            n = ring.len()
        ),
        if !carbonyls.is_empty() || !hydroxyls.is_empty() {
            "3. Numbering: the principal group is assigned ring position 1.".into()
        } else if double_bonds > 0 {
            format!("3. Numbering: the ring double bonds begin at locants {}.", double_locants.iter().map(usize::to_string).collect::<Vec<_>>().join(", "))
        } else {
            "3. Numbering: choose the direction that gives substituents the lowest locant set.".into()
        },
        format!("4. Unsaturation: {}.",
            if !carbonyls.is_empty() {
                "the ring contains one ketone carbonyl"
            } else if !hydroxyls.is_empty() {
                "the ring contains one hydroxyl group"
            } else if double_bonds > 0 {
                "the ring double bonds are included in the parent ending"
            } else {
                "all ring bonds are single"
            }
        ),
        if substituent_count == 0 {
            "5. Other groups: none.".into()
        } else {
            format!("5. Other groups: {prefix}; prefixes are alphabetized.")
        },
        format!("6. Assemble the supported name: {name}."),
    ];
    Ok(Some(NameResult {
        name,
        explanation: steps.join("\n")
            + "\n\nGenerated by Molecula's limited teaching rules; not a complete IUPAC/PIN name.",
        steps,
        highlighted_atoms: ring.into_iter().chain(side_atoms).collect(),
    }))
}

fn name_simple_ester(m: &Molecule, carbons: &[usize]) -> Result<Option<NameResult>, String> {
    let mut esters = Vec::new();
    for carbonyl in m.bonds.iter().filter(|bond| bond.order == 2) {
        let carbon = if m.atoms[carbonyl.a].element == "C" && m.atoms[carbonyl.b].element == "O" {
            carbonyl.a
        } else if m.atoms[carbonyl.b].element == "C" && m.atoms[carbonyl.a].element == "O" {
            carbonyl.b
        } else {
            continue;
        };
        for bond in m.bonds.iter().filter(|bond| bond.order == 1) {
            let oxygen = if bond.a == carbon && m.atoms[bond.b].element == "O" {
                bond.b
            } else if bond.b == carbon && m.atoms[bond.a].element == "O" {
                bond.a
            } else {
                continue;
            };
            if let Some(alkyl_carbon) = m.bonds.iter().find_map(|next| {
                if next.order != 1 {
                    return None;
                }
                let other = if next.a == oxygen {
                    Some(next.b)
                } else if next.b == oxygen {
                    Some(next.a)
                } else {
                    None
                }?;
                (other != carbon && m.atoms[other].element == "C").then_some(other)
            }) {
                esters.push((carbon, oxygen, alkyl_carbon));
            }
        }
    }
    if esters.is_empty() {
        return Ok(None);
    }
    if esters.len() != 1 {
        return Err("Only one ester group is supported by the built-in teaching namer.".into());
    }
    if m.atoms
        .iter()
        .any(|atom| atom.charge != 0 || atom.element == "N")
        || m.atoms.iter().filter(|atom| atom.element == "O").count() != 2
    {
        return Err("The ester namer currently supports neutral, unbranched C/H/O esters with one ester group.".into());
    }
    let (acyl_carbon, ester_oxygen, alkyl_carbon) = esters[0];
    let carbon_graph = carbon_neighbors(m);
    let component = |start: usize| {
        let mut found = BTreeSet::new();
        let mut pending = vec![start];
        while let Some(atom) = pending.pop() {
            if found.insert(atom) {
                pending.extend(carbon_graph[atom].iter().copied());
            }
        }
        found
    };
    let acyl = component(acyl_carbon);
    let alkyl = component(alkyl_carbon);
    if acyl.contains(&alkyl_carbon)
        || !acyl.is_disjoint(&alkyl)
        || acyl.union(&alkyl).count() != carbons.len()
    {
        return Err(
            "Branched, cyclic, or multiply connected ester carbon skeletons are not yet supported."
                .into(),
        );
    }
    let simple_path = |part: &BTreeSet<usize>, start: usize| {
        part.iter().all(|&atom| {
            carbon_graph[atom]
                .iter()
                .filter(|neighbor| part.contains(neighbor))
                .count()
                <= 2
        }) && {
            let mut visited = BTreeSet::new();
            let mut pending = vec![start];
            while let Some(atom) = pending.pop() {
                if visited.insert(atom) {
                    pending.extend(
                        carbon_graph[atom]
                            .iter()
                            .filter(|neighbor| part.contains(neighbor))
                            .copied(),
                    );
                }
            }
            visited == *part
                && part
                    .iter()
                    .filter(|&&atom| {
                        carbon_graph[atom]
                            .iter()
                            .filter(|neighbor| part.contains(neighbor))
                            .count()
                            == 1
                    })
                    .count()
                    == if part.len() == 1 { 0 } else { 2 }
        }
    };
    if !simple_path(&acyl, acyl_carbon) || !simple_path(&alkyl, alkyl_carbon) {
        return Err(
            "Branched ester carbon chains are not yet supported by the built-in teaching namer."
                .into(),
        );
    }
    if carbon_graph[acyl_carbon]
        .iter()
        .filter(|neighbor| acyl.contains(neighbor))
        .count()
        > 1
        || carbon_graph[alkyl_carbon]
            .iter()
            .filter(|neighbor| alkyl.contains(neighbor))
            .count()
            > 1
    {
        return Err(
            "The acyl carbonyl must be at the parent end and the alkoxy group must attach at a terminal carbon; internal alkyl attachment locants are not yet supported.".into(),
        );
    }
    if m.bonds.iter().any(|bond| {
        let first = m.atoms[bond.a].element.as_str();
        let second = m.atoms[bond.b].element.as_str();
        match (first, second, bond.order) {
            ("C", "C", 1) => {
                !(acyl.contains(&bond.a) && acyl.contains(&bond.b)
                    || alkyl.contains(&bond.a) && alkyl.contains(&bond.b))
            }
            ("C", "O", 2) | ("O", "C", 2) => !(bond.a == acyl_carbon || bond.b == acyl_carbon),
            ("C", "O", 1) | ("O", "C", 1) => {
                !((bond.a == acyl_carbon && bond.b == ester_oxygen)
                    || (bond.b == acyl_carbon && bond.a == ester_oxygen)
                    || (bond.a == ester_oxygen && bond.b == alkyl_carbon)
                    || (bond.b == ester_oxygen && bond.a == alkyl_carbon))
            }
            ("C", "H", 1) | ("H", "C", 1) => false,
            _ => true,
        }
    }) {
        return Err("Only saturated, straight-chain alkyl alkanoates are supported by the built-in teaching namer.".into());
    }
    let acyl_len = acyl.len();
    let alkyl_len = alkyl.len();
    let acyl_root = root(acyl_len)
        .ok_or("Ester parent chains longer than 12 carbons are not yet supported.")?;
    let alkyl_root = root(alkyl_len)
        .ok_or("Ester alkyl chains longer than 12 carbons are not yet supported.")?;
    let name = format!("{}yl {}anoate", alkyl_root, acyl_root);
    let steps = vec![
        "1. Functional group: identify the ester linkage, C(=O)–O–C.".into(),
        format!("2. Acid-derived parent: the carbonyl carbon belongs to the {acyl_len}-carbon acyl chain."),
        format!("3. Alkyl group: the oxygen is attached to a straight {alkyl_len}-carbon chain, named {alkyl_root}yl."),
        "4. Numbering: the carbonyl carbon is carbon 1 of the acid-derived parent.".into(),
        format!("5. Combine the alkyl group and alkanoate parent: {name}."),
    ];
    Ok(Some(NameResult {
        name: name.clone(),
        explanation: steps.join("\n")
            + "\n\nGenerated by Molecula's limited teaching rules; not a complete IUPAC/PIN name.",
        steps,
        highlighted_atoms: acyl.union(&alkyl).copied().collect(),
    }))
}

/// Generate a teachable name for an acyclic carbon parent with common carbonyl/alcohol groups.
pub fn name(m: &Molecule) -> Result<NameResult, String> {
    if m.atoms.is_empty() {
        return Err("Draw one molecule before naming it.".into());
    }
    if m.atoms.len() > 64 {
        return Err("The built-in naming subset is limited to 64 atoms.".into());
    }
    if m.atoms.iter().any(|a| {
        a.aromatic
            || a.r_group.is_some()
            || a.query.is_some()
            || a.radical > 0
            || a.isotope.is_some()
            || a.charge != 0
            || a.alias.is_some()
            || a.valence_override.is_some()
            || a.stereo.is_some()
    }) {
        return Err("This structure has aromatic, query, R-group, radical, isotope, charge, alias, valence-override, or stereochemical features outside the built-in naming subset.".into());
    }
    // Name exactly one connected component; a reaction canvas must be split into separate structures.
    let mut visited = vec![false; m.atoms.len()];
    let mut stack = vec![0usize];
    visited[0] = true;
    while let Some(i) = stack.pop() {
        for b in &m.bonds {
            let j = if b.a == i {
                Some(b.b)
            } else if b.b == i {
                Some(b.a)
            } else {
                None
            };
            if let Some(j) = j {
                if !visited[j] {
                    visited[j] = true;
                    stack.push(j);
                }
            }
        }
    }
    if visited.iter().any(|x| !*x) {
        return Err(
            "Name one connected molecule at a time; this canvas has multiple fragments.".into(),
        );
    }
    if m.atoms.iter().any(|a| {
        !matches!(
            a.element.as_str(),
            "C" | "H" | "O" | "N" | "F" | "Cl" | "Br" | "I"
        )
    }) {
        return Err(
            "The built-in naming subset supports C, H, O, N, F, Cl, Br, and I only.".into(),
        );
    }
    if m.bonds.iter().any(|b| {
        b.style != BondStyle::Normal
            || b.order == 4
            || b.order > 3
            || b.stereo.is_some()
            || b.query.is_some()
    }) {
        return Err("Special or aromatic bonds are outside the built-in naming subset.".into());
    }
    let cn = carbon_neighbors(m);
    let carbons: Vec<usize> = m
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(i, a)| (a.element == "C").then_some(i))
        .collect();
    if carbons.is_empty() {
        return Err("No carbon parent chain was found.".into());
    }
    if let Some(result) = name_simple_ester(m, &carbons)? {
        return Ok(result);
    }
    if let Some(result) = name_simple_carbocycle(m, &carbons)? {
        return Ok(result);
    }
    let mut seen = vec![false; m.atoms.len()];
    for &i in &carbons {
        if !seen[i] && has_cycle(&cn, i, None, &mut seen) {
            return Err("Ring systems are not supported by this naming subset.".into());
        }
    }
    let mut carbonyls = Vec::<(usize, usize, usize)>::new();
    let mut hydroxyls = Vec::new();
    let mut carboxylic_acids = Vec::new();
    let mut nitriles = Vec::new();
    for b in &m.bonds {
        if b.order == 3 {
            let center = if m.atoms[b.a].element == "C" && m.atoms[b.b].element == "N" {
                Some((b.a, b.b))
            } else if m.atoms[b.b].element == "C" && m.atoms[b.a].element == "N" {
                Some((b.b, b.a))
            } else {
                None
            };
            if let Some((carbon, nitrogen)) = center {
                nitriles.push((carbon, nitrogen));
            }
        }
        if b.order == 2 {
            let (c, o) = if m.atoms[b.a].element == "C" && m.atoms[b.b].element == "O" {
                (b.a, b.b)
            } else if m.atoms[b.b].element == "C" && m.atoms[b.a].element == "O" {
                (b.b, b.a)
            } else {
                continue;
            };
            let c_neighbors: Vec<usize> = cn[c].clone();
            carbonyls.push((c, o, c_neighbors.len()));
        }
        if b.order == 1 {
            let (c, o) = if m.atoms[b.a].element == "C" && m.atoms[b.b].element == "O" {
                (b.a, b.b)
            } else if m.atoms[b.b].element == "C" && m.atoms[b.a].element == "O" {
                (b.b, b.a)
            } else {
                continue;
            };
            let oxygen_degree = m.bonds.iter().filter(|x| x.a == o || x.b == o).count();
            if oxygen_degree == 1 && m.atoms[o].charge == 0 {
                hydroxyls.push((c, o));
            }
        }
    }
    for &(carbon, carbonyl_oxygen, _) in &carbonyls {
        if m.atoms[carbon].charge != 0
            || m.atoms[carbonyl_oxygen].charge != 0
            || m.bonds
                .iter()
                .filter(|bond| bond.a == carbonyl_oxygen || bond.b == carbonyl_oxygen)
                .count()
                != 1
        {
            return Err("Charged or nonstandard carbonyls are outside the naming subset.".into());
        }
        for bond in m
            .bonds
            .iter()
            .filter(|bond| bond.a == carbon || bond.b == carbon)
        {
            let other = if bond.a == carbon { bond.b } else { bond.a };
            if other == carbonyl_oxygen || m.atoms[other].element == "H" {
                continue;
            }
            let recognized_acid_oxygen = m.atoms[other].element == "O"
                && bond.order == 1
                && hydroxyls
                    .iter()
                    .any(|(center, oxygen)| *center == carbon && *oxygen == other);
            if m.atoms[other].element != "C" && !recognized_acid_oxygen {
                return Err("This carbonyl is attached to an unsupported heteroatom group.".into());
            }
        }
    }
    for &(c, _, _) in &carbonyls {
        let acid_oh = hydroxyls.iter().any(|(carbon, _)| *carbon == c);
        if acid_oh {
            carboxylic_acids.push(c);
        }
    }
    if nitriles.len() > 1 {
        return Err("Multiple nitrile groups are outside the built-in naming subset.".into());
    }
    if let Some(&(carbon, nitrogen)) = nitriles.first() {
        let nitrogen_degree = m
            .bonds
            .iter()
            .filter(|bond| bond.a == nitrogen || bond.b == nitrogen)
            .count();
        if nitrogen_degree != 1 || cn[carbon].len() > 1 {
            return Err("Only a terminal, neutral nitrile attached to one carbon parent chain is supported.".into());
        }
        if !carbonyls.is_empty() || !hydroxyls.is_empty() {
            return Err("Nitriles combined with aldehydes, ketones, alcohols or carboxylic acids are outside the built-in naming subset.".into());
        }
    }
    let aldehydes: Vec<usize> = carbonyls
        .iter()
        .filter_map(|(c, _, n)| (*n <= 1 && !carboxylic_acids.contains(c)).then_some(*c))
        .collect();
    if carboxylic_acids.len() > 1 {
        return Err("Multiple carboxylic acid groups are not yet supported.".into());
    }
    let ketones: Vec<usize> = carbonyls
        .iter()
        .filter_map(|(c, _, n)| (*n == 2).then_some(*c))
        .collect();
    if aldehydes.len() > 1 {
        return Err("Multiple aldehyde groups are not yet supported.".into());
    }
    let principal = if let Some(&c) = carboxylic_acids.first() {
        Principal::CarboxylicAcid(c)
    } else if let Some(&(c, _)) = nitriles.first() {
        Principal::Nitrile(c)
    } else if let Some(&c) = aldehydes.first() {
        Principal::Aldehyde(c)
    } else if let Some(&(c, _, _)) = carbonyls.iter().find(|(_, _, n)| *n == 2) {
        Principal::Ketone(c)
    } else if !hydroxyls.is_empty() {
        Principal::Alcohol(hydroxyls[0].0)
    } else {
        return Err("No supported aldehyde, ketone, or alcohol group was found.".into());
    };
    let principal_atoms: Vec<usize> = match principal {
        Principal::Ketone(_) => ketones.clone(),
        Principal::Alcohol(_) => hydroxyls.iter().map(|(carbon, _)| *carbon).collect(),
        Principal::CarboxylicAcid(carbon)
        | Principal::Nitrile(carbon)
        | Principal::Aldehyde(carbon) => vec![carbon],
    };
    if principal_atoms.len() > 4 {
        return Err(
            "More than four identical principal groups are outside the built-in naming subset."
                .into(),
        );
    }
    let mut paths = Vec::new();
    for &i in &carbons {
        enumerate_paths(&cn, i, None, &mut Vec::new(), &mut paths);
    }
    let mut candidates: Vec<Vec<usize>> = paths
        .into_iter()
        .filter(|p| principal_atoms.iter().all(|center| p.contains(center)))
        .collect();
    let max_len = candidates.iter().map(Vec::len).max().unwrap_or(0);
    candidates.retain(|p| p.len() == max_len);
    let mut oriented = Vec::new();
    for mut p in candidates {
        let orient = |p: &[usize]| -> (Vec<usize>, Vec<usize>, Vec<usize>) {
            let mut group_locs: Vec<usize> = principal_atoms
                .iter()
                .filter_map(|center| p.iter().position(|x| x == center).map(|i| i + 1))
                .collect();
            group_locs.sort_unstable();
            let unsats: Vec<usize> = (0..p.len().saturating_sub(1))
                .filter_map(|i| {
                    bond_order(m, p[i], p[i + 1])
                        .filter(|o| *o > 1)
                        .map(|_| i + 1)
                })
                .collect();
            let parent: BTreeSet<usize> = p.iter().copied().collect();
            let mut prefix_locants = Vec::new();
            for (locant, &carbon) in p.iter().enumerate() {
                for bond in &m.bonds {
                    let other = if bond.a == carbon {
                        Some(bond.b)
                    } else if bond.b == carbon {
                        Some(bond.a)
                    } else {
                        None
                    };
                    if let Some(other) = other {
                        if !parent.contains(&other)
                            && m.atoms[other].element != "H"
                            && !(matches!(principal, Principal::Alcohol(_))
                                && hydroxyls
                                    .iter()
                                    .any(|(center, oxygen)| *center == carbon && *oxygen == other))
                            && !matches!(principal, Principal::CarboxylicAcid(c) if c == carbon && m.atoms[other].element == "O")
                            && !matches!(principal, Principal::Nitrile(c) if c == carbon && m.atoms[other].element == "N")
                            && !matches!(principal, Principal::Aldehyde(c) if c == carbon && m.atoms[other].element == "O")
                        {
                            prefix_locants.push(locant + 1);
                        }
                    }
                }
            }
            prefix_locants.sort_unstable();
            (group_locs, unsats, prefix_locants)
        };
        if matches!(principal,Principal::CarboxylicAcid(c) | Principal::Nitrile(c) | Principal::Aldehyde(c) if p.last()==Some(&c))
        {
            p.reverse();
        } else if !matches!(
            principal,
            Principal::CarboxylicAcid(_) | Principal::Nitrile(_) | Principal::Aldehyde(_)
        ) {
            let rev: Vec<usize> = p.iter().rev().copied().collect();
            if orient(&rev) < orient(&p) {
                p = rev;
            }
        }
        let (loc, unsat, prefixes) = orient(&p);
        oriented.push((loc, unsat, prefixes, p));
    }
    oriented.sort_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| a.1.cmp(&b.1))
            .then_with(|| a.2.cmp(&b.2))
    });
    let chain = oriented
        .into_iter()
        .next()
        .ok_or("No parent chain includes the principal group")?
        .3;
    let n = chain.len();
    let stem = root(n).ok_or("Parent chains longer than 12 carbons are not yet supported.")?;
    let locants: BTreeMap<usize, usize> =
        chain.iter().enumerate().map(|(i, &a)| (a, i + 1)).collect();
    let mut groups = BTreeMap::<String, Vec<usize>>::new();
    for &(c, _, neighbors) in &carbonyls {
        if neighbors == 2 {
            if let Some(&l) = locants.get(&c) {
                if !(matches!(principal, Principal::Ketone(_)) && ketones.contains(&c)) {
                    prefix_group(&mut groups, "oxo", l);
                }
            }
        }
    }
    for &(c, _) in &hydroxyls {
        if matches!(principal, Principal::CarboxylicAcid(acid) if acid == c) {
            continue;
        }
        if let Some(&l) = locants.get(&c) {
            if !matches!(principal, Principal::Alcohol(_)) {
                prefix_group(&mut groups, "hydroxy", l);
            }
        }
    }
    let parent_set: BTreeSet<usize> = chain.iter().copied().collect();
    let mut covered = BTreeSet::new();
    let mut named_amino_atoms = BTreeSet::new();
    for (i, &c) in chain.iter().enumerate() {
        for &neighbor in &cn[c] {
            if parent_set.contains(&neighbor) || covered.contains(&neighbor) {
                continue;
            }
            let mut branch = Vec::new();
            let mut todo = vec![(neighbor, c)];
            while let Some((x, prev)) = todo.pop() {
                if parent_set.contains(&x) || covered.contains(&x) {
                    continue;
                }
                covered.insert(x);
                branch.push(x);
                if cn[x].iter().any(|j| {
                    !parent_set.contains(j) && *j != prev && bond_order(m, x, *j) != Some(1)
                }) {
                    return Err("Unsaturated or branched side chains are not yet supported.".into());
                }
                for &j in &cn[x] {
                    if j != prev && !parent_set.contains(&j) {
                        todo.push((j, x));
                    }
                }
            }
            let alkyl = name_alkyl_substituent(m, &branch, neighbor).ok_or_else(|| {
                "This branched or functionalized alkyl substituent is outside the built-in naming subset.".to_string()
            })?;
            prefix_group(&mut groups, &alkyl, i + 1);
        }
    }
    for (i, &c) in chain.iter().enumerate() {
        for b in &m.bonds {
            let x = if b.a == c {
                b.b
            } else if b.b == c {
                b.a
            } else {
                continue;
            };
            if m.atoms[x].element == "N" && b.order == 1 {
                let nitrogen_bonds: Vec<_> = m
                    .bonds
                    .iter()
                    .filter(|bond| bond.a == x || bond.b == x)
                    .collect();
                let valid_primary_amine = m.atoms[x].charge == 0
                    && nitrogen_bonds.iter().all(|bond| {
                        if bond.order != 1 {
                            return false;
                        }
                        let other = if bond.a == x { bond.b } else { bond.a };
                        other == c || m.atoms[other].element == "H"
                    });
                if !valid_primary_amine {
                    return Err(
                        "Only neutral primary amine substituents are supported by the built-in naming subset.".into(),
                    );
                }
                named_amino_atoms.insert(x);
                prefix_group(&mut groups, "amino", i + 1);
            } else if matches!(m.atoms[x].element.as_str(), "F" | "Cl" | "Br" | "I") && b.order == 1
            {
                let p = match m.atoms[x].element.as_str() {
                    "F" => "fluoro",
                    "Cl" => "chloro",
                    "Br" => "bromo",
                    _ => "iodo",
                };
                prefix_group(&mut groups, p, i + 1);
            }
        }
    }
    for (index, atom) in m.atoms.iter().enumerate() {
        match atom.element.as_str() {
            "H" => {}
            "C" if parent_set.contains(&index) || covered.contains(&index) => {}
            "N" if named_amino_atoms.contains(&index) => {}
            "N" if nitriles.iter().any(|(_, nitrogen)| *nitrogen == index) => {}
            "O" if carbonyls.iter().any(|(_, oxygen, _)| *oxygen == index)
                || hydroxyls.iter().any(|(_, oxygen)| *oxygen == index) => {}
            "F" | "Cl" | "Br" | "I"
                if m.bonds
                    .iter()
                    .filter(|bond| bond.a == index || bond.b == index)
                    .count()
                    == 1
                    && m.bonds.iter().any(|bond| {
                        let other = if bond.a == index {
                            Some(bond.b)
                        } else if bond.b == index {
                            Some(bond.a)
                        } else {
                            None
                        };
                        bond.order == 1 && other.is_some_and(|other| parent_set.contains(&other))
                    }) => {}
            _ => {
                return Err(format!(
                    "Atom {} ({}) is not part of a supported functional group or substituent.",
                    index + 1,
                    atom.element
                ));
            }
        }
    }
    let mut doubles = Vec::new();
    let mut triples = Vec::new();
    for i in 0..n - 1 {
        match bond_order(m, chain[i], chain[i + 1]) {
            Some(2) => doubles.push(i + 1),
            Some(3) => triples.push(i + 1),
            _ => {}
        }
    }
    let unsaturated = !doubles.is_empty() || !triples.is_empty();
    let unsat = if !unsaturated {
        format!("{stem}an")
    } else {
        let mut s = stem.to_owned();
        if !doubles.is_empty() {
            s += &format!(
                "-{}-{}",
                doubles
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(","),
                if doubles.len() > 1 { "dien" } else { "en" }
            );
        }
        if !triples.is_empty() {
            s += &format!(
                "-{}-{}",
                triples
                    .iter()
                    .map(usize::to_string)
                    .collect::<Vec<_>>()
                    .join(","),
                if triples.len() > 1 { "diyn" } else { "yn" }
            );
        }
        s
    };
    let base = match principal {
        Principal::CarboxylicAcid(_) => format!("{unsat}oic acid"),
        Principal::Nitrile(_) => format!("{unsat}enitrile"),
        Principal::Aldehyde(_) => {
            if unsaturated {
                format!("{unsat}al")
            } else {
                format!("{stem}anal")
            }
        }
        Principal::Ketone(_) if ketones.len() > 1 => {
            let suffix = match ketones.len() {
                2 => "dione",
                3 => "trione",
                4 => "tetrone",
                _ => unreachable!("larger principal-group counts are rejected above"),
            };
            let mut ketone_locants: Vec<usize> = ketones
                .iter()
                .filter_map(|center| locants.get(center).copied())
                .collect();
            ketone_locants.sort_unstable();
            let locant_list = ketone_locants
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",");
            if unsaturated {
                format!("{unsat}-{locant_list}-{suffix}")
            } else {
                format!("{stem}ane-{locant_list}-{suffix}")
            }
        }
        Principal::Ketone(c) => format!("{unsat}-{}-one", locants[&c]),
        Principal::Alcohol(_) if hydroxyls.len() > 1 => {
            let suffix = match hydroxyls.len() {
                2 => "diol",
                3 => "triol",
                4 => "tetrol",
                _ => unreachable!("larger principal-group counts are rejected above"),
            };
            let mut alcohol_locants: Vec<usize> = hydroxyls
                .iter()
                .filter_map(|(center, _)| locants.get(center).copied())
                .collect();
            alcohol_locants.sort_unstable();
            let locant_list = alcohol_locants
                .iter()
                .map(usize::to_string)
                .collect::<Vec<_>>()
                .join(",");
            if unsaturated {
                format!("{unsat}-{locant_list}-{suffix}")
            } else {
                format!("{stem}ane-{locant_list}-{suffix}")
            }
        }
        Principal::Alcohol(c) => format!("{unsat}-{}-ol", locants[&c]),
    };
    let prefix = format_prefixes(&mut groups);
    let name = if prefix.is_empty() {
        base
    } else {
        format!("{prefix}{base}")
    };
    let (group_text, principal_atom) = match principal {
        Principal::CarboxylicAcid(c) => ("carboxylic acid (-oic acid)", c),
        Principal::Nitrile(c) => ("nitrile (-nitrile)", c),
        Principal::Aldehyde(c) => ("aldehyde (-al)", c),
        Principal::Ketone(c) => ("ketone (-one)", c),
        Principal::Alcohol(c) => ("alcohol (-ol)", c),
    };
    let unsat_text = if doubles.is_empty() && triples.is_empty() {
        "no multiple bonds".into()
    } else {
        format!(
            "multiple-bond locants: {}",
            doubles
                .iter()
                .map(|x| format!("{x}= "))
                .chain(triples.iter().map(|x| format!("{x}≡ ")))
                .collect::<String>()
        )
    };
    let parent = format!(
        "{} carbons: {}",
        n,
        if unsaturated {
            "alkene/alkyne"
        } else {
            "alkane"
        }
    );
    let prefix_text = if prefix.is_empty() {
        "none".to_string()
    } else {
        prefix.clone()
    };
    let mut principal_locants: Vec<usize> = principal_atoms
        .iter()
        .filter_map(|center| locants.get(center).copied())
        .collect();
    principal_locants.sort_unstable();
    let principal_locant_text = principal_locants
        .iter()
        .map(usize::to_string)
        .collect::<Vec<_>>()
        .join(",");
    let steps=vec![
        format!("1. Principal group priority: {group_text} is the highest supported group found."),
        format!("2. Parent chain: choose the longest supported carbon chain containing the principal atom ({parent})."),
        format!("3. Numbering: principal-group locant(s) are {principal_locant_text}; then compare multiple-bond and substituent locants."),
        format!("4. Multiple bonds: {unsat_text}."),
        format!("5. Other groups: {prefix_text}; prefixes are alphabetized."),
        format!("6. Assemble the supported name: {name}."),
    ];
    let explanation = steps.join("\n")
        + "\n\nGenerated by Molecula's limited teaching rules; not a complete IUPAC/PIN name.";
    let mut highlighted_atoms = chain;
    highlighted_atoms.push(principal_atom);
    highlighted_atoms.sort_unstable();
    highlighted_atoms.dedup();
    Ok(NameResult {
        name,
        explanation,
        steps,
        highlighted_atoms,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Atom, Bond};
    use macroquad::prelude::vec2;
    fn bond(m: &mut Molecule, a: usize, b: usize, order: u8) {
        m.bonds.push(Bond {
            a,
            b,
            order,
            style: BondStyle::Normal,
            stereo: None,
            query: None,
        });
    }
    #[test]
    fn names_simple_cycloalkanes_and_cycloalkenes() {
        for (size, double_bond, expected) in
            [(6, None, "cyclohexane"), (5, Some((0, 1)), "cyclopentene")]
        {
            let mut m = Molecule::default();
            for _ in 0..size {
                m.atoms.push(Atom::new(vec2(0., 0.), "C"));
            }
            for i in 0..size {
                let order = if double_bond == Some((i, (i + 1) % size))
                    || double_bond == Some(((i + 1) % size, i))
                {
                    2
                } else {
                    1
                };
                bond(&mut m, i, (i + 1) % size, order);
            }
            let result = name(&m).unwrap();
            assert_eq!(result.name, expected);
            assert_eq!(result.highlighted_atoms.len(), size);
            assert!(result.explanation.contains("ring"));
        }
    }
    #[test]
    fn names_substituted_cycloalkanes_and_cycloalkenes_with_lowest_locants() {
        let mut substituted = Molecule::default();
        for _ in 0..6 {
            substituted.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        substituted.atoms.push(Atom::new(vec2(0., 0.), "C"));
        substituted.atoms.push(Atom::new(vec2(0., 0.), "Cl"));
        for index in 0..6 {
            bond(&mut substituted, index, (index + 1) % 6, 1);
        }
        bond(&mut substituted, 0, 6, 1);
        bond(&mut substituted, 2, 7, 1);
        let result = name(&substituted).unwrap();
        assert_eq!(result.name, "3-chloro-1-methylcyclohexane");
        assert!(result.highlighted_atoms.contains(&6));
        assert!(result.highlighted_atoms.contains(&7));

        let mut alkene = Molecule::default();
        for _ in 0..7 {
            alkene.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        for index in 0..6 {
            bond(
                &mut alkene,
                index,
                (index + 1) % 6,
                if index == 0 { 2 } else { 1 },
            );
        }
        bond(&mut alkene, 0, 6, 1);
        let result = name(&alkene).unwrap();
        assert_eq!(result.name, "1-methylcyclohex-1-ene");
    }
    #[test]
    fn names_simple_cyclic_ketone_and_alcohol() {
        for (oxygen_order, expected) in [(2, "cyclohexanone"), (1, "cyclohexanol")] {
            let mut m = Molecule::default();
            for _ in 0..6 {
                m.atoms.push(Atom::new(vec2(0., 0.), "C"));
            }
            m.atoms.push(Atom::new(vec2(0., 0.), "O"));
            for i in 0..6 {
                bond(&mut m, i, (i + 1) % 6, 1);
            }
            bond(&mut m, 0, 6, oxygen_order);
            let result = name(&m).unwrap();
            assert_eq!(result.name, expected);
            assert!(result.steps[0].contains(if oxygen_order == 2 {
                "ketone"
            } else {
                "alcohol"
            }));
        }
    }
    #[test]
    fn names_supported_diols_and_diketones_as_multiple_principal_suffixes() {
        let mut diol = Molecule::default();
        for _ in 0..6 {
            diol.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        for _ in 0..2 {
            diol.atoms.push(Atom::new(vec2(0., 0.), "O"));
        }
        for i in 0..5 {
            bond(&mut diol, i, i + 1, 1);
        }
        bond(&mut diol, 0, 6, 1);
        bond(&mut diol, 1, 7, 1);
        assert_eq!(name(&diol).unwrap().name, "hexane-1,2-diol");

        let mut diketone = Molecule::default();
        for _ in 0..6 {
            diketone.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        for _ in 0..2 {
            diketone.atoms.push(Atom::new(vec2(0., 0.), "O"));
        }
        for i in 0..5 {
            bond(&mut diketone, i, i + 1, 1);
        }
        bond(&mut diketone, 1, 6, 2);
        bond(&mut diketone, 3, 7, 2);
        assert_eq!(name(&diketone).unwrap().name, "hexane-2,4-dione");
    }
    #[test]
    fn names_substituted_cyclic_ketone_and_alcohol_with_principal_group_numbering() {
        for (oxygen_order, expected) in [
            (2, "2-methylcyclohexan-1-one"),
            (1, "2-methylcyclohexan-1-ol"),
        ] {
            let mut m = Molecule::default();
            for _ in 0..6 {
                m.atoms.push(Atom::new(vec2(0., 0.), "C"));
            }
            m.atoms.push(Atom::new(vec2(0., 0.), "O"));
            m.atoms.push(Atom::new(vec2(0., 0.), "C"));
            for i in 0..6 {
                bond(&mut m, i, (i + 1) % 6, 1);
            }
            bond(&mut m, 0, 6, oxygen_order);
            bond(&mut m, 1, 7, 1);
            assert_eq!(name(&m).unwrap().name, expected);
        }
    }
    #[test]
    fn names_cyclic_enols_and_enones_with_principal_group_numbering() {
        for (carbonyl_order, expected) in [(1, "cyclohex-2-en-1-ol"), (2, "cyclohex-2-en-1-one")] {
            let mut m = Molecule::default();
            for _ in 0..6 {
                m.atoms.push(Atom::new(vec2(0., 0.), "C"));
            }
            m.atoms.push(Atom::new(vec2(0., 0.), "O"));
            for i in 0..6 {
                bond(&mut m, i, (i + 1) % 6, if i == 0 { 2 } else { 1 });
            }
            bond(&mut m, 2, 6, carbonyl_order);
            let result = name(&m).unwrap();
            assert_eq!(result.name, expected);
            assert!(result.steps[0].contains(if carbonyl_order == 1 {
                "alcohol"
            } else {
                "ketone"
            }));
        }
    }
    #[test]
    fn names_cyclic_dienes_with_and_without_a_principal_suffix() {
        let mut diene = Molecule::default();
        for _ in 0..6 {
            diene.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        for i in 0..6 {
            bond(
                &mut diene,
                i,
                (i + 1) % 6,
                if i == 0 || i == 2 { 2 } else { 1 },
            );
        }
        assert_eq!(name(&diene).unwrap().name, "cyclohexa-1,3-diene");

        for (oxo_order, expected) in [
            (1, "cyclohexa-2,4-dien-1-ol"),
            (2, "cyclohexa-2,4-dien-1-one"),
        ] {
            let mut m = Molecule::default();
            for _ in 0..6 {
                m.atoms.push(Atom::new(vec2(0., 0.), "C"));
            }
            m.atoms.push(Atom::new(vec2(0., 0.), "O"));
            for i in 0..6 {
                let order = if i == 1 || i == 3 { 2 } else { 1 };
                bond(&mut m, i, (i + 1) % 6, order);
            }
            bond(&mut m, 0, 6, oxo_order);
            assert_eq!(name(&m).unwrap().name, expected);
        }
    }
    #[test]
    fn names_alkane_aldehyde_and_alkene_aldehyde_substituents() {
        let mut m = Molecule::default();
        for _ in 0..7 {
            m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        for i in 0..6 {
            bond(&mut m, i, i + 1, 1);
        }
        let o = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "O"));
        bond(&mut m, 0, o, 2);
        m.bonds
            .iter_mut()
            .find(|b| (b.a == 3 && b.b == 4) || (b.a == 4 && b.b == 3))
            .unwrap()
            .order = 2;
        let ethyl = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        let ethyl2 = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        bond(&mut m, 3, ethyl, 1);
        bond(&mut m, ethyl, ethyl2, 1);
        let oh = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "O"));
        bond(&mut m, 4, oh, 1);
        let oxo = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "O"));
        bond(&mut m, 2, oxo, 2);
        let result = name(&m).unwrap();
        assert_eq!(result.name, "4-ethyl-5-hydroxy-3-oxohept-4-enal");
        assert!(result.steps.iter().any(|s| s.contains("aldehyde")));
        assert!(result.steps.iter().any(|s| s.contains("4-ethyl")));
        assert_eq!(result.highlighted_atoms.len(), 7);
    }

    #[test]
    fn names_primary_amino_prefix_once_and_refuses_amide_as_aldehyde() {
        let mut amine = Molecule::default();
        for _ in 0..7 {
            amine.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        for i in 0..6 {
            bond(&mut amine, i, i + 1, 1);
        }
        let oxygen = amine.atoms.len();
        amine.atoms.push(Atom::new(vec2(0., 0.), "O"));
        bond(&mut amine, 0, oxygen, 2);
        let nitrogen = amine.atoms.len();
        amine.atoms.push(Atom::new(vec2(0., 0.), "N"));
        bond(&mut amine, 2, nitrogen, 1);
        assert_eq!(name(&amine).unwrap().name, "3-aminoheptanal");
        amine.atoms[nitrogen].charge = 1;
        assert!(name(&amine).unwrap_err().contains("charge"));

        let mut amide = Molecule::default();
        amide.atoms.push(Atom::new(vec2(0., 0.), "C"));
        amide.atoms.push(Atom::new(vec2(0., 0.), "O"));
        amide.atoms.push(Atom::new(vec2(0., 0.), "N"));
        bond(&mut amide, 0, 1, 2);
        bond(&mut amide, 0, 2, 1);
        assert!(name(&amide)
            .unwrap_err()
            .contains("unsupported heteroatom group"));
    }

    #[test]
    fn names_branched_saturated_alkyl_side_groups() {
        let mut m = Molecule::default();
        // Heptanal numbered from the aldehyde carbon, with a propan-2-yl group at C4.
        for _ in 0..7 {
            m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        for i in 0..6 {
            bond(&mut m, i, i + 1, 1);
        }
        let oxygen = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "O"));
        bond(&mut m, 0, oxygen, 2);

        // A three-carbon branch attached at its middle carbon.
        let middle = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        let left = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        let right = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        bond(&mut m, 3, middle, 1);
        bond(&mut m, middle, left, 1);
        bond(&mut m, middle, right, 1);

        let result = name(&m).unwrap();
        assert_eq!(result.name, "4-(propan-2-yl)heptanal");
        assert!(result.steps.iter().any(|step| step.contains("propan-2-yl")));
    }

    #[test]
    fn names_branched_butyl_substituent_and_refuses_functionalized_branch() {
        let mut m = Molecule::default();
        for _ in 0..7 {
            m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        for i in 0..6 {
            bond(&mut m, i, i + 1, 1);
        }
        let oxygen = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "O"));
        bond(&mut m, 0, oxygen, 2);

        // 2-methylpropan-1-yl at position 4.
        let root = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        let middle = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        let end = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        let methyl = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        bond(&mut m, 3, root, 1);
        bond(&mut m, root, middle, 1);
        bond(&mut m, middle, end, 1);
        bond(&mut m, middle, methyl, 1);
        assert_eq!(name(&m).unwrap().name, "4-(2-methylpropan-1-yl)heptanal");

        let hydroxy = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "O"));
        bond(&mut m, methyl, hydroxy, 1);
        assert!(name(&m)
            .unwrap_err()
            .contains("functionalized alkyl substituent"));
    }

    #[test]
    fn numbering_uses_lowest_substituent_locants_after_principal_group_and_unsaturation() {
        let mut m = Molecule::default();
        // Atom indices enumerate this hexane chain from the less favorable end.
        let chain: Vec<usize> = (0..6).rev().collect();
        for _ in 0..8 {
            m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        m.atoms.push(Atom::new(vec2(0., 0.), "O"));
        for pair in chain.windows(2) {
            bond(&mut m, pair[0], pair[1], 1);
        }
        // The alcohol sits at position 3 either way. Prefix locants decide the direction.
        bond(&mut m, chain[2], 8, 1);
        bond(&mut m, chain[1], 6, 1);
        bond(&mut m, chain[3], 7, 1);
        assert_eq!(name(&m).unwrap().name, "2,4-dimethylhexan-3-ol");
    }

    #[test]
    fn names_carboxylic_acid_and_alkene_acid_with_substituents() {
        let mut m = Molecule::default();
        for _ in 0..5 {
            m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        for i in 0..4 {
            bond(&mut m, i, i + 1, 1);
        }
        // Put the carboxyl carbon at the chain end; the chain must number from this end.
        let carbonyl_o = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "O"));
        bond(&mut m, 0, carbonyl_o, 2);
        let hydroxyl_o = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "O"));
        bond(&mut m, 0, hydroxyl_o, 1);
        m.bonds
            .iter_mut()
            .find(|b| (b.a == 2 && b.b == 3) || (b.a == 3 && b.b == 2))
            .unwrap()
            .order = 2;
        let branch = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        bond(&mut m, 3, branch, 1);

        let result = name(&m).unwrap();
        assert_eq!(result.name, "4-methylpent-3-enoic acid");
        assert!(result.steps[0].contains("carboxylic acid"));
        assert_eq!(result.highlighted_atoms.len(), 5);
    }

    #[test]
    fn names_saturated_carboxylic_acid() {
        let mut m = Molecule::default();
        for _ in 0..3 {
            m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        bond(&mut m, 0, 1, 1);
        bond(&mut m, 1, 2, 1);
        let carbonyl_o = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "O"));
        bond(&mut m, 0, carbonyl_o, 2);
        let hydroxyl_o = m.atoms.len();
        m.atoms.push(Atom::new(vec2(0., 0.), "O"));
        bond(&mut m, 0, hydroxyl_o, 1);
        assert_eq!(name(&m).unwrap().name, "propanoic acid");
    }

    #[test]
    fn names_terminal_nitriles_and_unsaturation_in_the_parent_chain() {
        let mut saturated = Molecule::default();
        for element in ["C", "C", "N"] {
            saturated.atoms.push(Atom::new(vec2(0., 0.), element));
        }
        bond(&mut saturated, 0, 1, 1);
        bond(&mut saturated, 1, 2, 3);
        let name_result = name(&saturated).unwrap();
        assert_eq!(name_result.name, "ethanenitrile");
        assert!(name_result.steps[0].contains("nitrile"));

        let mut unsaturated = Molecule::default();
        for element in ["C", "C", "C", "N"] {
            unsaturated.atoms.push(Atom::new(vec2(0., 0.), element));
        }
        bond(&mut unsaturated, 0, 1, 1);
        bond(&mut unsaturated, 1, 2, 2);
        bond(&mut unsaturated, 0, 3, 3);
        assert_eq!(name(&unsaturated).unwrap().name, "prop-2-enenitrile");

        let mut mixed = unsaturated.clone();
        mixed.atoms.push(Atom::new(vec2(0., 0.), "O"));
        bond(&mut mixed, 1, 4, 1);
        assert!(name(&mixed).is_err());
    }

    #[test]
    fn names_straight_chain_alkyl_alkanoate_esters() {
        let mut methyl_butanoate = Molecule::default();
        for _ in 0..5 {
            methyl_butanoate.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        methyl_butanoate.atoms.push(Atom::new(vec2(0., 0.), "O"));
        methyl_butanoate.atoms.push(Atom::new(vec2(0., 0.), "O"));
        for carbon in 0..3 {
            bond(&mut methyl_butanoate, carbon, carbon + 1, 1);
        }
        bond(&mut methyl_butanoate, 0, 5, 2);
        bond(&mut methyl_butanoate, 0, 6, 1);
        bond(&mut methyl_butanoate, 6, 4, 1);
        let result = name(&methyl_butanoate).unwrap();
        assert_eq!(result.name, "methyl butanoate");
        assert!(result.steps[0].contains("ester"));
        assert_eq!(result.highlighted_atoms.len(), 5);

        let mut ethyl_ethanoate = Molecule::default();
        for _ in 0..4 {
            ethyl_ethanoate.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        ethyl_ethanoate.atoms.push(Atom::new(vec2(0., 0.), "O"));
        ethyl_ethanoate.atoms.push(Atom::new(vec2(0., 0.), "O"));
        bond(&mut ethyl_ethanoate, 0, 1, 1);
        bond(&mut ethyl_ethanoate, 0, 4, 2);
        bond(&mut ethyl_ethanoate, 0, 5, 1);
        bond(&mut ethyl_ethanoate, 5, 2, 1);
        bond(&mut ethyl_ethanoate, 2, 3, 1);
        assert_eq!(name(&ethyl_ethanoate).unwrap().name, "ethyl ethanoate");
    }

    #[test]
    fn ester_namer_refuses_internal_alkyl_attachment() {
        let mut ester = Molecule::default();
        for _ in 0..5 {
            ester.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        ester.atoms.push(Atom::new(vec2(0., 0.), "O"));
        ester.atoms.push(Atom::new(vec2(0., 0.), "O"));
        bond(&mut ester, 0, 1, 1);
        bond(&mut ester, 0, 5, 2);
        bond(&mut ester, 0, 6, 1);
        bond(&mut ester, 6, 2, 1);
        bond(&mut ester, 2, 3, 1);
        bond(&mut ester, 2, 4, 1);
        assert!(name(&ester)
            .unwrap_err()
            .contains("internal alkyl attachment locants"));
    }

    #[test]
    fn ester_namer_refuses_branched_acyl_parent() {
        let mut ester = Molecule::default();
        for _ in 0..5 {
            ester.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        ester.atoms.push(Atom::new(vec2(0., 0.), "O"));
        ester.atoms.push(Atom::new(vec2(0., 0.), "O"));
        bond(&mut ester, 0, 1, 1);
        bond(&mut ester, 1, 2, 1);
        bond(&mut ester, 1, 3, 1);
        bond(&mut ester, 0, 5, 2);
        bond(&mut ester, 0, 6, 1);
        bond(&mut ester, 6, 4, 1);
        assert!(name(&ester)
            .unwrap_err()
            .contains("Branched ester carbon chains"));
    }
}
