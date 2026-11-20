use crate::{BondStyle, Molecule};
use std::collections::BTreeMap;

pub const SYMBOLS: [&str; 118] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
    "Fl", "Mc", "Lv", "Ts", "Og",
];
const ELECTRON_MASS_U: f64 = 0.000_548_579_909_065;

pub fn valid_element(s: &str) -> bool {
    SYMBOLS.contains(&s)
}
fn standard_mass(s: &str) -> Option<f64> {
    // CIAAW 2024 standard atomic weights. For interval weights, this table
    // uses the midpoint as a conventional single-value estimate. Explicit
    // isotope labels use the isotope mass below; elements without a standard
    // atomic weight remain unresolved rather than using a guessed mass number.
    Some(match s {
        "H" => 1.007975,
        "He" => 4.002602,
        "Li" => 6.9675,
        "Be" => 9.0121831,
        "B" => 10.8135,
        "C" => 12.0106,
        "N" => 14.006855,
        "O" => 15.9994,
        "F" => 18.998403,
        "Ne" => 20.1797,
        "Na" => 22.98976928,
        "Mg" => 24.3055,
        "Al" => 26.9815384,
        "Si" => 28.085,
        "P" => 30.973762,
        "S" => 32.0675,
        "Cl" => 35.4515,
        "Ar" => 39.8775,
        "Br" => 79.904,
        "I" => 126.90447,
        "K" => 39.0983,
        "Ca" => 40.078,
        "Sc" => 44.955907,
        "Ti" => 47.867,
        "V" => 50.9415,
        "Cr" => 51.9961,
        "Mn" => 54.938043,
        "Fe" => 55.845,
        "Co" => 58.933194,
        "Ni" => 58.6934,
        "Cu" => 63.546,
        "Zn" => 65.38,
        "Ga" => 69.723,
        "Ge" => 72.630,
        "As" => 74.921595,
        "Se" => 78.971,
        "Kr" => 83.798,
        "Rb" => 85.4678,
        "Sr" => 87.62,
        "Y" => 88.905838,
        "Zr" => 91.222,
        "Nb" => 92.90637,
        "Mo" => 95.95,
        "Ru" => 101.07,
        "Rh" => 102.90549,
        "Pd" => 106.42,
        "Ag" => 107.8682,
        "Cd" => 112.414,
        "In" => 114.818,
        "Sn" => 118.710,
        "Sb" => 121.760,
        "Te" => 127.60,
        "Xe" => 131.293,
        "Cs" => 132.90545196,
        "Ba" => 137.327,
        "La" => 138.90547,
        "Ce" => 140.116,
        "Pr" => 140.90766,
        "Nd" => 144.242,
        "Sm" => 150.36,
        "Eu" => 151.964,
        "Gd" => 157.249,
        "Tb" => 158.925354,
        "Dy" => 162.500,
        "Ho" => 164.930329,
        "Er" => 167.259,
        "Tm" => 168.934219,
        "Yb" => 173.045,
        "Lu" => 174.96669,
        "Hf" => 178.486,
        "Ta" => 180.94788,
        "W" => 183.84,
        "Re" => 186.207,
        "Os" => 190.23,
        "Ir" => 192.217,
        "Pt" => 195.084,
        "Au" => 196.96657,
        "Hg" => 200.592,
        "Tl" => 204.3835,
        "Pb" => 207.04,
        "Bi" => 208.98040,
        "Th" => 232.0377,
        "Pa" => 231.03588,
        "U" => 238.02891,
        _ => return None,
    })
}
fn isotope_mass(s: &str, n: u16) -> Option<f64> {
    Some(match (s, n) {
        ("H", 1) => 1.00782503223,
        ("H", 2) => 2.01410177812,
        ("H", 3) => 3.01604928199,
        ("C", 12) => 12.,
        ("C", 13) => 13.00335483507,
        ("N", 14) => 14.00307400443,
        ("N", 15) => 15.00010889888,
        ("O", 16) => 15.99491461957,
        ("O", 17) => 16.99913175650,
        ("O", 18) => 17.99915961286,
        ("S", 32) => 31.9720711744,
        ("Cl", 35) => 34.968852682,
        ("Br", 79) => 78.9183376,
        _ => return None,
    })
}
fn mono_mass(s: &str) -> Option<f64> {
    Some(match s {
        "H" => 1.00782503223,
        "C" => 12.,
        "N" => 14.00307400443,
        "O" => 15.99491461957,
        "F" => 18.998403163,
        "P" => 30.973761998,
        "S" => 31.9720711744,
        "Cl" => 34.968852682,
        "Br" => 78.9183376,
        "I" => 126.904468,
        _ => return None,
    })
}

#[derive(Debug, Clone)]
pub struct Analysis {
    pub formula: String,
    pub composition: Vec<(String, usize)>,
    pub formal_charge: i32,
    pub atom_count: usize,
    pub bond_count: usize,
    pub heavy_atom_count: usize,
    pub connected_components: usize,
    pub cycle_rank: usize,
    pub average_degree: f64,
    pub molar_mass: Option<f64>,
    pub exact_mass: Option<f64>,
    pub warnings: Vec<String>,
    pub implicit_h: Vec<Option<u8>>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct TopologicalDescriptors {
    pub heavy_atoms: usize,
    pub heavy_atom_bonds: usize,
    pub wiener_index: usize,
    pub diameter: usize,
    pub first_zagreb: usize,
    pub second_zagreb: usize,
    pub randic_connectivity: f64,
}

/// Calculate classical graph descriptors on the connected heavy-atom graph.
/// Explicit hydrogens are excluded. This is a topological calculation, not a
/// prediction of potency, shape, energy, or any experimental property.
pub fn topological_descriptors(m: &Molecule) -> Result<TopologicalDescriptors, String> {
    if m.atoms.is_empty() {
        return Err("Draw one molecule before calculating graph descriptors.".into());
    }
    if m.atoms.len() > 2048 {
        return Err("Local all-pairs graph descriptors are limited to 2,048 atoms.".into());
    }
    if m.atoms.iter().any(|atom| {
        !valid_element(&atom.element)
            || atom.query.is_some()
            || atom.r_group.is_some()
            || atom.alias.is_some()
    }) {
        return Err(
            "Graph descriptors require standard atoms without query, R-group, or alias labels."
                .into(),
        );
    }
    let heavy: Vec<usize> = m
        .atoms
        .iter()
        .enumerate()
        .filter_map(|(index, atom)| (atom.element != "H").then_some(index))
        .collect();
    if heavy.is_empty() {
        return Err("The molecule has no heavy atoms.".into());
    }
    let mut remap = vec![None; m.atoms.len()];
    for (local, &global) in heavy.iter().enumerate() {
        remap[global] = Some(local);
    }
    let mut adjacency = vec![Vec::<usize>::new(); heavy.len()];
    let mut seen_edges = std::collections::BTreeSet::new();
    for bond in &m.bonds {
        if bond.a >= m.atoms.len() || bond.b >= m.atoms.len() || bond.a == bond.b {
            return Err("The molecular graph contains an invalid bond reference.".into());
        }
        if !(1..=4).contains(&bond.order)
            || bond.query.is_some()
            || bond.style == BondStyle::Any
            || bond.style == BondStyle::Coordinate
        {
            return Err("Graph descriptors require defined covalent bonds; query and coordination bonds are unsupported.".into());
        }
        let (Some(a), Some(b)) = (remap[bond.a], remap[bond.b]) else {
            continue;
        };
        let edge = (a.min(b), a.max(b));
        if !seen_edges.insert(edge) {
            return Err("The molecular graph contains a duplicate heavy-atom bond.".into());
        }
        adjacency[a].push(b);
        adjacency[b].push(a);
    }
    if adjacency[0].is_empty() && heavy.len() > 1 {
        return Err("Graph descriptors require a single connected heavy-atom component.".into());
    }
    let mut degrees = Vec::with_capacity(heavy.len());
    for neighbors in &adjacency {
        degrees.push(neighbors.len());
    }
    let first_zagreb = degrees.iter().map(|degree| degree * degree).sum();
    let mut second_zagreb = 0usize;
    let mut randic_connectivity = 0.0;
    for &(a, b) in &seen_edges {
        second_zagreb += degrees[a] * degrees[b];
        if degrees[a] > 0 && degrees[b] > 0 {
            randic_connectivity += 1.0 / ((degrees[a] * degrees[b]) as f64).sqrt();
        }
    }
    let mut wiener_index = 0usize;
    let mut diameter = 0usize;
    for source in 0..heavy.len() {
        let mut distances = vec![usize::MAX; heavy.len()];
        distances[source] = 0;
        let mut queue = std::collections::VecDeque::from([source]);
        while let Some(current) = queue.pop_front() {
            for &neighbor in &adjacency[current] {
                if distances[neighbor] == usize::MAX {
                    distances[neighbor] = distances[current] + 1;
                    queue.push_back(neighbor);
                }
            }
        }
        if distances.contains(&usize::MAX) {
            return Err(
                "Graph descriptors require a single connected heavy-atom component.".into(),
            );
        }
        for &distance in distances.iter().skip(source + 1) {
            wiener_index += distance;
            diameter = diameter.max(distance);
        }
    }
    Ok(TopologicalDescriptors {
        heavy_atoms: heavy.len(),
        heavy_atom_bonds: seen_edges.len(),
        wiener_index,
        diameter,
        first_zagreb,
        second_zagreb,
        randic_connectivity,
    })
}

fn target_valence(element: &str, charge: i8, aromatic: bool) -> Option<i16> {
    if aromatic && element == "C" {
        return Some(3);
    }
    Some(match (element, charge) {
        ("C", 0) => 4,
        ("N", 0) => 3,
        ("N", 1) => 4,
        ("O", 0) => 2,
        ("O", -1) => 1,
        ("O", 1) => 3,
        ("S", 0) => 2,
        ("P", 0) => 3,
        ("B", 0) => 3,
        ("F" | "Cl" | "Br" | "I", 0) => 1,
        ("H", 0) => 1,
        _ => return None,
    })
}

pub fn analyze(m: &Molecule) -> Analysis {
    let mut warnings = Vec::new();
    let mut counts: BTreeMap<String, usize> = BTreeMap::new();
    let mut mass = Some(0.);
    let mut exact = Some(0.);
    let mut added_h = 0usize;
    let mut implicit = Vec::with_capacity(m.atoms.len());
    // Build per-atom degree and bond valence once. The inspector calls analyze
    // during rendering, so scanning every bond for every atom made large
    // structures quadratic and could freeze the UI.
    let mut degrees = vec![0i16; m.atoms.len()];
    let mut valences = vec![0i16; m.atoms.len()];
    for (bond_index, b) in m.bonds.iter().enumerate() {
        if b.stereo.is_some() && b.order != 2 {
            warnings.push(format!(
                "Bond {} has E/Z stereochemistry but is not a double bond",
                bond_index + 1
            ));
        }
        if matches!(
            b.style,
            BondStyle::Wedge | BondStyle::Hash | BondStyle::Wavy
        ) && b.order != 1
        {
            warnings.push(format!(
                "Bond {} has a wedge/hash/wavy display mark but is not a single bond",
                bond_index + 1
            ));
        }
        let contribution = if b.style == BondStyle::Coordinate {
            0
        } else if b.style == BondStyle::Aromatic {
            1
        } else {
            b.order as i16
        };
        for index in [b.a, b.b] {
            if let (Some(degree), Some(valence)) = (degrees.get_mut(index), valences.get_mut(index))
            {
                *degree += 1;
                *valence += contribution;
            }
        }
    }
    let mut adjacency = vec![Vec::new(); m.atoms.len()];
    for bond in &m.bonds {
        if bond.a < adjacency.len() && bond.b < adjacency.len() {
            adjacency[bond.a].push(bond.b);
            adjacency[bond.b].push(bond.a);
        }
    }
    let mut visited = vec![false; m.atoms.len()];
    let mut connected_components = 0usize;
    for root in 0..m.atoms.len() {
        if visited[root] {
            continue;
        }
        connected_components += 1;
        visited[root] = true;
        let mut stack = vec![root];
        while let Some(atom) = stack.pop() {
            for &neighbor in &adjacency[atom] {
                if !visited[neighbor] {
                    visited[neighbor] = true;
                    stack.push(neighbor);
                }
            }
        }
    }
    for (i, a) in m.atoms.iter().enumerate() {
        if !valid_element(&a.element) {
            warnings.push(format!(
                "Atom {} has an unknown or pseudo element: {}",
                i + 1,
                a.element
            ));
        }
        let val = valences[i];
        let degree = degrees[i];
        let target = a
            .valence_override
            .map(i16::from)
            .or_else(|| target_valence(&a.element, a.charge, a.aromatic));
        let ih = if let Some(x) = a.implicit_h_override {
            Some(x)
        } else if a.element == "H" {
            Some(0)
        } else if a.radical > 0 || a.query.is_some() || a.r_group.is_some() {
            None
        } else if a.valence_override.is_some() {
            target.map(|target| (target - val - a.explicit_h as i16).clamp(0, 8) as u8)
        } else if a.aromatic && a.element == "C" {
            Some((3 - degree).max(0) as u8)
        } else if a.aromatic && a.element == "N" {
            Some(0)
        } else {
            target
                .map(|target| (target - val - a.explicit_h as i16).clamp(0, 8) as u8)
                .or_else(|| {
                    if valid_element(&a.element) {
                        Some(0)
                    } else {
                        None
                    }
                })
        };
        if let Some(target) = target {
            if val + a.explicit_h as i16 > target {
                warnings.push(format!(
                    "Atom {} ({}): bond valence {} exceeds typical valence {}",
                    i + 1,
                    a.element,
                    val + a.explicit_h as i16,
                    target
                ));
            }
        }
        implicit.push(ih);
        *counts.entry(a.element.clone()).or_default() += 1;
        if let Some(h) = ih {
            added_h += h as usize;
            if h > 0 {
                *counts.entry("H".into()).or_default() += h as usize;
            }
        } else {
            mass = None;
            exact = None;
            warnings.push(format!(
                "Atom {}: implicit hydrogen count is unresolved",
                i + 1
            ));
        }
        if a.explicit_h > 0 {
            added_h += a.explicit_h as usize;
            *counts.entry("H".into()).or_default() += a.explicit_h as usize;
        }
        mass = mass.and_then(|x| {
            if let Some(n) = a.isotope {
                isotope_mass(&a.element, n)
            } else {
                standard_mass(&a.element)
            }
            .map(|v| x + v)
        });
        exact = exact.and_then(|x| {
            if let Some(n) = a.isotope {
                isotope_mass(&a.element, n)
            } else {
                mono_mass(&a.element)
            }
            .map(|v| x + v)
        });
    }
    mass = mass.and_then(|x| standard_mass("H").map(|v| x + v * added_h as f64));
    exact = exact.and_then(|x| mono_mass("H").map(|v| x + v * added_h as f64));
    let mut ordered = Vec::new();
    if let Some(n) = counts.remove("C") {
        ordered.push(("C".to_string(), n));
    }
    if let Some(n) = counts.remove("H") {
        ordered.push(("H".to_string(), n));
    }
    ordered.extend(counts);
    let formula = ordered
        .iter()
        .map(|(el, n)| {
            if *n == 1 {
                el.clone()
            } else {
                format!("{el}{n}")
            }
        })
        .collect::<String>();
    let composition = ordered;
    let formal_charge = m.atoms.iter().map(|a| a.charge as i32).sum();
    // Exact mass of an ionic formula includes the electron deficit/excess.
    // A positive charge removes electrons; a negative charge adds them.
    exact = exact.map(|neutral_mass| neutral_mass - formal_charge as f64 * ELECTRON_MASS_U);
    let atom_count = m.atoms.len();
    let bond_count = m.bonds.len();
    let cycle_rank = bond_count
        .saturating_add(connected_components)
        .saturating_sub(atom_count);
    let average_degree = if atom_count == 0 {
        0.0
    } else {
        degrees.iter().map(|degree| *degree as f64).sum::<f64>() / atom_count as f64
    };
    let heavy_atom_count = m.atoms.iter().filter(|a| a.element != "H").count();
    Analysis {
        formula,
        composition,
        formal_charge,
        atom_count,
        bond_count,
        heavy_atom_count,
        connected_components,
        cycle_rank,
        average_degree,
        molar_mass: mass,
        exact_mass: exact,
        warnings,
        implicit_h: implicit,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Atom, Bond, BondStyle};
    use macroquad::prelude::vec2;
    #[test]
    fn ethanol_formula_and_mass() {
        let mut m = Molecule::default();
        for x in ["C", "C", "O"] {
            m.atoms.push(Atom::new(vec2(0., 0.), x));
        }
        for (a, b) in [(0, 1), (1, 2)] {
            m.bonds.push(Bond {
                a,
                b,
                order: 1,
                style: BondStyle::Normal,
                stereo: None,
                query: None,
            });
        }
        let a = analyze(&m);
        assert_eq!(a.formula, "C2H6O");
        assert!((a.molar_mass.unwrap() - 46.069).abs() < 0.01);
    }

    #[test]
    fn topological_descriptors_match_known_path_graph_values_and_ignore_hydrogens() {
        let mut molecule = Molecule::default();
        for element in ["C", "C", "O", "H"] {
            molecule.atoms.push(Atom::new(vec2(0., 0.), element));
        }
        for (a, b) in [(0, 1), (1, 2), (0, 3)] {
            molecule.bonds.push(Bond {
                a,
                b,
                order: 1,
                style: BondStyle::Normal,
                stereo: None,
                query: None,
            });
        }
        let result = topological_descriptors(&molecule).unwrap();
        assert_eq!(result.heavy_atoms, 3);
        assert_eq!(result.heavy_atom_bonds, 2);
        assert_eq!(result.wiener_index, 4);
        assert_eq!(result.diameter, 2);
        assert_eq!(result.first_zagreb, 6);
        assert_eq!(result.second_zagreb, 4);
        assert!((result.randic_connectivity - 2.0_f64.sqrt()).abs() < 1e-12);
    }

    #[test]
    fn topological_descriptors_match_known_benzene_cycle_graph_values() {
        let mut molecule = Molecule::default();
        for _ in 0..6 {
            molecule.atoms.push(Atom::new(vec2(0., 0.), "C"));
        }
        for index in 0..6 {
            molecule.bonds.push(Bond {
                a: index,
                b: (index + 1) % 6,
                order: 4,
                style: BondStyle::Aromatic,
                stereo: None,
                query: None,
            });
        }
        let result = topological_descriptors(&molecule).unwrap();
        assert_eq!(result.wiener_index, 27);
        assert_eq!(result.diameter, 3);
        assert_eq!(result.first_zagreb, 24);
        assert_eq!(result.second_zagreb, 24);
        assert!((result.randic_connectivity - 3.0).abs() < 1e-12);
    }

    #[test]
    fn topological_descriptors_refuse_disconnected_molecules() {
        let mut molecule = Molecule::default();
        molecule.atoms.push(Atom::new(vec2(0., 0.), "C"));
        molecule.atoms.push(Atom::new(vec2(1., 0.), "O"));
        assert!(topological_descriptors(&molecule)
            .unwrap_err()
            .contains("connected"));
    }

    #[test]
    fn standard_atomic_weights_cover_stable_metals_and_keep_gaps_explicit() {
        assert_eq!(standard_mass("Ni"), Some(58.6934));
        assert_eq!(standard_mass("Ge"), Some(72.630));
        assert_eq!(standard_mass("U"), Some(238.02891));
        assert_eq!(standard_mass("Tc"), None);
        assert_eq!(standard_mass("Og"), None);
    }

    #[test]
    fn isotope_label_uses_isotopic_mass_instead_of_standard_weight() {
        let mut molecule = Molecule::default();
        let mut carbon = Atom::new(vec2(0., 0.), "C");
        carbon.isotope = Some(13);
        carbon.implicit_h_override = Some(0);
        molecule.atoms.push(carbon);
        let result = analyze(&molecule);
        assert_eq!(result.formula, "C");
        assert!((result.molar_mass.unwrap() - 13.00335483507).abs() < 1e-10);
        assert!((result.exact_mass.unwrap() - 13.00335483507).abs() < 1e-10);
    }

    #[test]
    fn exact_ionic_mass_accounts_for_electron_gain_and_loss() {
        let mut cation = Molecule::default();
        let mut carbon = Atom::new(vec2(0., 0.), "C");
        carbon.implicit_h_override = Some(0);
        carbon.charge = 1;
        cation.atoms.push(carbon);
        let mut anion = cation.clone();
        anion.atoms[0].charge = -1;

        assert!((analyze(&cation).exact_mass.unwrap() - (12. - ELECTRON_MASS_U)).abs() < 1e-12);
        assert!((analyze(&anion).exact_mass.unwrap() - (12. + ELECTRON_MASS_U)).abs() < 1e-12);
    }
    #[test]
    fn invalid_valence_is_reported() {
        let mut m = Molecule::default();
        m.atoms.push(Atom::new(vec2(0., 0.), "C"));
        for i in 1..=5 {
            m.atoms.push(Atom::new(vec2(i as f32, 0.), "H"));
            m.bonds.push(Bond {
                a: 0,
                b: i,
                order: 1,
                style: BondStyle::Normal,
                stereo: None,
                query: None,
            });
        }
        assert!(analyze(&m).warnings.iter().any(|x| x.contains("exceeds")));
    }

    #[test]
    fn inconsistent_double_bond_stereo_and_wedge_orders_are_reported() {
        let mut molecule = Molecule::default();
        molecule.atoms.push(Atom::new(vec2(0., 0.), "C"));
        molecule.atoms.push(Atom::new(vec2(40., 0.), "C"));
        molecule.atoms.push(Atom::new(vec2(80., 0.), "C"));
        molecule.bonds.extend([
            Bond {
                a: 0,
                b: 1,
                order: 1,
                style: BondStyle::Normal,
                stereo: Some("E".into()),
                query: None,
            },
            Bond {
                a: 1,
                b: 2,
                order: 2,
                style: BondStyle::Wedge,
                stereo: None,
                query: None,
            },
        ]);

        let analysis = analyze(&molecule);
        assert!(analysis
            .warnings
            .iter()
            .any(|warning| warning.contains("E/Z stereochemistry but is not a double bond")));
        assert!(analysis
            .warnings
            .iter()
            .any(|warning| warning.contains("display mark but is not a single bond")));
    }

    #[test]
    fn atom_valence_override_controls_hydrogen_estimate_and_warning() {
        let mut mol = Molecule::default();
        let mut carbon = Atom::new(vec2(0., 0.), "C");
        carbon.valence_override = Some(3);
        mol.atoms.push(carbon);
        mol.atoms.push(Atom::new(vec2(1., 0.), "C"));
        mol.bonds.push(crate::Bond {
            a: 0,
            b: 1,
            order: 1,
            style: BondStyle::Normal,
            stereo: None,
            query: None,
        });
        let analysis = analyze(&mol);
        assert_eq!(analysis.implicit_h[0], Some(2));
        mol.atoms[0].valence_override = Some(0);
        let invalid = analyze(&mol);
        assert_eq!(invalid.implicit_h[0], Some(0));
        assert!(invalid
            .warnings
            .iter()
            .any(|warning| warning.contains("typical valence 0")));
    }

    #[test]
    fn graph_topology_counts_components_and_independent_cycles() {
        let mut m = Molecule::default();
        m.atoms
            .extend((0..4).map(|i| Atom::new(vec2(i as f32, 0.), "C")));
        for (a, b) in [(0, 1), (1, 2), (2, 0)] {
            m.bonds.push(Bond {
                a,
                b,
                order: 1,
                style: BondStyle::Normal,
                stereo: None,
                query: None,
            });
        }
        let result = analyze(&m);
        assert_eq!(result.connected_components, 2);
        assert_eq!(result.cycle_rank, 1);
        assert_eq!(result.average_degree, 1.5);
    }

    #[test]
    fn analyzes_large_chain_without_quadratic_bond_scans() {
        let n = 20_000usize;
        let mut m = Molecule::default();
        m.atoms
            .extend((0..n).map(|i| Atom::new(vec2(i as f32, 0.), "C")));
        m.bonds.extend((1..n).map(|i| Bond {
            a: i - 1,
            b: i,
            order: 1,
            style: BondStyle::Normal,
            stereo: None,
            query: None,
        }));
        let a = analyze(&m);
        assert_eq!(a.atom_count, n);
        assert_eq!(a.bond_count, n - 1);
        assert_eq!(a.implicit_h[0], Some(3));
        assert_eq!(a.implicit_h[n / 2], Some(2));
        assert_eq!(a.implicit_h[n - 1], Some(3));
    }
}
