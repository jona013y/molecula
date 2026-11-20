use crate::{chemistry, Graphic, Molecule};
use std::collections::{BTreeMap, BTreeSet};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Role {
    Reactant,
    Product,
    Agent,
}

pub struct ReactionReport {
    pub step: u16,
    pub conditions_above: String,
    pub conditions_below: String,
    pub reactants: Vec<String>,
    pub products: Vec<String>,
    pub agents: Vec<String>,
    pub balanced: Option<bool>,
    pub charge_left: i32,
    pub charge_right: i32,
    pub warnings: Vec<String>,
    pub bond_changes: Vec<String>,
    pub changed_map_pairs: Vec<(u16, u16)>,
    pub element_balance: Vec<(String, i32, i32)>,
    pub isotope_balance: Vec<(String, i32, i32)>,
    pub components: Vec<ReactionComponentRow>,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ReactionComponentRow {
    pub role: &'static str,
    pub role_source: &'static str,
    pub coefficient: u8,
    pub formula: String,
    pub formal_charge: i32,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoefficientProposal {
    pub atom_indices: Vec<usize>,
    pub coefficient: u8,
    pub formula: String,
}

fn components(m: &Molecule) -> Vec<Vec<usize>> {
    let mut seen = vec![false; m.atoms.len()];
    let mut result = Vec::new();
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
                let other = if b.a == i {
                    Some(b.b)
                } else if b.b == i {
                    Some(b.a)
                } else {
                    None
                };
                if let Some(j) = other {
                    if !seen[j] {
                        seen[j] = true;
                        stack.push(j);
                    }
                }
            }
        }
        result.push(group);
    }
    result
}

pub fn analyze_all(m: &Molecule, graphics: &[Graphic]) -> Result<Vec<ReactionReport>, String> {
    let count=graphics.iter().filter(|g| matches!(g,Graphic::Arrow{kind,..} if kind=="reaction"||kind=="equilibrium"||kind=="retrosynthesis")).count();
    if count == 0 {
        return Err("Draw a reaction arrow first".into());
    }
    (1..=count)
        .map(|step| analyze_step(m, graphics, step as u16))
        .collect()
}

pub fn analyze_step(
    m: &Molecule,
    graphics: &[Graphic],
    step: u16,
) -> Result<ReactionReport, String> {
    let (start, end, conditions_above, conditions_below) = graphics
        .iter()
        .filter_map(|g| match g {
            Graphic::Arrow {
                start,
                end,
                kind,
                above,
                below,
            } if kind == "reaction" || kind == "equilibrium" || kind == "retrosynthesis" => {
                Some((*start, *end, above.clone(), below.clone()))
            }
            _ => None,
        })
        .nth(step.saturating_sub(1) as usize)
        .ok_or("Draw a reaction, equilibrium, or retrosynthesis arrow first")?;
    let axis = end - start;
    if axis.length_squared() < 100. {
        return Err("Reaction arrow is too short".into());
    }
    let mut out = ReactionReport {
        step,
        conditions_above,
        conditions_below,
        reactants: Vec::new(),
        products: Vec::new(),
        agents: Vec::new(),
        balanced: None,
        charge_left: 0,
        charge_right: 0,
        warnings: Vec::new(),
        bond_changes: Vec::new(),
        changed_map_pairs: Vec::new(),
        element_balance: Vec::new(),
        isotope_balance: Vec::new(),
        components: Vec::new(),
    };
    let mut left = BTreeMap::<String, usize>::new();
    let mut right = BTreeMap::<String, usize>::new();
    let mut left_isotopes = BTreeMap::<String, usize>::new();
    let mut right_isotopes = BTreeMap::<String, usize>::new();
    let mut left_maps = BTreeSet::new();
    let mut right_maps = BTreeSet::new();
    let mut left_bonds = BTreeMap::<(u16, u16), u8>::new();
    let mut right_bonds = BTreeMap::<(u16, u16), u8>::new();
    let mut left_elements = BTreeMap::<u16, String>::new();
    let mut right_elements = BTreeMap::<u16, String>::new();
    let mut reliable = true;
    for group in components(m) {
        let explicit_steps: BTreeSet<u16> =
            group.iter().map(|&i| m.atoms[i].reaction_step).collect();
        if explicit_steps.len() == 1 && !explicit_steps.contains(&step) {
            continue;
        }
        if explicit_steps.len() > 1 {
            reliable = false;
            out.warnings
                .push("Reaction component has inconsistent step assignments".into());
            continue;
        }
        let explicit_roles: BTreeSet<String> = group
            .iter()
            .filter_map(|&i| m.atoms[i].reaction_role.clone())
            .collect();
        let coefficients: BTreeSet<u8> = group.iter().map(|&i| m.atoms[i].coefficient).collect();
        let coefficient = if coefficients.len() == 1 {
            coefficients.first().copied().unwrap_or(1) as usize
        } else {
            reliable = false;
            out.warnings
                .push("Reaction component has inconsistent coefficients; treating it as 1".into());
            1
        };
        let center = group
            .iter()
            .map(|&i| m.atoms[i].pos)
            .sum::<macroquad::prelude::Vec2>()
            / group.len() as f32;
        let t = (center - start).dot(axis) / axis.length_squared();
        let inferred = if t < 0.45 {
            Role::Reactant
        } else if t > 0.55 {
            Role::Product
        } else {
            Role::Agent
        };
        let role = if explicit_roles.len() > 1 {
            reliable = false;
            out.warnings
                .push("Reaction component has conflicting explicit roles".into());
            inferred
        } else if let Some(role) = explicit_roles.first() {
            match role.as_str() {
                "reactant" => Role::Reactant,
                "product" => Role::Product,
                "agent" => Role::Agent,
                _ => inferred,
            }
        } else {
            inferred
        };
        let role_name = match role {
            Role::Reactant => "reactant",
            Role::Product => "product",
            Role::Agent => "agent",
        };
        let role_source = if explicit_roles.len() == 1 {
            "explicit"
        } else {
            "inferred"
        };
        let mut sub = Molecule::default();
        let mut remap = BTreeMap::new();
        for &i in &group {
            remap.insert(i, sub.atoms.len());
            sub.atoms.push(m.atoms[i].clone());
        }
        for b in &m.bonds {
            if let (Some(&a), Some(&c)) = (remap.get(&b.a), remap.get(&b.b)) {
                let mut bond = b.clone();
                bond.a = a;
                bond.b = c;
                sub.bonds.push(bond);
            }
        }
        let a = chemistry::analyze(&sub);
        if !a.warnings.is_empty() || a.implicit_h.iter().any(Option::is_none) {
            reliable = false;
            out.warnings.extend(a.warnings);
        }
        let formula = if a.formula.is_empty() {
            "?".into()
        } else {
            a.formula.clone()
        };
        let component_summary = match a.molar_mass {
            Some(mass) => format!("{formula} ({mass:.3} g/mol)"),
            None => format!("{formula} (molar mass unavailable)"),
        };
        let component_summary = if coefficient > 1 {
            format!("{coefficient} × {component_summary}")
        } else {
            component_summary
        };
        out.components.push(ReactionComponentRow {
            role: role_name,
            role_source,
            coefficient: coefficient as u8,
            formula: formula.clone(),
            formal_charge: a.formal_charge,
        });
        match role {
            Role::Reactant | Role::Product => {
                let (counts, isotopes, maps, names, charge) = if role == Role::Reactant {
                    (
                        &mut left,
                        &mut left_isotopes,
                        &mut left_maps,
                        &mut out.reactants,
                        &mut out.charge_left,
                    )
                } else {
                    (
                        &mut right,
                        &mut right_isotopes,
                        &mut right_maps,
                        &mut out.products,
                        &mut out.charge_right,
                    )
                };
                names.push(component_summary);
                *charge += a.formal_charge * coefficient as i32;
                for (el, n) in a.composition {
                    *counts.entry(el).or_default() += n * coefficient;
                }
                for atom in &sub.atoms {
                    if let Some(isotope) = atom.isotope {
                        *isotopes
                            .entry(format!("{isotope}{}", atom.element))
                            .or_default() += coefficient;
                    }
                }
                let elements = if role == Role::Reactant {
                    &mut left_elements
                } else {
                    &mut right_elements
                };
                for atom in &sub.atoms {
                    if let Some(map) = atom.atom_map {
                        if !maps.insert(map) {
                            out.warnings
                                .push(format!("Atom map {map} is repeated on one side"));
                        }
                        elements.insert(map, atom.element.clone());
                    }
                }
                let mapped_bonds = if role == Role::Reactant {
                    &mut left_bonds
                } else {
                    &mut right_bonds
                };
                for bond in &sub.bonds {
                    if let (Some(a), Some(b)) =
                        (sub.atoms[bond.a].atom_map, sub.atoms[bond.b].atom_map)
                    {
                        let pair = (a.min(b), a.max(b));
                        mapped_bonds.insert(pair, bond.order);
                    }
                }
            }
            Role::Agent => out.agents.push(component_summary),
        }
    }
    if out.reactants.is_empty() || out.products.is_empty() {
        out.warnings
            .push("Place structures on both sides of the arrow".into());
    } else if reliable {
        out.balanced = Some(
            left == right && left_isotopes == right_isotopes && out.charge_left == out.charge_right,
        );
    }
    for element in left.keys().chain(right.keys()).collect::<BTreeSet<_>>() {
        out.element_balance.push((
            element.clone(),
            left.get(element).copied().unwrap_or(0) as i32,
            right.get(element).copied().unwrap_or(0) as i32,
        ));
    }
    for isotope in left_isotopes
        .keys()
        .chain(right_isotopes.keys())
        .collect::<BTreeSet<_>>()
    {
        out.isotope_balance.push((
            isotope.clone(),
            left_isotopes.get(isotope).copied().unwrap_or(0) as i32,
            right_isotopes.get(isotope).copied().unwrap_or(0) as i32,
        ));
    }
    if !left_maps.is_empty() || !right_maps.is_empty() {
        for map in left_maps.symmetric_difference(&right_maps) {
            out.warnings
                .push(format!("Atom map {map} occurs on only one side"));
        }
        for (map, element) in &left_elements {
            if let Some(other) = right_elements.get(map) {
                if element != other {
                    out.warnings.push(format!(
                        "Atom map {map} changes element: {element} → {other}"
                    ));
                }
            }
        }
        for pair in left_bonds
            .keys()
            .chain(right_bonds.keys())
            .collect::<BTreeSet<_>>()
        {
            let before = left_bonds.get(pair).copied().unwrap_or(0);
            let after = right_bonds.get(pair).copied().unwrap_or(0);
            if before == after {
                continue;
            }
            let label = if before == 0 {
                format!("{}-{} formed (order {after})", pair.0, pair.1)
            } else if after == 0 {
                format!("{}-{} broken (order {before})", pair.0, pair.1)
            } else {
                format!("{}-{} order {before} → {after}", pair.0, pair.1)
            };
            out.bond_changes.push(label);
            out.changed_map_pairs.push(*pair);
        }
    }
    out.warnings.dedup();
    Ok(out)
}

/// Calculate positive integer stoichiometric coefficients for a single reaction step.
/// The returned atom indices let the editor apply each coefficient to its entire component.
pub fn propose_balancing_coefficients(
    m: &Molecule,
    graphics: &[Graphic],
    step: u16,
) -> Result<Vec<CoefficientProposal>, String> {
    let (start, end) = graphics
        .iter()
        .filter_map(|graphic| match graphic {
            Graphic::Arrow {
                start, end, kind, ..
            } if kind == "reaction" || kind == "equilibrium" || kind == "retrosynthesis" => {
                Some((*start, *end))
            }
            _ => None,
        })
        .nth(step.saturating_sub(1) as usize)
        .ok_or("Draw a reaction arrow before balancing")?;
    let axis = end - start;
    if axis.length_squared() < 100. {
        return Err("Reaction arrow is too short to classify its sides".into());
    }

    struct Column {
        atom_indices: Vec<usize>,
        role: Role,
        composition: BTreeMap<String, i64>,
        formula: String,
    }
    let mut columns = Vec::new();
    for group in components(m) {
        let steps: BTreeSet<u16> = group.iter().map(|&i| m.atoms[i].reaction_step).collect();
        if steps.len() > 1 {
            return Err("A component has conflicting reaction-step assignments".into());
        }
        if steps.first().is_some_and(|&assigned| assigned != step) {
            continue;
        }
        let roles: BTreeSet<String> = group
            .iter()
            .filter_map(|&i| m.atoms[i].reaction_role.clone())
            .collect();
        if roles.len() > 1 {
            return Err("A component has conflicting reactant/product/agent roles".into());
        }
        let center = group
            .iter()
            .map(|&i| m.atoms[i].pos)
            .sum::<macroquad::prelude::Vec2>()
            / group.len() as f32;
        let t = (center - start).dot(axis) / axis.length_squared();
        let inferred = if t < 0.45 {
            Role::Reactant
        } else if t > 0.55 {
            Role::Product
        } else {
            Role::Agent
        };
        let role = roles
            .first()
            .map(|role| match role.as_str() {
                "reactant" => Role::Reactant,
                "product" => Role::Product,
                "agent" => Role::Agent,
                _ => inferred,
            })
            .unwrap_or(inferred);
        if role == Role::Agent {
            continue;
        }
        let mut sub = Molecule::default();
        let mut remap = BTreeMap::new();
        for &index in &group {
            let mut atom = m.atoms[index].clone();
            atom.coefficient = 1;
            remap.insert(index, sub.atoms.len());
            sub.atoms.push(atom);
        }
        for bond in &m.bonds {
            if let (Some(&a), Some(&b)) = (remap.get(&bond.a), remap.get(&bond.b)) {
                let mut copied = bond.clone();
                copied.a = a;
                copied.b = b;
                sub.bonds.push(copied);
            }
        }
        let analysis = chemistry::analyze(&sub);
        if !analysis.warnings.is_empty() || analysis.implicit_h.iter().any(Option::is_none) {
            return Err(format!(
                "Cannot balance component {}: chemistry is unresolved ({})",
                columns.len() + 1,
                analysis.warnings.join("; ")
            ));
        }
        let mut composition: BTreeMap<String, i64> = analysis
            .composition
            .into_iter()
            .map(|(element, count)| (element, count as i64))
            .collect();
        for atom in &sub.atoms {
            if let Some(isotope) = atom.isotope {
                *composition
                    .entry(format!("__isotope_{}_{}", atom.element, isotope))
                    .or_default() += 1;
            }
        }
        composition.insert("__formal_charge".into(), analysis.formal_charge as i64);
        columns.push(Column {
            atom_indices: group,
            role,
            composition,
            formula: analysis.formula,
        });
    }
    if !columns.iter().any(|column| column.role == Role::Reactant)
        || !columns.iter().any(|column| column.role == Role::Product)
    {
        return Err("Mark or place reactants and products on opposite sides first".into());
    }
    if columns.len() < 2 {
        return Err("A reaction needs at least one reactant and one product".into());
    }

    let rows: BTreeSet<String> = columns
        .iter()
        .flat_map(|column| column.composition.keys().cloned())
        .collect();
    let matrix: Vec<Vec<f64>> = rows
        .iter()
        .map(|row| {
            columns
                .iter()
                .map(|column| {
                    let count = *column.composition.get(row).unwrap_or(&0) as f64;
                    if column.role == Role::Reactant {
                        count
                    } else {
                        -count
                    }
                })
                .collect()
        })
        .collect();
    let mut rref = matrix.clone();
    let mut pivot_columns = Vec::new();
    let mut pivot_row = 0;
    for column in 0..columns.len() {
        let Some(found) = (pivot_row..rref.len()).find(|&row| rref[row][column].abs() > 1e-9)
        else {
            continue;
        };
        rref.swap(pivot_row, found);
        let pivot = rref[pivot_row][column];
        for value in &mut rref[pivot_row] {
            *value /= pivot;
        }
        let normalized = rref[pivot_row].clone();
        for (row_index, row_values) in rref.iter_mut().enumerate() {
            if row_index == pivot_row {
                continue;
            }
            let factor = row_values[column];
            if factor.abs() > 1e-9 {
                for (value, pivot_value) in row_values.iter_mut().zip(&normalized) {
                    *value -= factor * pivot_value;
                }
            }
        }
        pivot_columns.push(column);
        pivot_row += 1;
        if pivot_row == rref.len() {
            break;
        }
    }
    let free_columns: Vec<usize> = (0..columns.len())
        .filter(|column| !pivot_columns.contains(column))
        .collect();
    if free_columns.len() != 1 {
        return Err(format!(
            "The coefficient system has {} free dimensions; combine duplicate or underdetermined components before automatic balancing",
            free_columns.len()
        ));
    }
    let free = free_columns[0];
    let mut solution = vec![0.0f64; columns.len()];
    solution[free] = 1.0;
    for (row, &pivot) in pivot_columns.iter().enumerate() {
        solution[pivot] = -rref[row][free];
    }
    if solution.iter().all(|&value| value < -1e-9) {
        for value in &mut solution {
            *value = -*value;
        }
    }
    if solution.iter().any(|&value| value <= 1e-9) {
        return Err("No strictly positive coefficient exists for every component".into());
    }

    let mut rational = Vec::with_capacity(solution.len());
    for value in solution {
        let mut found = None;
        for denominator in 1i128..=10_000 {
            let numerator = (value * denominator as f64).round() as i128;
            if (numerator as f64 / denominator as f64 - value).abs() < 1e-9 {
                let divisor = gcd(numerator, denominator);
                found = Some((numerator / divisor, denominator / divisor));
                break;
            }
        }
        rational
            .push(found.ok_or("Could not represent the balance as small integer coefficients")?);
    }
    let common_denominator = rational
        .iter()
        .try_fold(1i128, |acc, &(_, denominator)| lcm(acc, denominator))
        .ok_or("Balanced coefficients are too large")?;
    let mut integer = Vec::with_capacity(rational.len());
    for &(numerator, denominator) in &rational {
        integer.push(
            numerator
                .checked_mul(common_denominator / denominator)
                .ok_or("Balanced coefficients are too large")?,
        );
    }
    let divisor = integer.iter().copied().fold(0, gcd);
    for coefficient in &mut integer {
        *coefficient /= divisor;
    }
    if integer
        .iter()
        .any(|&coefficient| !(1..=99).contains(&coefficient))
    {
        return Err(
            "The smallest positive balance exceeds the editor coefficient limit (99)".into(),
        );
    }
    for matrix_row in &matrix {
        let total: f64 = matrix_row
            .iter()
            .zip(&integer)
            .map(|(&count, &coefficient)| count * coefficient as f64)
            .sum();
        if total.abs() > 1e-7 {
            return Err(
                "The proposed integer coefficients failed exact element/charge verification".into(),
            );
        }
    }
    Ok(columns
        .into_iter()
        .zip(integer)
        .map(|(column, coefficient)| CoefficientProposal {
            atom_indices: column.atom_indices,
            coefficient: coefficient as u8,
            formula: column.formula,
        })
        .collect())
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        let remainder = a % b;
        a = b;
        b = remainder;
    }
    a.max(1)
}

fn lcm(a: i128, b: i128) -> Option<i128> {
    a.checked_div(gcd(a, b))?.checked_mul(b)
}

pub fn report_text(r: &ReactionReport) -> String {
    let side = |v: &Vec<String>| {
        if v.is_empty() {
            "(none)".into()
        } else {
            v.join(" + ")
        }
    };
    let balance = match r.balanced {
        Some(true) => "Balanced by elements and formal charge",
        Some(false) => "Not balanced by elements and formal charge",
        None => "Balance unavailable: unresolved chemistry or incomplete sides",
    };
    let mut s = format!(
        "Reaction step {}\nReactants: {}\nProducts: {}\nAgents: {}\nFormal charge: {} → {}\n{}\nRoles may be inferred from geometry unless explicitly assigned.",
        r.step,
        side(&r.reactants),
        side(&r.products),
        side(&r.agents),
        r.charge_left,
        r.charge_right,
        balance
    );
    if !r.conditions_above.trim().is_empty() || !r.conditions_below.trim().is_empty() {
        s.push_str("\nReaction conditions:");
        if !r.conditions_above.trim().is_empty() {
            s.push_str(&format!("\nAbove arrow: {}", r.conditions_above.trim()));
        }
        if !r.conditions_below.trim().is_empty() {
            s.push_str(&format!("\nBelow arrow: {}", r.conditions_below.trim()));
        }
    }
    if !r.bond_changes.is_empty() {
        s.push_str("\nMapped bond changes:\n");
        s.push_str(&r.bond_changes.join("\n"));
    }
    if !r.element_balance.is_empty() {
        s.push_str("\nStoichiometry table (reactants | products | Δ):\n");
        for (element, left, right) in &r.element_balance {
            s.push_str(&format!(
                "{element}: {left} | {right} | {:+}\n",
                right - left
            ));
        }
    }
    if !r.isotope_balance.is_empty() {
        s.push_str("\nIsotope balance (reactants | products | Δ):\n");
        for (isotope, left, right) in &r.isotope_balance {
            s.push_str(&format!(
                "{isotope}: {left} | {right} | {:+}\n",
                right - left
            ));
        }
    }
    if !r.components.is_empty() {
        s.push_str(
            "\nReaction component table (role | coefficient | formula | charge | role source):\n",
        );
        for component in &r.components {
            s.push_str(&format!(
                "{} | {} | {} | {:+} | {}\n",
                component.role,
                component.coefficient,
                component.formula,
                component.formal_charge,
                component.role_source
            ));
        }
    }
    if !r.warnings.is_empty() {
        s.push_str("\nWarnings:\n");
        s.push_str(&r.warnings.join("\n"));
    }
    s
}

pub fn report_all_text(reports: &[ReactionReport]) -> String {
    reports
        .iter()
        .map(report_text)
        .collect::<Vec<_>>()
        .join("\n\n──────────\n\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Atom;
    use macroquad::prelude::vec2;
    #[test]
    fn balances_identical_sides_and_detects_mismatch() {
        let mut m = Molecule::default();
        m.atoms.push(Atom::new(vec2(100., 100.), "C"));
        m.atoms.push(Atom::new(vec2(400., 100.), "C"));
        let arrow = vec![Graphic::Arrow {
            start: vec2(200., 100.),
            end: vec2(300., 100.),
            kind: "reaction".into(),
            above: String::new(),
            below: String::new(),
        }];
        assert_eq!(analyze_step(&m, &arrow, 1).unwrap().balanced, Some(true));
        let report = analyze_step(&m, &arrow, 1).unwrap();
        assert_eq!(report.components.len(), 2);
        assert!(report
            .components
            .iter()
            .all(|row| row.role_source == "inferred"));
        assert!(report_text(&report).contains(
            "Reaction component table (role | coefficient | formula | charge | role source):"
        ));
        assert!(
            !report_text(&analyze_step(&m, &arrow, 1).unwrap()).contains("Reaction conditions:")
        );
        m.atoms[1].element = "O".into();
        assert_eq!(analyze_step(&m, &arrow, 1).unwrap().balanced, Some(false));
    }
    #[test]
    fn requires_a_reaction_arrow() {
        assert!(analyze_all(&Molecule::default(), &[]).is_err());
    }
    #[test]
    fn detects_broken_mapped_bond() {
        use crate::{Bond, BondStyle};
        let mut m = Molecule::default();
        for (x, map) in [(100., 1), (120., 2), (400., 1), (430., 2)] {
            let mut a = Atom::new(vec2(x, 100.), "C");
            a.atom_map = Some(map);
            m.atoms.push(a);
        }
        m.bonds.push(Bond {
            a: 0,
            b: 1,
            order: 1,
            style: BondStyle::Normal,
            stereo: None,
            query: None,
        });
        let arrow = vec![Graphic::Arrow {
            start: vec2(200., 100.),
            end: vec2(300., 100.),
            kind: "reaction".into(),
            above: String::new(),
            below: String::new(),
        }];
        let r = analyze_step(&m, &arrow, 1).unwrap();
        assert!(r.bond_changes.iter().any(|s| s.contains("broken")));
    }

    #[test]
    fn explicit_roles_coefficients_and_steps_drive_balance_analysis() {
        let mut m = Molecule::default();
        for (x, role, step, coefficient) in [
            (80., "reactant", 1, 2),
            (420., "product", 1, 1),
            (500., "reactant", 2, 1),
            (760., "product", 2, 1),
        ] {
            let mut a = Atom::new(vec2(x, 120.), "C");
            a.reaction_role = Some(role.into());
            a.reaction_step = step;
            a.coefficient = coefficient;
            m.atoms.push(a);
        }
        let arrows = vec![
            Graphic::Arrow {
                start: vec2(220., 120.),
                end: vec2(300., 120.),
                kind: "reaction".into(),
                above: "heat, reflux".into(),
                below: "ethanol".into(),
            },
            Graphic::Arrow {
                start: vec2(610., 120.),
                end: vec2(680., 120.),
                kind: "reaction".into(),
                above: "hv".into(),
                below: "water".into(),
            },
        ];
        let reports = analyze_all(&m, &arrows).unwrap();
        assert_eq!(reports.len(), 2);
        assert_eq!(reports[0].step, 1);
        assert_eq!(reports[0].balanced, Some(false));
        assert_eq!(reports[0].conditions_above, "heat, reflux");
        assert_eq!(reports[0].conditions_below, "ethanol");
        assert!(report_text(&reports[0]).contains("Above arrow: heat, reflux"));
        assert!(reports[0].reactants[0].starts_with("2 ×"));
        assert_eq!(reports[0].components.len(), 2);
        assert_eq!(reports[0].components[0].role, "reactant");
        assert_eq!(reports[0].components[0].role_source, "explicit");
        assert_eq!(reports[0].components[0].coefficient, 2);
        assert_eq!(reports[0].components[0].formula, "CH4");
        assert_eq!(reports[0].components[0].formal_charge, 0);
        assert!(report_text(&reports[0]).contains("reactant | 2 | CH4 | +0 | explicit"));
        assert_eq!(reports[1].step, 2);
        assert_eq!(reports[1].balanced, Some(true));
        assert_eq!(reports[1].conditions_above, "hv");
        assert_eq!(reports[1].conditions_below, "water");
        assert!(report_all_text(&reports).contains("Above arrow: hv"));
    }

    #[test]
    fn reaction_balance_preserves_isotope_counts_and_solver_enforces_them() {
        let mut m = Molecule::default();
        let mut labeled = Atom::new(vec2(80., 120.), "C");
        labeled.isotope = Some(13);
        labeled.reaction_role = Some("reactant".into());
        let mut natural = Atom::new(vec2(420., 120.), "C");
        natural.reaction_role = Some("product".into());
        m.atoms.extend([labeled, natural]);
        let arrow = vec![Graphic::Arrow {
            start: vec2(220., 120.),
            end: vec2(300., 120.),
            kind: "reaction".into(),
            above: String::new(),
            below: String::new(),
        }];

        let report = analyze_step(&m, &arrow, 1).unwrap();
        assert_eq!(
            report.element_balance,
            vec![("C".into(), 1, 1), ("H".into(), 4, 4)]
        );
        assert_eq!(report.isotope_balance, vec![("13C".into(), 1, 0)]);
        assert_eq!(report.balanced, Some(false));
        assert!(report_text(&report).contains("13C: 1 | 0 | -1"));
        assert!(propose_balancing_coefficients(&m, &arrow, 1).is_err());

        m.atoms[1].isotope = Some(13);
        let conserved = analyze_step(&m, &arrow, 1).unwrap();
        assert_eq!(conserved.isotope_balance, vec![("13C".into(), 1, 1)]);
        assert_eq!(conserved.balanced, Some(true));
    }

    #[test]
    fn proposes_smallest_integer_coefficients_for_water_formation() {
        use crate::{Bond, BondStyle};
        let mut m = Molecule::default();
        let mut add = |element: &str, x: f32, role: &str| {
            let mut atom = Atom::new(vec2(x, 100.), element);
            atom.reaction_role = Some(role.into());
            m.atoms.push(atom);
        };
        add("H", 80., "reactant");
        add("H", 90., "reactant");
        add("O", 120., "reactant");
        add("O", 130., "reactant");
        add("O", 420., "product");
        add("H", 410., "product");
        add("H", 430., "product");
        for (a, b, order) in [(0, 1, 1), (2, 3, 2), (4, 5, 1), (4, 6, 1)] {
            m.bonds.push(Bond {
                a,
                b,
                order,
                style: BondStyle::Normal,
                stereo: None,
                query: None,
            });
        }
        let arrow = vec![Graphic::Arrow {
            start: vec2(220., 100.),
            end: vec2(300., 100.),
            kind: "reaction".into(),
            above: String::new(),
            below: String::new(),
        }];
        let proposal = propose_balancing_coefficients(&m, &arrow, 1).unwrap();
        let by_formula: BTreeMap<_, _> = proposal
            .iter()
            .map(|item| (item.formula.as_str(), item.coefficient))
            .collect();
        assert_eq!(by_formula["H2"], 2);
        assert_eq!(by_formula["O2"], 1);
        assert_eq!(by_formula["H2O"], 2);
    }

    #[test]
    fn refuses_underdetermined_duplicate_component_equation() {
        use crate::{Bond, BondStyle};
        let mut m = Molecule::default();
        for (x, role) in [
            (80., "reactant"),
            (90., "reactant"),
            (420., "product"),
            (430., "product"),
            (470., "product"),
            (480., "product"),
        ] {
            let mut atom = Atom::new(vec2(x, 100.), "H");
            atom.reaction_role = Some(role.into());
            m.atoms.push(atom);
        }
        for (a, b) in [(0, 1), (2, 3), (4, 5)] {
            m.bonds.push(Bond {
                a,
                b,
                order: 1,
                style: BondStyle::Normal,
                stereo: None,
                query: None,
            });
        }
        let arrow = vec![Graphic::Arrow {
            start: vec2(220., 100.),
            end: vec2(300., 100.),
            kind: "reaction".into(),
            above: String::new(),
            below: String::new(),
        }];
        assert!(propose_balancing_coefficients(&m, &arrow, 1)
            .unwrap_err()
            .contains("free dimensions"));
    }

    #[test]
    fn balance_solver_includes_formal_charge_as_a_conservation_row() {
        use crate::{Bond, BondStyle};
        let mut m = Molecule::default();
        for (x, role) in [
            (80., "reactant"),
            (90., "reactant"),
            (420., "product"),
            (430., "product"),
        ] {
            let mut atom = Atom::new(vec2(x, 100.), "H");
            atom.reaction_role = Some(role.into());
            m.atoms.push(atom);
        }
        m.atoms[0].charge = 1;
        for (a, b) in [(0, 1), (2, 3)] {
            m.bonds.push(Bond {
                a,
                b,
                order: 1,
                style: BondStyle::Normal,
                stereo: None,
                query: None,
            });
        }
        let arrow = vec![Graphic::Arrow {
            start: vec2(220., 100.),
            end: vec2(300., 100.),
            kind: "reaction".into(),
            above: String::new(),
            below: String::new(),
        }];
        assert!(propose_balancing_coefficients(&m, &arrow, 1).is_err());
    }
}
