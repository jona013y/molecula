//! Bounded SMARTS import and export for the query features implemented locally.
//! This is a documented subset, not a full SMARTS implementation.

use crate::{chemistry, shortcuts, Atom, Bond, BondStyle, Molecule};
use std::collections::{BTreeSet, HashMap};

fn bond_token(bond: &Bond) -> Result<&'static str, String> {
    if bond.stereo.is_some() {
        return Err("SMARTS export of stereochemical bond predicates is not supported yet.".into());
    }
    if let Some(query) = bond.query.as_deref() {
        if query.split(';').any(|term| term.starts_with("topology-")) {
            return Err("SMARTS export cannot preserve Molecula's ring-only/chain-only bond topology predicate.".into());
        }
        let mut tokens = BTreeSet::new();
        for term in query.split(';').filter(|term| !term.is_empty()) {
            tokens.insert(match term {
                "any" => "~",
                "single-or-double" => "-,=",
                "single-or-aromatic" => "-,:",
                "double-or-aromatic" => "=,:",
                _ => {
                    return Err(format!(
                        "SMARTS export does not support bond query `{term}`."
                    ))
                }
            });
        }
        if tokens.len() == 1 {
            return Ok(*tokens.first().unwrap());
        }
        return Err(
            "SMARTS export cannot combine multiple independent bond query alternatives.".into(),
        );
    }
    match bond.style {
        BondStyle::Normal => match bond.order {
            1 => Ok("-"),
            2 => Ok("="),
            3 => Ok("#"),
            4 => Ok(":"),
            _ => Err(format!(
                "SMARTS export does not support bond order {}.",
                bond.order
            )),
        },
        BondStyle::Aromatic if bond.order == 4 => Ok(":"),
        BondStyle::Any => Ok("~"),
        BondStyle::Aromatic => Err("Aromatic bond data has an inconsistent order.".into()),
        BondStyle::Wedge | BondStyle::Hash | BondStyle::Wavy => {
            Err("SMARTS export cannot preserve wedge/hash/wavy stereobond display marks.".into())
        }
        BondStyle::Coordinate => Err("SMARTS export does not represent coordination bonds.".into()),
    }
}

fn atomic_number(element: &str) -> Result<usize, String> {
    chemistry::SYMBOLS
        .iter()
        .position(|candidate| *candidate == element)
        .map(|index| index + 1)
        .ok_or_else(|| format!("Unknown element symbol `{element}`."))
}

fn element_symbol(element: &str, aromatic: bool) -> Result<String, String> {
    if aromatic {
        return Ok(match element {
            "B" => "b",
            "C" => "c",
            "N" => "n",
            "O" => "o",
            "P" => "p",
            "S" => "s",
            _ => {
                return Err(format!(
                    "SMARTS aromatic atoms of element {element} are unsupported."
                ))
            }
        }
        .into());
    }
    Ok(element.to_string())
}

fn atom_token(atom: &crate::Atom) -> Result<String, String> {
    if atom.alias.is_some() || atom.r_group.is_some() {
        return Err("SMARTS export cannot preserve aliases or R-group labels.".into());
    }
    if atom.stereo.is_some() || atom.radical != 0 {
        return Err(
            "SMARTS export cannot preserve atom stereochemistry or radical annotations.".into(),
        );
    }
    if atom.implicit_h_override.is_some() || atom.valence_override.is_some() || atom.explicit_h != 0
    {
        return Err(
            "SMARTS export cannot preserve editor-specific hydrogen or valence overrides.".into(),
        );
    }

    let query_terms: Vec<&str> = atom
        .query
        .as_deref()
        .unwrap_or("")
        .split(';')
        .filter(|term| !term.is_empty())
        .collect();
    let mut predicates = Vec::new();
    let mut element_selector = None;
    let mut aromatic_query = atom.aromatic;
    let mut aliphatic_query = false;
    let wildcard_element_test = query_terms.iter().any(|term| {
        matches!(
            shortcuts::parse_query_predicate(term),
            Some(
                shortcuts::QueryPredicate::AnyAtom
                    | shortcuts::QueryPredicate::Heteroatom
                    | shortcuts::QueryPredicate::ElementList(_)
                    | shortcuts::QueryPredicate::ExcludedElementList(_)
            )
        )
    });
    for term in query_terms {
        let predicate = shortcuts::parse_query_predicate(term)
            .ok_or_else(|| format!("SMARTS export does not support atom query `{term}`."))?;
        let text = match predicate {
            shortcuts::QueryPredicate::AnyAtom => "*".into(),
            shortcuts::QueryPredicate::Heteroatom => "!#6&!#1".into(),
            shortcuts::QueryPredicate::ElementList(elements) => {
                let numbers = elements
                    .iter()
                    .map(|element| atomic_number(element).map(|n| format!("#{n}")))
                    .collect::<Result<Vec<_>, _>>()?;
                if element_selector.replace(numbers.join(",")).is_some() {
                    return Err(
                        "SMARTS export refuses conflicting element selectors on one atom.".into(),
                    );
                }
                continue;
            }
            shortcuts::QueryPredicate::ExcludedElementList(elements) => {
                let numbers = elements
                    .iter()
                    .map(|element| atomic_number(element).map(|n| format!("!#{n}")))
                    .collect::<Result<Vec<_>, _>>()?;
                if element_selector.replace(numbers.join("&")).is_some() {
                    return Err(
                        "SMARTS export refuses conflicting element selectors on one atom.".into(),
                    );
                }
                continue;
            }
            shortcuts::QueryPredicate::Aromatic => {
                aromatic_query = true;
                continue;
            }
            shortcuts::QueryPredicate::Aliphatic => {
                aliphatic_query = true;
                "A".into()
            }
            shortcuts::QueryPredicate::HydrogenCount(count) => format!("H{count}"),
            shortcuts::QueryPredicate::Connectivity(count) => format!("X{count}"),
            shortcuts::QueryPredicate::Degree(count) => format!("D{count}"),
            shortcuts::QueryPredicate::RingCount(Some(0)) => "!R".into(),
            shortcuts::QueryPredicate::RingCount(Some(_)) => {
                return Err("SMARTS ring-count semantics can differ from Molecula's local cycle-basis count; only `.R0` is exported.".into());
            }
            shortcuts::QueryPredicate::MinimumRingSize(size) => format!("r{size}"),
            shortcuts::QueryPredicate::Valence(count) => format!("v{count}"),
            shortcuts::QueryPredicate::Unsaturated
            | shortcuts::QueryPredicate::ImplicitHydrogenCount(_)
            | shortcuts::QueryPredicate::RingCount(None)
            | shortcuts::QueryPredicate::RingBonds(_)
            | shortcuts::QueryPredicate::Substitutions(_) => {
                return Err(format!(
                    "SMARTS export cannot preserve Molecula atom predicate `{term}` exactly."
                ));
            }
        };
        predicates.push(text);
    }
    if aromatic_query && aliphatic_query {
        return Err(
            "An atom cannot be both aromatic and aliphatic in the supported SMARTS subset.".into(),
        );
    }
    let has_element_selector = element_selector.is_some();
    let symbol = if has_element_selector || wildcard_element_test {
        None
    } else {
        Some(element_symbol(&atom.element, aromatic_query)?)
    };
    let mut base_terms = Vec::new();
    if let Some(selector) = element_selector.take() {
        base_terms.push(selector);
    } else if let Some(symbol) = symbol {
        base_terms.push(symbol);
    }
    base_terms.extend(predicates);
    if aromatic_query && (has_element_selector || wildcard_element_test) {
        base_terms.push("a".into());
    }

    let has_special_attributes =
        atom.isotope.is_some() || atom.charge != 0 || atom.atom_map.is_some();
    let ordinary_unbracketed = atom.query.is_none()
        && !has_special_attributes
        && !atom.aromatic
        && matches!(
            atom.element.as_str(),
            "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
        );
    if ordinary_unbracketed {
        return element_symbol(&atom.element, false);
    }

    if let Some(isotope) = atom.isotope {
        if atom.query.is_some() || base_terms.is_empty() {
            return Err(
                "SMARTS export of isotopically constrained query atoms is not supported.".into(),
            );
        }
        base_terms[0] = format!("{isotope}{}", base_terms[0]);
    }
    if atom.charge != 0 {
        if wildcard_element_test {
            return Err("SMARTS export cannot preserve a formal-charge constraint combined with a wildcard/list query atom in Molecula's matcher.".into());
        }
        let charge = if atom.charge.unsigned_abs() == 1 {
            if atom.charge > 0 { "+" } else { "-" }.to_string()
        } else {
            format!(
                "{}{magnitude}",
                if atom.charge > 0 { "+" } else { "-" },
                magnitude = atom.charge.unsigned_abs()
            )
        };
        if atom.query.is_none() {
            base_terms[0].push_str(&charge);
        } else {
            base_terms.push(charge);
        }
    }
    if base_terms.is_empty() {
        return Err("Could not encode an empty SMARTS atom expression.".into());
    }
    let map = atom
        .atom_map
        .map(|number| format!(":{number}"))
        .unwrap_or_default();
    Ok(format!("[{}{}]", base_terms.join(";"), map))
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct BondSpec {
    order: u8,
    style: BondStyle,
    query: Option<String>,
}

impl BondSpec {
    fn from_smarts(token: &str) -> Result<Self, String> {
        let (order, style, query) = match token {
            "-" => (1, BondStyle::Normal, None),
            "=" => (2, BondStyle::Normal, None),
            "#" => (3, BondStyle::Normal, None),
            ":" => (4, BondStyle::Aromatic, None),
            "~" => (1, BondStyle::Any, None),
            "-,=" => (1, BondStyle::Normal, Some("single-or-double")),
            "-,:" => (1, BondStyle::Normal, Some("single-or-aromatic")),
            "=,:" => (2, BondStyle::Normal, Some("double-or-aromatic")),
            _ => return Err(format!("Unsupported SMARTS bond expression `{token}`.")),
        };
        Ok(Self {
            order,
            style,
            query: query.map(str::to_string),
        })
    }

    fn default_smarts() -> Self {
        Self {
            order: 1,
            style: BondStyle::Normal,
            query: Some("single-or-aromatic".into()),
        }
    }

    fn into_bond(self, a: usize, b: usize) -> Bond {
        Bond {
            a,
            b,
            order: self.order,
            style: self.style,
            stereo: None,
            query: self.query,
        }
    }
}

fn element_from_atomic_number(number: usize) -> Result<&'static str, String> {
    chemistry::SYMBOLS
        .get(
            number
                .checked_sub(1)
                .ok_or("Atomic number must be positive")?,
        )
        .copied()
        .ok_or_else(|| format!("SMARTS atomic number #{number} is outside the periodic table."))
}

fn parse_charge_term(term: &str) -> Option<i8> {
    let sign = term.chars().next()?;
    if sign != '+' && sign != '-' {
        return None;
    }
    let magnitude = if term.len() == 1 {
        1
    } else {
        let parsed: u8 = term[1..].parse().ok()?;
        if parsed == 0 {
            return None;
        }
        parsed
    };
    let signed = if sign == '+' {
        i16::from(magnitude)
    } else {
        -i16::from(magnitude)
    };
    i8::try_from(signed).ok()
}

fn parse_atomic_selector(term: &str) -> Result<Option<(String, String)>, String> {
    if term == "*" {
        return Ok(Some(("C".into(), "*".into())));
    }
    if term == "!#6&!#1" {
        return Ok(Some(("C".into(), "Q".into())));
    }
    if term.starts_with('#') {
        let elements = term
            .split(',')
            .map(|number| {
                let number = number
                    .strip_prefix('#')
                    .ok_or("Malformed SMARTS atomic-number list")?
                    .parse::<usize>()
                    .map_err(|_| "Malformed SMARTS atomic number")?;
                element_from_atomic_number(number).map(str::to_string)
            })
            .collect::<Result<Vec<_>, _>>()?;
        if elements.is_empty() {
            return Err("Empty SMARTS element list.".into());
        }
        return Ok(Some(("C".into(), elements.join(","))));
    }
    if term.starts_with("!#") {
        let mut elements = Vec::new();
        for field in term.split('&') {
            let number = field
                .strip_prefix("!#")
                .ok_or("Malformed excluded SMARTS element list")?
                .parse::<usize>()
                .map_err(|_| "Malformed excluded SMARTS atomic number")?;
            elements.push(element_from_atomic_number(number)?.to_string());
        }
        if elements.as_slice() == ["C", "H"] || elements.as_slice() == ["H", "C"] {
            return Ok(Some(("C".into(), "Q".into())));
        }
        return Ok(Some(("C".into(), format!("!{}", elements.join(",")))));
    }
    Ok(None)
}

fn parse_bracket_atom(content: &str, position: macroquad::prelude::Vec2) -> Result<Atom, String> {
    if content.is_empty() {
        return Err("Empty SMARTS bracket atom.".into());
    }
    let (expression, map) = if let Some((left, right)) = content.rsplit_once(':') {
        if right.is_empty() || !right.chars().all(|ch| ch.is_ascii_digit()) {
            return Err("Malformed SMARTS atom-map suffix.".into());
        }
        let map: u16 = right
            .parse()
            .map_err(|_| "SMARTS atom-map number is outside the supported range")?;
        (left, Some(map))
    } else {
        (content, None)
    };
    let terms: Vec<&str> = expression.split(';').collect();
    if terms.iter().any(|term| term.is_empty()) {
        return Err("Empty SMARTS atom predicate.".into());
    }
    let mut first = terms[0];
    let isotope_digits = first.chars().take_while(|ch| ch.is_ascii_digit()).count();
    let isotope = if isotope_digits > 0
        && first
            .as_bytes()
            .get(isotope_digits)
            .is_some_and(u8::is_ascii_alphabetic)
    {
        let value: u16 = first[..isotope_digits]
            .parse()
            .map_err(|_| "SMARTS isotope mass is outside the supported range")?;
        first = &first[isotope_digits..];
        Some(value)
    } else {
        None
    };

    let mut charge = 0i8;
    if let Some((symbol, suffix)) = first.rsplit_once(['+', '-']) {
        if !symbol.is_empty() && !symbol.starts_with('#') && !symbol.starts_with('!') {
            let sign = first[symbol.len()..symbol.len() + 1].to_string();
            let suffix = format!("{sign}{suffix}");
            if let Some(parsed) = parse_charge_term(&suffix) {
                first = symbol;
                charge = parsed;
            }
        }
    }
    if let Some(parsed) = parse_charge_term(first) {
        first = "*";
        charge = parsed;
    }

    let (element, mut query) = if first == "a" {
        ("C".to_string(), Some("*;.a".to_string()))
    } else if first == "A" {
        ("C".to_string(), Some("*;.A".to_string()))
    } else if let Some((element, selector)) = parse_atomic_selector(first)? {
        (element, Some(selector))
    } else {
        let aromatic = first.chars().next().is_some_and(char::is_lowercase);
        let symbol = if aromatic {
            match first {
                "b" => "B",
                "c" => "C",
                "n" => "N",
                "o" => "O",
                "p" => "P",
                "s" => "S",
                _ => return Err(format!("Unsupported aromatic SMARTS atom `{first}`.")),
            }
        } else {
            first
        };
        if !chemistry::valid_element(symbol) {
            return Err(format!("Unsupported SMARTS element `{symbol}`."));
        }
        (
            symbol.to_string(),
            if aromatic { Some(".a".into()) } else { None },
        )
    };

    let mut extra = Vec::new();
    for term in terms.iter().skip(1).copied() {
        if let Some(parsed) = parse_charge_term(term) {
            charge = parsed;
            continue;
        }
        let predicate = match term {
            "a" => ".a".to_string(),
            "A" => ".A".to_string(),
            "!R" => ".R0".to_string(),
            _ if term.starts_with('H') && term[1..].chars().all(|ch| ch.is_ascii_digit()) => {
                format!(".H{}", &term[1..])
            }
            _ if term.starts_with('X') && term[1..].chars().all(|ch| ch.is_ascii_digit()) => {
                format!(".X{}", &term[1..])
            }
            _ if term.starts_with('D') && term[1..].chars().all(|ch| ch.is_ascii_digit()) => {
                format!(".D{}", &term[1..])
            }
            _ if term.starts_with('r') && term[1..].chars().all(|ch| ch.is_ascii_digit()) => {
                format!(".r{}", &term[1..])
            }
            _ if term.starts_with('v') && term[1..].chars().all(|ch| ch.is_ascii_digit()) => {
                format!(".v{}", &term[1..])
            }
            _ if term.starts_with('R') && term[1..].chars().all(|ch| ch == '0') => {
                format!(".R{}", &term[1..])
            }
            _ => return Err(format!("Unsupported SMARTS atom predicate `{term}`.")),
        };
        extra.push(predicate);
    }
    if !extra.is_empty() {
        let existing = query.take().unwrap_or_default();
        let mut predicates: Vec<String> = existing
            .split(';')
            .filter(|part| !part.is_empty())
            .map(str::to_string)
            .collect();
        predicates.extend(extra);
        query = Some(predicates.join(";"));
    }
    if query
        .as_deref()
        .is_some_and(|value| value.contains("*") || value.contains(',') || value == "Q")
        && charge != 0
    {
        return Err("Molecula cannot preserve formal-charge constraints combined with a wildcard/list SMARTS atom.".into());
    }
    let aromatic = query
        .as_deref()
        .is_some_and(|value| value.split(';').any(|term| term == ".a"));
    let mut atom = Atom::new(position, &element);
    if let Some(terms) = query {
        let mut merged: Option<String> = None;
        for term in terms.split(';').filter(|term| !term.is_empty()) {
            merged = Some(shortcuts::merge_query_terms(merged.as_deref(), term)?);
        }
        atom.query = merged;
    }
    atom.aromatic = aromatic && atom.query.is_none();
    atom.isotope = isotope;
    atom.charge = charge;
    atom.atom_map = map;
    Ok(atom)
}

/// Parse the supported SMARTS subset used by the editor and exporter.
/// General recursive SMARTS, stereochemistry, and unsupported predicates are
/// rejected with an error; this never falls back to a guessed structure.
pub fn parse(input: &str) -> Result<Molecule, String> {
    let input = input.trim();
    if input.is_empty() {
        return Err("The SMARTS input is empty.".into());
    }
    if input.len() > 100_000 {
        return Err("SMARTS input is limited to 100,000 bytes.".into());
    }
    let mut molecule = Molecule::default();
    let mut current = None;
    let mut pending_bond: Option<BondSpec> = None;
    let mut branches = Vec::new();
    let mut rings = HashMap::<String, (usize, Option<BondSpec>)>::new();
    let mut last_was_dot = false;
    let bytes = input.as_bytes();
    let mut cursor = 0usize;
    while cursor < bytes.len() {
        let ch = bytes[cursor] as char;
        if ch.is_ascii_whitespace() {
            return Err("SMARTS whitespace is only allowed around the complete pattern.".into());
        }
        if ch == '[' {
            let close = input[cursor + 1..]
                .find(']')
                .map(|offset| cursor + 1 + offset)
                .ok_or("Unclosed SMARTS bracket atom")?;
            let atom =
                parse_bracket_atom(&input[cursor + 1..close], macroquad::prelude::vec2(0., 0.))?;
            let index = molecule.atoms.len();
            molecule.atoms.push(atom);
            if molecule.atoms.len() > 512 {
                return Err("SMARTS input is limited to 512 atoms.".into());
            }
            last_was_dot = false;
            if let Some(previous) = current {
                let spec = pending_bond.take().unwrap_or_else(BondSpec::default_smarts);
                molecule.bonds.push(spec.into_bond(previous, index));
            } else {
                if pending_bond.is_some() {
                    return Err("A SMARTS bond expression cannot begin a component.".into());
                }
            }
            current = Some(index);
            cursor = close + 1;
            continue;
        }
        if ch.is_ascii_uppercase() || ch.is_ascii_lowercase() {
            let end = if cursor + 1 < bytes.len()
                && bytes[cursor + 1].is_ascii_lowercase()
                && ch.is_ascii_uppercase()
            {
                cursor + 2
            } else {
                cursor + 1
            };
            let symbol = &input[cursor..end];
            let (element, aromatic) = if ch.is_ascii_lowercase() {
                let element = match symbol {
                    "b" => "B",
                    "c" => "C",
                    "n" => "N",
                    "o" => "O",
                    "p" => "P",
                    "s" => "S",
                    _ => return Err(format!("Unsupported aromatic SMARTS atom `{symbol}`.")),
                };
                (element, true)
            } else {
                (symbol, false)
            };
            if !chemistry::valid_element(element) {
                return Err(format!("Unsupported SMARTS element `{element}`."));
            }
            let mut atom = Atom::new(macroquad::prelude::vec2(0., 0.), element);
            atom.aromatic = aromatic;
            let index = molecule.atoms.len();
            molecule.atoms.push(atom);
            if molecule.atoms.len() > 512 {
                return Err("SMARTS input is limited to 512 atoms.".into());
            }
            last_was_dot = false;
            if let Some(previous) = current {
                let spec = pending_bond.take().unwrap_or_else(BondSpec::default_smarts);
                molecule.bonds.push(spec.into_bond(previous, index));
            } else {
                if pending_bond.is_some() {
                    return Err("A SMARTS bond expression cannot begin a component.".into());
                }
            }
            current = Some(index);
            cursor = end;
            continue;
        }
        if matches!(ch, '-' | '=' | '#' | ':' | '~' | ',') {
            if current.is_none() || pending_bond.is_some() {
                return Err("Misplaced or repeated SMARTS bond expression.".into());
            }
            let start = cursor;
            cursor += 1;
            while cursor < bytes.len()
                && matches!(bytes[cursor] as char, '-' | '=' | '#' | ':' | '~' | ',')
            {
                cursor += 1;
            }
            pending_bond = Some(BondSpec::from_smarts(&input[start..cursor])?);
            continue;
        }
        if ch.is_ascii_digit() || ch == '%' {
            let label = if ch == '%' {
                if cursor + 2 >= bytes.len()
                    || !bytes[cursor + 1..cursor + 3].iter().all(u8::is_ascii_digit)
                {
                    return Err("SMARTS percent ring labels must contain two digits.".into());
                }
                let label = input[cursor + 1..cursor + 3].to_string();
                cursor += 3;
                label
            } else {
                cursor += 1;
                ch.to_string()
            };
            let atom = current.ok_or("SMARTS ring label has no current atom")?;
            if let Some((other, opened_spec)) = rings.remove(&label) {
                if other == atom {
                    return Err("A SMARTS ring closure cannot connect an atom to itself.".into());
                }
                if let (Some(first), Some(second)) = (&opened_spec, &pending_bond) {
                    if first != second {
                        return Err(
                            "The two ends of a SMARTS ring closure specify different bonds.".into(),
                        );
                    }
                }
                let spec = pending_bond
                    .take()
                    .or(opened_spec)
                    .unwrap_or_else(BondSpec::default_smarts);
                molecule.bonds.push(spec.into_bond(other, atom));
            } else {
                rings.insert(label, (atom, pending_bond.take()));
            }
            continue;
        }
        match ch {
            '(' => {
                if pending_bond.is_some() {
                    return Err("A SMARTS branch must begin with an atom.".into());
                }
                branches.push(current.ok_or("SMARTS branch has no parent atom")?);
                cursor += 1;
            }
            ')' => {
                if pending_bond.is_some() {
                    return Err("A SMARTS branch ends with an incomplete bond.".into());
                }
                current = Some(
                    branches
                        .pop()
                        .ok_or("Unmatched SMARTS closing parenthesis")?,
                );
                cursor += 1;
            }
            '.' => {
                if pending_bond.is_some()
                    || !branches.is_empty()
                    || current.is_none()
                    || last_was_dot
                {
                    return Err(
                        "A SMARTS component separator cannot occur inside a branch or bond.".into(),
                    );
                }
                current = None;
                last_was_dot = true;
                cursor += 1;
            }
            _ => {
                return Err(format!(
                    "Unsupported SMARTS token `{ch}` at byte {}.",
                    cursor + 1
                ))
            }
        }
    }
    if !branches.is_empty() {
        return Err("Unclosed SMARTS branch.".into());
    }
    if !rings.is_empty() {
        return Err("A SMARTS ring closure label appears only once.".into());
    }
    if pending_bond.is_some() {
        return Err("SMARTS input ends with an incomplete bond.".into());
    }
    if last_was_dot {
        return Err("SMARTS input ends with an empty component.".into());
    }
    if molecule.atoms.len() > 512 {
        return Err("SMARTS input is limited to 512 atoms.".into());
    }
    let mut unique_bonds = BTreeSet::new();
    for bond in &molecule.bonds {
        let edge = (bond.a.min(bond.b), bond.a.max(bond.b));
        if bond.a == bond.b || !unique_bonds.insert(edge) {
            return Err("SMARTS input contains a self-loop or duplicate bond.".into());
        }
    }
    layout_imported_molecule(&mut molecule);
    Ok(molecule)
}

/// Gives imported query graphs a readable, deterministic 2D arrangement.
/// Acyclic fragments start as bond-length trees, simple rings as polygons, then
/// a bounded spring relaxation resolves fused closures and atom overlaps.
pub(crate) fn layout_imported_molecule(molecule: &mut Molecule) {
    use macroquad::prelude::vec2;

    let count = molecule.atoms.len();
    if count == 0 {
        return;
    }
    let mut adjacency = vec![Vec::new(); count];
    for bond in &molecule.bonds {
        adjacency[bond.a].push(bond.b);
        adjacency[bond.b].push(bond.a);
    }
    for neighbors in &mut adjacency {
        neighbors.sort_unstable();
    }

    let mut seen = vec![false; count];
    let mut placed = vec![vec2(0., 0.); count];
    let mut next_component_x = 0.0f32;
    for root in 0..count {
        if seen[root] {
            continue;
        }
        let mut members = Vec::new();
        let mut queue = vec![root];
        seen[root] = true;
        while let Some(atom) = queue.pop() {
            members.push(atom);
            for &neighbor in &adjacency[atom] {
                if !seen[neighbor] {
                    seen[neighbor] = true;
                    queue.push(neighbor);
                }
            }
        }
        members.sort_unstable();
        let local_index: HashMap<usize, usize> = members
            .iter()
            .enumerate()
            .map(|(local, &global)| (global, local))
            .collect();
        let edges: Vec<(usize, usize)> = molecule
            .bonds
            .iter()
            .filter_map(|bond| Some((*local_index.get(&bond.a)?, *local_index.get(&bond.b)?)))
            .collect();
        let mut positions = vec![vec2(0., 0.); members.len()];

        let is_simple_ring = members.len() >= 3
            && edges.len() == members.len()
            && members.iter().all(|&atom| adjacency[atom].len() == 2);
        if is_simple_ring {
            let mut order = vec![root];
            let mut previous = usize::MAX;
            let mut current = root;
            loop {
                let next = adjacency[current]
                    .iter()
                    .copied()
                    .find(|neighbor| *neighbor != previous)
                    .unwrap_or(root);
                if next == root {
                    break;
                }
                if order.contains(&next) {
                    break;
                }
                order.push(next);
                previous = current;
                current = next;
            }
            if order.len() == members.len() {
                let radius = crate::BOND
                    / (std::f32::consts::PI / members.len() as f32)
                        .sin()
                        .max(0.01)
                    / 2.0;
                for (index, &atom) in order.iter().enumerate() {
                    let angle = -std::f32::consts::FRAC_PI_2
                        + std::f32::consts::TAU * index as f32 / members.len() as f32;
                    positions[local_index[&atom]] =
                        vec2(angle.cos() * radius, angle.sin() * radius);
                }
            }
        } else {
            let mut tree_seen = vec![false; members.len()];
            let root_local = local_index[&root];
            tree_seen[root_local] = true;
            let mut stack = vec![(root, None, 0.0f32)];
            while let Some((atom, parent, heading)) = stack.pop() {
                let mut children: Vec<usize> = adjacency[atom]
                    .iter()
                    .copied()
                    .filter(|neighbor| Some(*neighbor) != parent)
                    .filter(|neighbor| !tree_seen[local_index[neighbor]])
                    .collect();
                children.sort_unstable();
                let child_count = children.len();
                for (ordinal, child) in children.into_iter().enumerate() {
                    let direction = if parent.is_none() {
                        std::f32::consts::TAU * ordinal as f32 / child_count.max(1) as f32
                    } else if child_count == 1 {
                        heading
                    } else {
                        let step = (std::f32::consts::TAU / child_count as f32)
                            .min(2.0 * std::f32::consts::PI / 3.0);
                        heading + step * ordinal as f32
                            - step * (child_count.saturating_sub(1) as f32 / 2.0)
                    };
                    let child_local = local_index[&child];
                    positions[child_local] = positions[local_index[&atom]]
                        + vec2(direction.cos(), direction.sin()) * crate::BOND;
                    tree_seen[child_local] = true;
                    stack.push((child, Some(atom), direction));
                }
            }
        }

        // Only cyclic polygraphs need relaxation. Tree layouts and simple
        // regular rings already have exact bond lengths and clearer angles.
        if edges.len() >= members.len() && !is_simple_ring {
            for step in 0..100 {
                let mut forces = vec![vec2(0., 0.); positions.len()];
                for i in 0..positions.len() {
                    for j in i + 1..positions.len() {
                        let mut vector = positions[i] - positions[j];
                        let mut distance = vector.length();
                        if distance < 0.01 {
                            vector = vec2((i + 1) as f32, (j + 1) as f32);
                            distance = vector.length();
                        }
                        let force = vector / distance * (crate::BOND * crate::BOND / distance);
                        forces[i] += force;
                        forces[j] -= force;
                    }
                }
                for &(a, b) in &edges {
                    let vector = positions[a] - positions[b];
                    let distance = vector.length().max(0.01);
                    let force = vector / distance * (distance * distance / crate::BOND);
                    forces[a] -= force;
                    forces[b] += force;
                }
                let temperature = 8.0 * (1.0 - step as f32 / 100.0);
                for (position, force) in positions.iter_mut().zip(forces) {
                    let magnitude = force.length();
                    if magnitude > 0.001 {
                        *position += force / magnitude * magnitude.min(temperature);
                    }
                }
            }
        }

        let min_x = positions
            .iter()
            .map(|position| position.x)
            .fold(f32::INFINITY, f32::min);
        let max_x = positions
            .iter()
            .map(|position| position.x)
            .fold(f32::NEG_INFINITY, f32::max);
        let shift = next_component_x - min_x;
        for (local, &global) in members.iter().enumerate() {
            placed[global] = positions[local] + vec2(shift, 0.0);
        }
        next_component_x = max_x + shift + 3.0 * crate::BOND;
    }

    let min = placed
        .iter()
        .copied()
        .fold(vec2(f32::INFINITY, f32::INFINITY), |acc, p| {
            vec2(acc.x.min(p.x), acc.y.min(p.y))
        });
    let max = placed
        .iter()
        .copied()
        .fold(vec2(f32::NEG_INFINITY, f32::NEG_INFINITY), |acc, p| {
            vec2(acc.x.max(p.x), acc.y.max(p.y))
        });
    let center = (min + max) * 0.5;
    for (atom, position) in molecule.atoms.iter_mut().zip(placed) {
        atom.pos = position - center;
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
    ring_annotations: &[Vec<(usize, String)>],
    bonds: &[Bond],
) -> Result<String, String> {
    let mut out = atoms[atom].clone();
    for (label, token) in &ring_annotations[atom] {
        out.push_str(token);
        out.push_str(&ring_label(*label));
    }
    let descendants = &children[atom];
    if let Some((&(continuation, continuation_edge), branches)) = descendants.split_last() {
        for &(child, edge) in branches {
            out.push('(');
            out.push_str(bond_token(&bonds[edge])?);
            out.push_str(&render_tree(
                child,
                atoms,
                children,
                ring_annotations,
                bonds,
            )?);
            out.push(')');
        }
        out.push_str(bond_token(&bonds[continuation_edge])?);
        out.push_str(&render_tree(
            continuation,
            atoms,
            children,
            ring_annotations,
            bonds,
        )?);
    }
    Ok(out)
}

/// Serialize the representable SMARTS subset without an external toolkit.
/// Refuses every editor predicate or annotation that cannot be represented
/// faithfully by this serializer.
pub fn export(molecule: &Molecule) -> Result<String, String> {
    if molecule.atoms.is_empty() {
        return Err("Draw a structure before exporting SMARTS.".into());
    }
    if molecule.atoms.len() > 512 {
        return Err("Local SMARTS export is limited to 512 atoms.".into());
    }
    let atom_tokens = molecule
        .atoms
        .iter()
        .map(atom_token)
        .collect::<Result<Vec<_>, _>>()?;
    let mut adjacency = vec![Vec::<(usize, usize)>::new(); molecule.atoms.len()];
    let mut seen_edges = BTreeSet::new();
    for (edge, bond) in molecule.bonds.iter().enumerate() {
        if bond.a >= molecule.atoms.len() || bond.b >= molecule.atoms.len() || bond.a == bond.b {
            return Err("The structure contains an invalid bond reference.".into());
        }
        let _ = bond_token(bond)?;
        let key = (bond.a.min(bond.b), bond.a.max(bond.b));
        if !seen_edges.insert(key) {
            return Err("SMARTS export refuses duplicate bonds.".into());
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
        while let Some((current, next_index)) = stack.pop() {
            if next_index >= adjacency[current].len() {
                continue;
            }
            stack.push((current, next_index + 1));
            let (neighbor, edge) = adjacency[current][next_index];
            if !visited[neighbor] {
                visited[neighbor] = true;
                children[current].push((neighbor, edge));
                parent_edges.insert(edge);
                stack.push((neighbor, 0));
            }
        }
    }
    let ring_edges: Vec<usize> = (0..molecule.bonds.len())
        .filter(|edge| !parent_edges.contains(edge))
        .collect();
    if ring_edges.len() > 99 {
        return Err("SMARTS ring closures are limited to 99 labels in this exporter.".into());
    }
    let mut ring_annotations = vec![Vec::<(usize, String)>::new(); molecule.atoms.len()];
    for (index, &bond_index) in ring_edges.iter().enumerate() {
        let label = index + 1;
        let token = bond_token(&molecule.bonds[bond_index])?.to_string();
        let bond = &molecule.bonds[bond_index];
        ring_annotations[bond.a].push((label, token.clone()));
        ring_annotations[bond.b].push((label, token));
    }
    for annotations in &mut ring_annotations {
        annotations.sort_by_key(|(label, _)| *label);
    }
    let mut components = Vec::new();
    for root in roots {
        components.push(render_tree(
            root,
            &atom_tokens,
            &children,
            &ring_annotations,
            &molecule.bonds,
        )?);
    }
    Ok(components.join("."))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Atom, Bond};
    use macroquad::prelude::vec2;

    fn atom(molecule: &mut Molecule, element: &str) -> usize {
        let index = molecule.atoms.len();
        molecule.atoms.push(Atom::new(vec2(0., 0.), element));
        index
    }

    fn bond(molecule: &mut Molecule, a: usize, b: usize, order: u8) {
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
    }

    #[test]
    fn serializes_linear_branching_and_ring_graphs() {
        let mut chain = Molecule::default();
        let c1 = atom(&mut chain, "C");
        let c2 = atom(&mut chain, "C");
        let oxygen = atom(&mut chain, "O");
        bond(&mut chain, c1, c2, 1);
        bond(&mut chain, c2, oxygen, 1);
        assert_eq!(export(&chain).unwrap(), "C-C-O");

        let mut branch = Molecule::default();
        for _ in 0..4 {
            atom(&mut branch, "C");
        }
        for end in [0, 2, 3] {
            bond(&mut branch, 1, end, 1);
        }
        assert_eq!(export(&branch).unwrap(), "C-C(-C)-C");

        let mut ring = Molecule::default();
        for _ in 0..3 {
            atom(&mut ring, "C");
        }
        for index in 0..3 {
            bond(&mut ring, index, (index + 1) % 3, 1);
        }
        assert_eq!(export(&ring).unwrap(), "C-1-C-C-1");
    }

    #[test]
    fn serializes_supported_atom_and_bond_queries_with_maps() {
        let mut molecule = Molecule::default();
        let mut carbon = Atom::new(vec2(0., 0.), "C");
        carbon.query = Some("*;.H1".into());
        carbon.atom_map = Some(7);
        molecule.atoms.push(carbon);
        molecule.atoms.push(Atom::new(vec2(1., 0.), "O"));
        let bond = Bond {
            a: 0,
            b: 1,
            order: 0,
            style: BondStyle::Normal,
            stereo: None,
            query: Some("single-or-double".into()),
        };
        molecule.bonds.push(bond);
        assert_eq!(export(&molecule).unwrap(), "[*;H1:7]-,=O");

        molecule.atoms[0].query = Some("Au,Ag".into());
        molecule.atoms[0].atom_map = None;
        assert_eq!(export(&molecule).unwrap(), "[#79,#47]-,=O");
    }

    #[test]
    fn serializes_isotopes_formal_charges_and_atom_maps() {
        let mut molecule = Molecule::default();
        let mut nitrogen = Atom::new(vec2(0., 0.), "N");
        nitrogen.isotope = Some(15);
        nitrogen.charge = 1;
        nitrogen.atom_map = Some(4);
        molecule.atoms.push(nitrogen);
        assert_eq!(export(&molecule).unwrap(), "[15N+:4]");
    }

    #[test]
    fn refuses_query_features_whose_semantics_cannot_be_preserved() {
        let mut molecule = Molecule::default();
        let mut carbon = Atom::new(vec2(0., 0.), "C");
        carbon.query = Some(".s2".into());
        molecule.atoms.push(carbon);
        assert!(export(&molecule).unwrap_err().contains("cannot preserve"));
        molecule.atoms[0].query = Some(".R2".into());
        assert!(export(&molecule).unwrap_err().contains("cycle-basis"));
    }

    #[test]
    fn parses_supported_patterns_and_round_trips_exported_subset() {
        for pattern in [
            "C-C-O",
            "C-C(-C)-C",
            "C-1-C-C-1",
            "[*;H1:7]-,=O",
            "[#79,#47]",
            "[15N+:4]",
            "[c]:[c]:[c]",
            "C.O",
        ] {
            let molecule = parse(pattern).unwrap_or_else(|error| panic!("{pattern}: {error}"));
            assert_eq!(export(&molecule).unwrap(), pattern, "{pattern}");
        }
    }

    #[test]
    fn rejects_unsupported_or_malformed_smarts_without_guessing() {
        for pattern in [
            "",
            "C(",
            "C-",
            "C1CC",
            "[C@@H]",
            "[$(C=O)]",
            "C..O",
            "C-1-C-1",
            "[C;H]",
            "[C;H0;H1]",
        ] {
            assert!(parse(pattern).is_err(), "unexpectedly accepted {pattern}");
        }
    }

    #[test]
    fn lays_out_imported_chains_rings_and_disconnected_components() {
        let chain = parse("C-C-C-C").unwrap();
        for bond in &chain.bonds {
            let length = chain.atoms[bond.a].pos.distance(chain.atoms[bond.b].pos);
            assert!(
                (length - crate::BOND).abs() < 2.0,
                "chain bond length {length}"
            );
        }

        let ring = parse("C-1-C-C-C-C-C-1").unwrap();
        assert_eq!(ring.bonds.len(), 6);
        for bond in &ring.bonds {
            let length = ring.atoms[bond.a].pos.distance(ring.atoms[bond.b].pos);
            assert!(
                (length - crate::BOND).abs() < 2.0,
                "ring bond length {length}"
            );
        }

        let fragments = parse("C.O").unwrap();
        assert!(fragments.atoms[0].pos.distance(fragments.atoms[1].pos) > 2.0 * crate::BOND);
    }
}
