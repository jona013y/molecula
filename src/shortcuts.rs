use crate::chemistry::valid_element;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QueryPredicate {
    AnyAtom,
    Heteroatom,
    ElementList(Vec<String>),
    ExcludedElementList(Vec<String>),
    Aromatic,
    Aliphatic,
    Unsaturated,
    HydrogenCount(u8),
    ImplicitHydrogenCount(u8),
    Connectivity(u8),
    Degree(u8),
    RingCount(Option<u8>),
    MinimumRingSize(u8),
    Substitutions(Option<u8>),
    Valence(u8),
    RingBonds(Option<u8>),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum QuickAction {
    Bond(u8),
    BondQuery(&'static str),
    Element(String),
    Charge(i8),
    RGroup(u16),
    AtomMap(u16),
    AssignUniqueAtomMaps,
    Isotope(Option<u16>),
    Alias(Option<String>),
    Query(String),
}
#[derive(Default)]
pub struct QuickInput {
    pub pending: String,
    last: f64,
}
impl QuickInput {
    pub fn push(&mut self, ch: char, time: f64) -> Vec<QuickAction> {
        let mut out = Vec::new();
        if !self.pending.is_empty() && time - self.last > 0.65 {
            out.extend(self.flush());
        }
        if ch.is_whitespace() {
            out.extend(self.flush());
            return out;
        }
        let candidate = format!("{}{}", self.pending, ch);
        if is_prefix(&candidate) {
            self.pending = candidate;
            self.last = time;
        } else {
            out.extend(self.flush());
            let single = ch.to_string();
            if is_prefix(&single) {
                self.pending = single;
                self.last = time;
            }
        }
        out
    }
    pub fn tick(&mut self, time: f64) -> Vec<QuickAction> {
        if !self.pending.is_empty() && time - self.last > 0.65 {
            self.flush()
        } else {
            vec![]
        }
    }
    pub fn flush(&mut self) -> Vec<QuickAction> {
        let text = std::mem::take(&mut self.pending);
        parse(&text).into_iter().collect()
    }
}
fn is_prefix(s: &str) -> bool {
    if s.is_empty() {
        return true;
    }
    if [
        "1", "2", "3", "4", "5", "6", "7", "0", "12", "14", "24", "+", "-", "*", "Q",
    ]
    .iter()
    .any(|x| x.starts_with(s))
    {
        return true;
    }
    if query_prefix(s) {
        return true;
    }
    if s.starts_with('R') && s[1..].chars().all(|c| c.is_ascii_digit()) {
        return true;
    }
    if s == "M+" || "M+".starts_with(s) || s == "M=" || "M=".starts_with(s) {
        return true;
    }
    if s.starts_with('M') && s[1..].chars().all(|c| c.is_ascii_digit()) {
        return true;
    }
    if let Some(mass) = s.strip_prefix('@') {
        return mass.chars().all(|c| c.is_ascii_digit());
    }
    if "label:".starts_with(s) {
        return true;
    }
    if let Some(label) = s.strip_prefix("label:") {
        return label.len() <= 24
            && label
                .chars()
                .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-');
    }
    if valid_element(s) {
        return true;
    }
    crate::chemistry::SYMBOLS.iter().any(|x| x.starts_with(s))
}
fn element_list_prefix(s: &str) -> bool {
    let s = s.strip_prefix('!').unwrap_or(s);
    if s.is_empty() {
        return true;
    }
    let values: Vec<_> = s.split(',').collect();
    values.iter().enumerate().all(|(index, value)| {
        if value.is_empty() {
            return index + 1 == values.len();
        }
        crate::chemistry::SYMBOLS
            .iter()
            .any(|symbol| symbol.starts_with(value))
    })
}
fn query_prefix(s: &str) -> bool {
    if s == "!" {
        return true;
    }
    if let Some(body) = s.strip_prefix('.') {
        return body.len() <= 7 && body.chars().all(|c| c.is_ascii_alphanumeric() || c == '*');
    }
    (s.contains(',') || s.starts_with('!')) && element_list_prefix(s)
}

pub fn parse_query_predicate(s: &str) -> Option<QueryPredicate> {
    Some(match s {
        "any atom" | "*" => QueryPredicate::AnyAtom,
        "heteroatom" | "Q" => QueryPredicate::Heteroatom,
        ".a" => QueryPredicate::Aromatic,
        ".A" => QueryPredicate::Aliphatic,
        ".u" => QueryPredicate::Unsaturated,
        _ if s.starts_with(".H") => QueryPredicate::HydrogenCount(s[2..].parse().ok()?),
        _ if s.starts_with(".h") => QueryPredicate::ImplicitHydrogenCount(s[2..].parse().ok()?),
        _ if s.starts_with(".X") => QueryPredicate::Connectivity(s[2..].parse().ok()?),
        _ if s.starts_with(".D") => QueryPredicate::Degree(s[2..].parse().ok()?),
        _ if s.starts_with(".R") => {
            let value = &s[2..];
            QueryPredicate::RingCount(if value == "*" {
                None
            } else {
                Some(value.parse().ok()?)
            })
        }
        _ if s.starts_with(".s") => {
            let value = &s[2..];
            QueryPredicate::Substitutions(if value == "*" {
                None
            } else {
                Some(value.parse().ok()?)
            })
        }
        _ if s.starts_with(".v") => QueryPredicate::Valence(s[2..].parse().ok()?),
        _ if s.starts_with(".rb") => {
            let value = &s[3..];
            QueryPredicate::RingBonds(if value == "*" {
                None
            } else {
                Some(value.parse().ok()?)
            })
        }
        _ if s.starts_with(".r") => QueryPredicate::MinimumRingSize(s[2..].parse().ok()?),
        _ if s.starts_with('!') || s.contains(',') => {
            let excluded = s.starts_with('!');
            let list = s.strip_prefix('!').unwrap_or(s);
            let mut elements = Vec::new();
            for element in list.split(',') {
                if !valid_element(element) || elements.iter().any(|seen| seen == element) {
                    return None;
                }
                elements.push(element.to_string());
            }
            if excluded {
                QueryPredicate::ExcludedElementList(elements)
            } else {
                QueryPredicate::ElementList(elements)
            }
        }
        _ => return None,
    })
}

fn query_family(predicate: &QueryPredicate) -> Option<&'static str> {
    Some(match predicate {
        QueryPredicate::Aromatic | QueryPredicate::Aliphatic => "aromaticity",
        QueryPredicate::HydrogenCount(_) => "hydrogen-count",
        QueryPredicate::ImplicitHydrogenCount(_) => "implicit-hydrogen-count",
        QueryPredicate::Connectivity(_) => "connectivity",
        QueryPredicate::Degree(_) => "degree",
        QueryPredicate::RingCount(_) => "ring-count",
        QueryPredicate::MinimumRingSize(_) => "minimum-ring-size",
        QueryPredicate::Substitutions(_) => "substitutions",
        QueryPredicate::Valence(_) => "valence",
        QueryPredicate::RingBonds(_) => "ring-bonds",
        _ => return None,
    })
}

pub fn merge_query_terms(current: Option<&str>, incoming: &str) -> Result<String, String> {
    let incoming_predicate = parse_query_predicate(incoming)
        .ok_or_else(|| format!("Unsupported atom query predicate: {incoming}"))?;
    let mut terms = current
        .filter(|value| !value.trim().is_empty())
        .map(|value| value.split(';').map(str::to_string).collect::<Vec<_>>())
        .unwrap_or_default();
    for term in &terms {
        if let Some(predicate) = parse_query_predicate(term) {
            if predicate == incoming_predicate {
                return Ok(terms.join(";"));
            }
            if query_family(&predicate).is_some()
                && query_family(&predicate) == query_family(&incoming_predicate)
            {
                return Err(format!(
                    "Conflicting query predicates: {term} and {incoming}"
                ));
            }
        }
    }
    let requested_elements = match &incoming_predicate {
        QueryPredicate::ElementList(elements) => Some(elements),
        _ => None,
    };
    let excluded_elements = match &incoming_predicate {
        QueryPredicate::ExcludedElementList(elements) => Some(elements),
        _ => None,
    };
    if let Some(elements) = requested_elements {
        if terms.iter().any(|term| {
            matches!(parse_query_predicate(term), Some(QueryPredicate::ExcludedElementList(ref excluded)) if elements.iter().all(|element| excluded.contains(element)))
        }) {
            return Err("Element query list is fully excluded by an existing predicate".into());
        }
    }
    if let Some(elements) = excluded_elements {
        if terms.iter().any(|term| {
            matches!(parse_query_predicate(term), Some(QueryPredicate::ElementList(ref included)) if included.iter().all(|element| elements.contains(element)))
        }) {
            return Err("Exclusion list removes every element in the existing query list".into());
        }
    }
    terms.push(incoming.to_string());
    Ok(terms.join(";"))
}
pub fn parse(s: &str) -> Option<QuickAction> {
    Some(match s {
        "1" => QuickAction::Bond(1),
        "2" => QuickAction::Bond(2),
        "3" => QuickAction::Bond(3),
        "4" => QuickAction::Bond(4),
        "5" => QuickAction::Bond(5),
        "6" => QuickAction::Bond(6),
        "7" => QuickAction::Bond(7),
        "0" => QuickAction::BondQuery("any"),
        "12" => QuickAction::BondQuery("single-or-double"),
        "14" => QuickAction::BondQuery("single-or-aromatic"),
        "24" => QuickAction::BondQuery("double-or-aromatic"),
        _ if s.starts_with('@') => {
            let mass = s[1..].parse::<u16>().ok()?;
            if mass == 0 {
                QuickAction::Isotope(None)
            } else if mass <= 300 {
                QuickAction::Isotope(Some(mass))
            } else {
                return None;
            }
        }
        "label:-" => QuickAction::Alias(None),
        _ if s.starts_with("label:") => {
            let label = &s[6..];
            if label.is_empty()
                || label.len() > 24
                || !label
                    .chars()
                    .all(|c| c.is_ascii_alphanumeric() || c == '_' || c == '-')
            {
                return None;
            }
            QuickAction::Alias(Some(label.to_string()))
        }
        "+" => QuickAction::Charge(1),
        "-" => QuickAction::Charge(-1),
        "*" => QuickAction::Query("any atom".into()),
        "Q" => QuickAction::Query("heteroatom".into()),
        _ if s.starts_with('R') && s.len() > 1 => {
            let n = s[1..].parse().ok()?;
            if n == 0 || n > 32767 {
                return None;
            }
            QuickAction::RGroup(n)
        }
        "M+" | "M=" => QuickAction::AssignUniqueAtomMaps,
        _ if s.starts_with('M') && s.len() > 1 => QuickAction::AtomMap(s[1..].parse().ok()?),
        _ if parse_query_predicate(s).is_some() => QuickAction::Query(s.into()),
        _ if valid_element(s) => QuickAction::Element(s.into()),
        _ => return None,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn multi_key_bond_does_not_commit_first_key() {
        let mut q = QuickInput::default();
        assert!(q.push('1', 0.).is_empty());
        assert!(q.push('2', 0.2).is_empty());
        assert_eq!(
            q.tick(0.9),
            vec![QuickAction::BondQuery("single-or-double")]
        );
    }
    #[test]
    fn composite_bond_sequences_wait_for_a_complete_token() {
        for (token, action) in [
            ("12", QuickAction::BondQuery("single-or-double")),
            ("14", QuickAction::BondQuery("single-or-aromatic")),
            ("24", QuickAction::BondQuery("double-or-aromatic")),
        ] {
            let mut input = QuickInput::default();
            for (index, ch) in token.chars().enumerate() {
                assert!(input.push(ch, index as f64 * 0.2).is_empty());
            }
            assert!(input.tick(0.7).is_empty());
            assert_eq!(input.tick(1.0), vec![action], "token {token}");
        }

        let mut input = QuickInput::default();
        assert!(input.push('1', 0.0).is_empty());
        assert!(input.tick(0.5).is_empty());
        assert_eq!(input.tick(0.7), vec![QuickAction::Bond(1)]);
    }
    #[test]
    fn atom_map_unique_shortcuts_wait_for_complete_tokens() {
        for token in ["M+", "M="] {
            let mut input = QuickInput::default();
            for (index, ch) in token.chars().enumerate() {
                assert!(input.push(ch, index as f64 * 0.1).is_empty());
            }
            assert_eq!(input.tick(0.9), vec![QuickAction::AssignUniqueAtomMaps]);
        }
        assert_eq!(parse("M23"), Some(QuickAction::AtomMap(23)));
        assert_eq!(parse("M0"), Some(QuickAction::AtomMap(0)));
    }
    #[test]
    fn isotope_mass_input_waits_for_complete_valid_token() {
        let mut input = QuickInput::default();
        for (index, ch) in "@13".chars().enumerate() {
            assert!(input.push(ch, index as f64 * 0.1).is_empty());
        }
        assert!(input.tick(0.7).is_empty());
        assert_eq!(input.tick(0.9), vec![QuickAction::Isotope(Some(13))]);
        let mut leading_zero = QuickInput::default();
        for (index, ch) in "@013".chars().enumerate() {
            assert!(leading_zero.push(ch, index as f64 * 0.1).is_empty());
        }
        assert_eq!(leading_zero.tick(1.0), vec![QuickAction::Isotope(Some(13))]);
        assert_eq!(parse("@0"), Some(QuickAction::Isotope(None)));
        assert_eq!(parse("@301"), None);
        assert_eq!(parse("@"), None);
    }
    #[test]
    fn alias_input_waits_for_a_complete_label_and_supports_clear() {
        let mut input = QuickInput::default();
        for (index, ch) in "label:Ph".chars().enumerate() {
            assert!(input.push(ch, index as f64 * 0.05).is_empty());
        }
        assert_eq!(input.tick(1.1), vec![QuickAction::Alias(Some("Ph".into()))]);
        assert_eq!(parse("label:-"), Some(QuickAction::Alias(None)));
        assert_eq!(parse("label:"), None);
        assert_eq!(parse("label:very_long_label_name_overflow"), None);
    }
    #[test]
    fn cl_is_not_two_elements() {
        let mut q = QuickInput::default();
        assert!(q.push('C', 0.).is_empty());
        assert!(q.push('l', 0.1).is_empty());
        assert_eq!(q.tick(0.8), vec![QuickAction::Element("Cl".into())]);
    }
    #[test]
    fn incomplete_group_is_not_action() {
        assert_eq!(parse("R"), None);
        assert_eq!(parse("R32768"), None);
    }

    #[test]
    fn parses_search_predicates_and_element_lists_without_guessing() {
        for token in [
            ".a",
            ".A",
            ".u",
            ".H0",
            ".H4",
            ".h1",
            ".X3",
            ".D2",
            ".R0",
            ".R1",
            ".r5",
            ".s*",
            ".s2",
            ".v3",
            ".rb*",
            ".rb2",
            "Au,Ag,Pt",
            "!Au,Ag,Pt",
            "Q",
            "*",
        ] {
            assert!(parse_query_predicate(token).is_some(), "{token}");
            assert!(
                matches!(parse(token), Some(QuickAction::Query(_))),
                "{token}"
            );
        }
        assert_eq!(parse(".H"), None);
        assert_eq!(parse(".X-1"), None);
        assert_eq!(parse(".Rfoo"), None);
        assert_eq!(parse("Au,,Pt"), None);
        assert_eq!(parse("Fe"), Some(QuickAction::Element("Fe".into())));
    }

    #[test]
    fn quick_input_buffers_element_lists_as_a_single_query() {
        let mut q = QuickInput::default();
        for (index, ch) in "Au,Ag,Pt".chars().enumerate() {
            assert!(q.push(ch, index as f64 * 0.05).is_empty());
        }
        assert_eq!(q.tick(1.1), vec![QuickAction::Query("Au,Ag,Pt".into())]);
    }

    #[test]
    fn query_predicates_compose_idempotently_and_reject_conflicts() {
        let query = merge_query_terms(None, ".a").unwrap();
        let query = merge_query_terms(Some(&query), ".H0").unwrap();
        assert_eq!(query, ".a;.H0");
        assert_eq!(merge_query_terms(Some(&query), ".H0").unwrap(), query);
        assert!(merge_query_terms(Some(&query), ".A").is_err());
        assert!(merge_query_terms(Some(&query), ".H2").is_err());
        let elements = merge_query_terms(None, "Au,Ag").unwrap();
        assert!(merge_query_terms(Some(&elements), "!Au,Ag").is_err());
    }
}
