//! Small, explicit substructure matcher for Molecula's supported query subset.
//!
//! This is deliberately not advertised as full SMARTS: recursive queries and
//! stereochemical matching beyond explicit stored CIP/E/Z labels are rejected. Ring counts use a deterministic local
//! minimum cycle basis, whose tie-breaking can differ from other toolkits.
use crate::{chemistry, shortcuts, Bond, BondStyle, Molecule};
use std::collections::{BTreeSet, HashSet};

const MAX_MATCHES: usize = 1000;
const MAX_RING_SEARCH_WORK: usize = 20_000_000;
const MAX_RING_CANDIDATES: usize = 100_000;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MatchReport {
    /// Unique target atom sets. Symmetry-equivalent permutations count once.
    pub matches: Vec<Vec<usize>>,
    pub query_atoms: usize,
    pub target_atoms: usize,
}

fn bond_order(bond: &Bond) -> u8 {
    if bond.style == BondStyle::Aromatic {
        4
    } else {
        bond.order
    }
}

fn bond_matches(query: &Bond, target: &Bond, target_is_ring: bool) -> bool {
    let Some(query_text) = query.query.as_deref() else {
        return if query.style == BondStyle::Any {
            true
        } else if query.style == BondStyle::Aromatic || query.order == 4 {
            bond_order(target) == 4
        } else {
            query.order == target.order && target.style != BondStyle::Any
        };
    };
    let mut has_order_predicate = false;
    for predicate in query_text.split(';') {
        let matches = match predicate {
            "any" => {
                has_order_predicate = true;
                true
            }
            "single-or-double" => {
                has_order_predicate = true;
                matches!(bond_order(target), 1 | 2)
            }
            "single-or-aromatic" => {
                has_order_predicate = true;
                matches!(bond_order(target), 1 | 4)
            }
            "double-or-aromatic" => {
                has_order_predicate = true;
                matches!(bond_order(target), 2 | 4)
            }
            "topology-ring" => target_is_ring,
            "topology-chain" => !target_is_ring,
            _ => return false,
        };
        if !matches {
            return false;
        }
    }
    has_order_predicate
        || query_text
            .split(';')
            .any(|term| matches!(term, "topology-ring" | "topology-chain"))
}

fn bond_stereo_matches(query: &Bond, target: &Bond) -> Result<bool, String> {
    let Some(stereo) = query.stereo.as_deref() else {
        return Ok(true);
    };
    if query.order != 2 || !matches!(stereo, "E" | "Z" | "cis" | "trans") {
        return Err(
            "Only explicit E/Z or cis/trans labels on double-bond queries are supported".into(),
        );
    }
    Ok(target.order == 2 && target.stereo.as_deref() == Some(stereo))
}

fn bond_between(molecule: &Molecule, a: usize, b: usize) -> Option<(usize, &Bond)> {
    molecule
        .bonds
        .iter()
        .enumerate()
        .find(|(_, bond)| (bond.a == a && bond.b == b) || (bond.a == b && bond.b == a))
}

/// Count incident bonds which belong to at least one cycle, using an iterative
/// Tarjan bridge pass. Iteration avoids stack overflow on long drawn chains.
fn ring_bond_counts(
    adjacency: &[Vec<(usize, usize)>],
    target: &HashSet<usize>,
) -> (Vec<u8>, HashSet<usize>) {
    let n = adjacency.len();
    let unseen = usize::MAX;
    let mut discovered = vec![unseen; n];
    let mut low = vec![0usize; n];
    let mut parent = vec![unseen; n];
    let mut parent_edge = vec![unseen; n];
    let mut next_neighbor = vec![0usize; n];
    let mut bridges = HashSet::new();
    let mut clock = 0usize;
    for &root in target {
        if discovered[root] != unseen {
            continue;
        }
        discovered[root] = clock;
        low[root] = clock;
        clock += 1;
        let mut stack = vec![root];
        while let Some(&v) = stack.last() {
            if next_neighbor[v] < adjacency[v].len() {
                let (w, edge) = adjacency[v][next_neighbor[v]];
                next_neighbor[v] += 1;
                if !target.contains(&w) || edge == parent_edge[v] {
                    continue;
                }
                if discovered[w] == unseen {
                    parent[w] = v;
                    parent_edge[w] = edge;
                    discovered[w] = clock;
                    low[w] = clock;
                    clock += 1;
                    stack.push(w);
                } else {
                    low[v] = low[v].min(discovered[w]);
                }
            } else {
                stack.pop();
                if parent[v] != unseen {
                    let p = parent[v];
                    low[p] = low[p].min(low[v]);
                    if low[v] > discovered[p] {
                        bridges.insert(parent_edge[v]);
                    }
                }
            }
        }
    }
    let mut counts = vec![0u8; n];
    for &v in target {
        counts[v] = adjacency[v]
            .iter()
            .filter(|(w, edge)| target.contains(w) && !bridges.contains(edge))
            .count()
            .min(u8::MAX as usize) as u8;
    }
    let ring_edges = adjacency
        .iter()
        .enumerate()
        .filter(|(atom, _)| target.contains(atom))
        .flat_map(|(_, edges)| edges.iter())
        .filter_map(|(neighbor, edge)| {
            (target.contains(neighbor) && !bridges.contains(edge)).then_some(*edge)
        })
        .collect();
    (counts, ring_edges)
}

fn xor_sorted(left: &[usize], right: &[usize]) -> Vec<usize> {
    let mut out = Vec::with_capacity(left.len() + right.len());
    let (mut i, mut j) = (0, 0);
    while i < left.len() && j < right.len() {
        match left[i].cmp(&right[j]) {
            std::cmp::Ordering::Less => {
                out.push(left[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                out.push(right[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                i += 1;
                j += 1;
            }
        }
    }
    out.extend_from_slice(&left[i..]);
    out.extend_from_slice(&right[j..]);
    out
}

/// Build a deterministic minimum cycle basis from Horton candidates and count
/// how many selected basis cycles contain each atom. This gives `.R0`, `.R1`,
/// etc. a concrete local meaning without claiming another toolkit's SSSR.
fn ring_basis_membership(
    adjacency: &[Vec<(usize, usize)>],
    bonds: &[Bond],
    target: &HashSet<usize>,
) -> Result<Vec<u8>, String> {
    let mut edge_count = 0usize;
    for &v in target {
        edge_count += adjacency[v]
            .iter()
            .filter(|(w, _)| target.contains(w))
            .count();
    }
    edge_count /= 2;
    if edge_count == 0 {
        return Ok(vec![0; adjacency.len()]);
    }

    let mut components = 0usize;
    let mut seen = HashSet::new();
    for &root in target {
        if !seen.insert(root) {
            continue;
        }
        components += 1;
        let mut todo = vec![root];
        while let Some(v) = todo.pop() {
            for &(w, _) in &adjacency[v] {
                if target.contains(&w) && seen.insert(w) {
                    todo.push(w);
                }
            }
        }
    }
    let rank = edge_count + components - target.len();
    if rank == 0 {
        return Ok(vec![0; adjacency.len()]);
    }

    let mut candidates = BTreeSet::<Vec<usize>>::new();
    let mut work = 0usize;
    for &root in target {
        work = work.saturating_add(target.len()); // BFS array initialization cost
        if work > MAX_RING_SEARCH_WORK {
            return Err("Ring-count query exceeds the local cycle-basis work budget for this large drawing; reduce the target size".into());
        }
        let mut parent = vec![usize::MAX; adjacency.len()];
        let mut parent_edge = vec![usize::MAX; adjacency.len()];
        let mut distance = vec![usize::MAX; adjacency.len()];
        let mut queue = std::collections::VecDeque::new();
        distance[root] = 0;
        queue.push_back(root);
        while let Some(v) = queue.pop_front() {
            for &(w, edge) in &adjacency[v] {
                if !target.contains(&w) {
                    continue;
                }
                work = work.saturating_add(1);
                if distance[w] == usize::MAX {
                    distance[w] = distance[v] + 1;
                    parent[w] = v;
                    parent_edge[w] = edge;
                    queue.push_back(w);
                }
            }
        }
        for &v in target {
            if distance[v] == usize::MAX {
                continue;
            }
            for &(w, edge) in &adjacency[v] {
                if !target.contains(&w)
                    || v >= w
                    || edge == parent_edge[v]
                    || edge == parent_edge[w]
                {
                    continue;
                }
                let mut cycle = HashSet::new();
                let mut node = v;
                while parent_edge[node] != usize::MAX {
                    let path_edge = parent_edge[node];
                    if !cycle.insert(path_edge) {
                        cycle.remove(&path_edge);
                    }
                    node = parent[node];
                }
                let mut node = w;
                while parent_edge[node] != usize::MAX {
                    let path_edge = parent_edge[node];
                    if !cycle.insert(path_edge) {
                        cycle.remove(&path_edge);
                    }
                    node = parent[node];
                }
                if !cycle.insert(edge) {
                    cycle.remove(&edge);
                }
                if cycle.len() < 3 {
                    continue;
                }
                let mut cycle: Vec<usize> = cycle.into_iter().collect();
                cycle.sort_unstable();
                work = work.saturating_add(cycle.len());
                if work > MAX_RING_SEARCH_WORK {
                    return Err("Ring-count query exceeds the local cycle-basis work budget for this large drawing; reduce the target size".into());
                }
                candidates.insert(cycle);
                if candidates.len() > MAX_RING_CANDIDATES {
                    return Err(
                        "Ring-count query generated too many cycle candidates for this drawing"
                            .into(),
                    );
                }
            }
        }
        if work > MAX_RING_SEARCH_WORK {
            return Err("Ring-count query exceeds the local cycle-basis work budget for this large drawing; reduce the target size".into());
        }
    }

    let mut ordered: Vec<Vec<usize>> = candidates.into_iter().collect();
    ordered.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));
    let mut basis = std::collections::HashMap::<usize, Vec<usize>>::new();
    let mut memberships = vec![0u8; adjacency.len()];
    let mut selected_cycles = 0usize;
    for cycle in ordered {
        let original = cycle.clone();
        let mut reduced = cycle;
        while let Some(&pivot) = reduced.last() {
            if let Some(row) = basis.get(&pivot) {
                reduced = xor_sorted(&reduced, row);
            } else {
                basis.insert(pivot, reduced);
                let mut atoms = BTreeSet::new();
                for edge in original {
                    if let Some(bond) = bonds.get(edge) {
                        atoms.insert(bond.a);
                        atoms.insert(bond.b);
                    }
                }
                for atom in atoms {
                    memberships[atom] = memberships[atom].saturating_add(1);
                }
                selected_cycles += 1;
                break;
            }
        }
        if selected_cycles == rank {
            return Ok(memberships);
        }
    }
    Err("Could not construct a complete local minimum cycle basis; ring-count matching was not performed".into())
}

/// Exact smallest simple cycle through each atom. A cycle through `v` is the
/// shortest path between a pair of distinct neighbors after removing `v`, +2.
fn minimum_ring_sizes(
    adjacency: &[Vec<(usize, usize)>],
    target: &HashSet<usize>,
) -> Result<Vec<Option<u8>>, String> {
    let mut work = 0usize;
    let mut results = vec![None; adjacency.len()];
    for &root in target {
        let neighbors: Vec<usize> = adjacency[root]
            .iter()
            .map(|(neighbor, _)| *neighbor)
            .filter(|neighbor| target.contains(neighbor))
            .collect();
        let pairs = neighbors
            .len()
            .saturating_mul(neighbors.len().saturating_sub(1))
            / 2;
        work = work.saturating_add(
            pairs.saturating_mul(target.len() + adjacency.iter().map(Vec::len).sum::<usize>()),
        );
        if work > MAX_RING_SEARCH_WORK {
            return Err("Minimum-ring query exceeds the local search budget for this large drawing; reduce the target size".into());
        }
        for a in 0..neighbors.len() {
            for b in a + 1..neighbors.len() {
                let start = neighbors[a];
                let end = neighbors[b];
                let mut distance = vec![usize::MAX; adjacency.len()];
                let mut queue = std::collections::VecDeque::new();
                distance[start] = 0;
                queue.push_back(start);
                while let Some(v) = queue.pop_front() {
                    if v == end {
                        break;
                    }
                    for &(w, _) in &adjacency[v] {
                        if w == root || !target.contains(&w) || distance[w] != usize::MAX {
                            continue;
                        }
                        distance[w] = distance[v] + 1;
                        queue.push_back(w);
                    }
                }
                if distance[end] != usize::MAX {
                    let size = distance[end] + 2;
                    let size = size.min(u8::MAX as usize) as u8;
                    results[root] = Some(results[root].map_or(size, |old: u8| old.min(size)));
                }
            }
        }
    }
    Ok(results)
}

struct TargetAtomFeatures {
    aromatic: bool,
    aliphatic: bool,
    degree: usize,
    valence: u16,
    implicit_h: u8,
    total_h: u16,
    unsaturated: bool,
    minimum_ring_size: Option<u8>,
    ring_count: u8,
    ring_bonds: u8,
    substitutions: u8,
    query_substitutions: u8,
}

fn atom_matches(
    query: &crate::Atom,
    target: &crate::Atom,
    features: &TargetAtomFeatures,
) -> Result<bool, String> {
    if let Some(stereo) = query.stereo.as_deref() {
        if !matches!(stereo, "R" | "S") {
            return Err(
                "Only explicit R/S CIP labels are supported for stereochemical atom matching"
                    .into(),
            );
        }
        if target.stereo.as_deref() != Some(stereo) {
            return Ok(false);
        }
    }
    if query.r_group.is_some() || query.alias.is_some() {
        return Err("R-group and alias query atoms are not supported by the local matcher".into());
    }
    let terms: Vec<&str> = query
        .query
        .as_deref()
        .unwrap_or("")
        .split(';')
        .filter(|s| !s.is_empty())
        .collect();
    let has_wildcard_element_test = terms.iter().any(|term| {
        matches!(
            shortcuts::parse_query_predicate(term),
            Some(
                shortcuts::QueryPredicate::ElementList(_)
                    | shortcuts::QueryPredicate::ExcludedElementList(_)
                    | shortcuts::QueryPredicate::AnyAtom
                    | shortcuts::QueryPredicate::Heteroatom
            )
        )
    });
    if !has_wildcard_element_test {
        if !chemistry::valid_element(&query.element) {
            return Err(format!("Unsupported query atom label: {}", query.element));
        }
        if query.element != target.element {
            return Ok(false);
        }
        if query.charge != target.charge
            || query
                .isotope
                .is_some_and(|isotope| target.isotope != Some(isotope))
        {
            return Ok(false);
        }
    }
    for term in terms {
        let predicate = shortcuts::parse_query_predicate(term)
            .ok_or_else(|| format!("Unsupported query predicate: {term}"))?;
        let matches = match predicate {
            shortcuts::QueryPredicate::AnyAtom => true,
            shortcuts::QueryPredicate::Heteroatom => target.element != "C" && target.element != "H",
            shortcuts::QueryPredicate::ElementList(elements) => elements.contains(&target.element),
            shortcuts::QueryPredicate::ExcludedElementList(elements) => {
                !elements.contains(&target.element)
            }
            shortcuts::QueryPredicate::Aromatic => features.aromatic,
            shortcuts::QueryPredicate::Aliphatic => features.aliphatic,
            shortcuts::QueryPredicate::Unsaturated => features.unsaturated,
            shortcuts::QueryPredicate::HydrogenCount(n) => features.total_h == n as u16,
            shortcuts::QueryPredicate::ImplicitHydrogenCount(n) => features.implicit_h == n,
            shortcuts::QueryPredicate::Connectivity(n) => {
                features.degree + features.total_h as usize == n as usize
            }
            shortcuts::QueryPredicate::Degree(n) => features.degree == n as usize,
            shortcuts::QueryPredicate::Valence(n) => features.valence == n as u16,
            shortcuts::QueryPredicate::RingCount(None) => true,
            shortcuts::QueryPredicate::RingCount(Some(0)) => features.ring_bonds == 0,
            shortcuts::QueryPredicate::RingCount(Some(count)) => features.ring_count == count,
            shortcuts::QueryPredicate::MinimumRingSize(size) => {
                features.minimum_ring_size == Some(size)
            }
            shortcuts::QueryPredicate::RingBonds(None) => features.ring_bonds > 0,
            shortcuts::QueryPredicate::RingBonds(Some(size)) => features.ring_bonds == size,
            shortcuts::QueryPredicate::Substitutions(None) => {
                features.substitutions == features.query_substitutions
            }
            shortcuts::QueryPredicate::Substitutions(Some(count)) => {
                features.substitutions == count
            }
        };
        if !matches {
            return Ok(false);
        }
    }
    Ok(true)
}

/// Match the selected connected graph against all other connected components
/// in the same document. The selected query fragment is excluded from targets.
pub fn find_selected_matches(
    molecule: &Molecule,
    selected: &[usize],
) -> Result<MatchReport, String> {
    let selected: BTreeSet<usize> = selected.iter().copied().collect();
    if selected.is_empty() {
        return Err("Select a connected query fragment first".into());
    }
    if selected.iter().any(|&i| i >= molecule.atoms.len()) {
        return Err("Selection contains an invalid atom index".into());
    }
    let mut edge_keys = HashSet::new();
    for bond in &molecule.bonds {
        if bond.a >= molecule.atoms.len() || bond.b >= molecule.atoms.len() {
            return Err("Molecule contains an invalid bond endpoint".into());
        }
        if bond.a == bond.b || !edge_keys.insert((bond.a.min(bond.b), bond.a.max(bond.b))) {
            return Err("Molecule contains a self-loop or duplicate bond".into());
        }
        if !(1..=4).contains(&bond.order) {
            return Err("Molecule contains an unsupported bond order".into());
        }
    }
    let mut query_bonds = Vec::new();
    let mut q_adj = vec![Vec::<usize>::new(); molecule.atoms.len()];
    for (index, bond) in molecule.bonds.iter().enumerate() {
        if selected.contains(&bond.a) && selected.contains(&bond.b) {
            if !matches!(
                bond.style,
                BondStyle::Normal | BondStyle::Aromatic | BondStyle::Any
            ) {
                return Err("Stereochemical bond matching is not implemented".into());
            }
            q_adj[bond.a].push(bond.b);
            q_adj[bond.b].push(bond.a);
            query_bonds.push((index, bond.a, bond.b));
        }
    }
    let start = *selected.first().unwrap();
    let mut visited = HashSet::new();
    let mut stack = vec![start];
    visited.insert(start);
    while let Some(atom) = stack.pop() {
        for &next in &q_adj[atom] {
            if visited.insert(next) {
                stack.push(next);
            }
        }
    }
    if visited.len() != selected.len() {
        return Err("Selected query atoms must form one connected fragment".into());
    }
    if molecule
        .bonds
        .iter()
        .any(|bond| selected.contains(&bond.a) != selected.contains(&bond.b))
    {
        return Err("Select the complete query molecule; draw target molecules as separate disconnected components".into());
    }

    let mut adjacency = vec![Vec::<(usize, usize)>::new(); molecule.atoms.len()];
    let mut bond_valence = vec![0u16; molecule.atoms.len()];
    let mut unsaturated = vec![false; molecule.atoms.len()];
    for (bi, bond) in molecule.bonds.iter().enumerate() {
        if bond.a >= molecule.atoms.len() || bond.b >= molecule.atoms.len() {
            return Err("Molecule contains an invalid bond endpoint".into());
        }
        adjacency[bond.a].push((bond.b, bi));
        adjacency[bond.b].push((bond.a, bi));
        let order = match bond.style {
            BondStyle::Coordinate => 0,
            BondStyle::Aromatic => 1,
            _ => bond_order(bond) as u16,
        };
        bond_valence[bond.a] += order.min(3);
        bond_valence[bond.b] += order.min(3);
        if order > 1 {
            unsaturated[bond.a] = true;
            unsaturated[bond.b] = true;
        }
    }
    let analysis = chemistry::analyze(molecule);
    let mut target_components = Vec::<Vec<usize>>::new();
    let mut seen = HashSet::new();
    for root in 0..molecule.atoms.len() {
        if selected.contains(&root) || !seen.insert(root) {
            continue;
        }
        let mut component = vec![root];
        let mut todo = vec![root];
        while let Some(i) = todo.pop() {
            for &(j, _) in &adjacency[i] {
                if !selected.contains(&j) && seen.insert(j) {
                    component.push(j);
                    todo.push(j);
                }
            }
        }
        target_components.push(component);
    }
    let target_count: usize = target_components.iter().map(Vec::len).sum();
    if target_count == 0 {
        return Err("Draw at least one separate target molecule beside the selected query".into());
    }

    let target_atom_set: HashSet<usize> = target_components.iter().flatten().copied().collect();
    let (ring_bond_count, ring_edges) = ring_bond_counts(&adjacency, &target_atom_set);
    let needs_ring_basis = selected.iter().any(|&index| {
        molecule.atoms[index]
            .query
            .as_deref()
            .unwrap_or("")
            .split(';')
            .any(|term| {
                matches!(
                    shortcuts::parse_query_predicate(term),
                    Some(shortcuts::QueryPredicate::RingCount(Some(count))) if count > 0
                )
            })
    });
    let ring_count = if needs_ring_basis {
        ring_basis_membership(&adjacency, &molecule.bonds, &target_atom_set)?
    } else {
        vec![0; molecule.atoms.len()]
    };
    let needs_minimum_ring = selected.iter().any(|&index| {
        molecule.atoms[index]
            .query
            .as_deref()
            .unwrap_or("")
            .split(';')
            .any(|term| {
                matches!(
                    shortcuts::parse_query_predicate(term),
                    Some(shortcuts::QueryPredicate::MinimumRingSize(_))
                )
            })
    });
    let minimum_ring_size = if needs_minimum_ring {
        minimum_ring_sizes(&adjacency, &target_atom_set)?
    } else {
        vec![None; molecule.atoms.len()]
    };

    let q_atoms: Vec<usize> = selected.iter().copied().collect();
    let mut candidates = std::collections::HashMap::<usize, Vec<usize>>::new();
    for &qi in &q_atoms {
        let mut list = Vec::new();
        for component in &target_components {
            for &ti in component {
                let implicit_h_value = analysis.implicit_h.get(ti).copied().flatten();
                let needs_hydrogen_count = molecule.atoms[qi]
                    .query
                    .as_deref()
                    .unwrap_or("")
                    .split(';')
                    .any(|term| {
                        matches!(term, ".H" | ".h")
                            || term.starts_with(".H")
                            || term.starts_with(".h")
                            || term.starts_with(".X")
                            || term.starts_with(".v")
                    });
                if needs_hydrogen_count && implicit_h_value.is_none() {
                    return Err(format!("Cannot evaluate hydrogen/valence predicate for target atom {} because its implicit hydrogen count is unresolved", ti + 1));
                }
                let implicit_h = implicit_h_value.unwrap_or(0);
                let bonded_h = adjacency[ti]
                    .iter()
                    .filter(|(neighbor, _)| molecule.atoms[*neighbor].element == "H")
                    .count() as u16;
                let total_h = molecule.atoms[ti].explicit_h as u16 + implicit_h as u16 + bonded_h;
                let heavy_degree = adjacency[ti]
                    .iter()
                    .filter(|(neighbor, _)| molecule.atoms[*neighbor].element != "H")
                    .count();
                let isotope_h_neighbors = adjacency[ti]
                    .iter()
                    .filter(|(neighbor, _)| {
                        molecule.atoms[*neighbor].element == "H"
                            && molecule.atoms[*neighbor].isotope.is_some()
                    })
                    .count();
                let query_substitutions = adjacency[qi]
                    .iter()
                    .filter(|(neighbor, _)| {
                        selected.contains(neighbor)
                            && (molecule.atoms[*neighbor].element != "H"
                                || molecule.atoms[*neighbor].isotope.is_some())
                    })
                    .count();
                let has_aromatic_bond = adjacency[ti].iter().any(|(_, bond_index)| {
                    let bond = &molecule.bonds[*bond_index];
                    bond.style == BondStyle::Aromatic || bond.order == 4
                });
                let features = TargetAtomFeatures {
                    aromatic: molecule.atoms[ti].aromatic || has_aromatic_bond,
                    aliphatic: !molecule.atoms[ti].aromatic && !has_aromatic_bond,
                    degree: heavy_degree,
                    valence: bond_valence[ti] + implicit_h as u16,
                    implicit_h,
                    total_h,
                    unsaturated: unsaturated[ti],
                    minimum_ring_size: minimum_ring_size[ti],
                    ring_count: ring_count[ti],
                    ring_bonds: ring_bond_count[ti],
                    substitutions: (heavy_degree + isotope_h_neighbors).min(u8::MAX as usize) as u8,
                    query_substitutions: query_substitutions.min(u8::MAX as usize) as u8,
                };
                if atom_matches(&molecule.atoms[qi], &molecule.atoms[ti], &features)? {
                    list.push(ti);
                }
            }
        }
        if list.is_empty() {
            return Ok(MatchReport {
                matches: vec![],
                query_atoms: selected.len(),
                target_atoms: target_count,
            });
        }
        candidates.insert(qi, list);
    }
    let mut order = q_atoms.clone();
    order.sort_by_key(|qi| (candidates[qi].len(), std::cmp::Reverse(q_adj[*qi].len())));
    let mut mapping = vec![usize::MAX; molecule.atoms.len()];
    let mut used = HashSet::new();
    let mut unique = BTreeSet::<Vec<usize>>::new();

    // Recursion state is explicit here so each backtracking call has no hidden
    // global mutation; the argument list reflects the search invariants.
    #[allow(clippy::too_many_arguments)]
    fn visit(
        depth: usize,
        order: &[usize],
        candidates: &std::collections::HashMap<usize, Vec<usize>>,
        mapping: &mut [usize],
        used: &mut HashSet<usize>,
        selected: &BTreeSet<usize>,
        target_atoms: &HashSet<usize>,
        ring_edges: &HashSet<usize>,
        query_bonds: &[(usize, usize, usize)],
        molecule: &Molecule,
        unique: &mut BTreeSet<Vec<usize>>,
    ) -> Result<(), String> {
        if unique.len() >= MAX_MATCHES {
            return Err(format!(
                "Match limit ({MAX_MATCHES}) reached; narrow the query"
            ));
        }
        if depth == order.len() {
            let mut atom_set: Vec<usize> = selected.iter().map(|&q| mapping[q]).collect();
            atom_set.sort_unstable();
            unique.insert(atom_set);
            return Ok(());
        }
        let qi = order[depth];
        for &ti in &candidates[&qi] {
            if !target_atoms.contains(&ti) || !used.insert(ti) {
                continue;
            }
            mapping[qi] = ti;
            let mut compatible = true;
            for &(bi, qa, qb) in query_bonds {
                let (other_q, this_q) = if qa == qi {
                    (qb, qi)
                } else if qb == qi {
                    (qa, qi)
                } else {
                    continue;
                };
                if mapping[other_q] != usize::MAX {
                    let target_bond = bond_between(molecule, mapping[this_q], mapping[other_q]);
                    if target_bond.is_none_or(|(bond_index, bond)| {
                        !bond_matches(&molecule.bonds[bi], bond, ring_edges.contains(&bond_index))
                    }) {
                        compatible = false;
                        break;
                    }
                    if let Some((_, target_bond)) = target_bond {
                        if !bond_stereo_matches(&molecule.bonds[bi], target_bond)? {
                            compatible = false;
                            break;
                        }
                    }
                }
            }
            if compatible {
                visit(
                    depth + 1,
                    order,
                    candidates,
                    mapping,
                    used,
                    selected,
                    target_atoms,
                    ring_edges,
                    query_bonds,
                    molecule,
                    unique,
                )?;
            }
            mapping[qi] = usize::MAX;
            used.remove(&ti);
        }
        Ok(())
    }
    visit(
        0,
        &order,
        &candidates,
        &mut mapping,
        &mut used,
        &selected,
        &target_atom_set,
        &ring_edges,
        &query_bonds,
        molecule,
        &mut unique,
    )?;
    Ok(MatchReport {
        matches: unique.into_iter().collect(),
        query_atoms: selected.len(),
        target_atoms: target_count,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Atom;
    use macroquad::prelude::vec2;

    fn atom(element: &str, query: Option<&str>) -> Atom {
        let mut a = Atom::new(vec2(0., 0.), element);
        a.query = query.map(str::to_owned);
        a
    }
    fn bond(a: usize, b: usize, order: u8, query: Option<&str>) -> Bond {
        Bond {
            a,
            b,
            order,
            style: BondStyle::Normal,
            stereo: None,
            query: query.map(str::to_owned),
        }
    }

    #[test]
    fn finds_query_fragment_in_disconnected_target_and_deduplicates_symmetry() {
        let mut m = Molecule::default();
        // Query C-O fragment; target is ethanol, so the terminal C-O occurrence matches.
        m.atoms.extend([
            atom("C", None),
            atom("O", None),
            atom("C", None),
            atom("C", None),
            atom("O", None),
        ]);
        m.bonds.extend([
            bond(0, 1, 1, None),
            bond(2, 3, 1, None),
            bond(3, 4, 1, None),
        ]);
        let report = find_selected_matches(&m, &[0, 1]).unwrap();
        assert_eq!(report.matches, vec![vec![3, 4]]);
    }

    #[test]
    fn evaluates_supported_element_and_hydrogen_query_predicates() {
        let mut m = Molecule::default();
        let query = atom("C", Some(".H3"));
        m.atoms.extend([query, atom("C", None), atom("C", None)]);
        m.bonds.push(bond(1, 2, 1, None));
        assert_eq!(
            find_selected_matches(&m, &[0]).unwrap().matches,
            vec![vec![1], vec![2]]
        );
    }

    #[test]
    fn atom_lists_wildcards_and_heteroatoms_match_their_defined_sets() {
        let matches_for = |predicate: &str| {
            let mut molecule = Molecule::default();
            molecule.atoms.push(atom("C", Some(predicate)));
            molecule.atoms.extend(
                ["C", "N", "O", "H", "Cl"]
                    .into_iter()
                    .map(|e| atom(e, None)),
            );
            find_selected_matches(&molecule, &[0]).unwrap().matches
        };
        assert_eq!(matches_for("N,O"), vec![vec![2], vec![3]]);
        assert_eq!(matches_for("!O"), vec![vec![1], vec![2], vec![4], vec![5]]);
        assert_eq!(matches_for("heteroatom"), vec![vec![2], vec![3], vec![5]]);
        assert_eq!(
            matches_for("any atom"),
            vec![vec![1], vec![2], vec![3], vec![4], vec![5]]
        );
    }

    #[test]
    fn explicit_cip_query_labels_match_only_identically_labelled_atoms() {
        let mut molecule = Molecule::default();
        let mut query = atom("C", None);
        query.stereo = Some("R".into());
        molecule.atoms.push(query);
        let mut r_target = atom("C", None);
        r_target.stereo = Some("R".into());
        let mut s_target = atom("C", None);
        s_target.stereo = Some("S".into());
        molecule.atoms.extend([r_target, s_target, atom("C", None)]);
        assert_eq!(
            find_selected_matches(&molecule, &[0]).unwrap().matches,
            vec![vec![1]]
        );

        molecule.atoms[0].stereo = Some("ABS".into());
        assert!(find_selected_matches(&molecule, &[0])
            .unwrap_err()
            .contains("explicit R/S CIP labels"));
    }

    #[test]
    fn explicit_double_bond_stereo_query_matches_only_the_same_label() {
        let mut molecule = Molecule::default();
        molecule.atoms.extend((0..8).map(|_| atom("C", None)));
        molecule.bonds.extend([
            bond(0, 1, 2, None),
            bond(2, 3, 2, None),
            bond(4, 5, 2, None),
            bond(6, 7, 2, None),
        ]);
        molecule.bonds[0].stereo = Some("E".into());
        molecule.bonds[1].stereo = Some("E".into());
        molecule.bonds[2].stereo = Some("Z".into());
        assert_eq!(
            find_selected_matches(&molecule, &[0, 1]).unwrap().matches,
            vec![vec![2, 3]]
        );

        molecule.bonds[0].stereo = Some("cis".into());
        assert!(find_selected_matches(&molecule, &[0, 1])
            .unwrap()
            .matches
            .is_empty());
    }

    #[test]
    fn ring_membership_uses_deterministic_minimum_cycle_basis_for_fused_rings() {
        let mut m = Molecule::default();
        m.atoms.push(atom("C", Some(".R2")));
        m.atoms.extend((0..10).map(|_| atom("C", None)));
        for vertex in 0..10 {
            m.bonds
                .push(bond(vertex + 1, (vertex + 1) % 10 + 1, 1, None));
        }
        m.bonds.push(bond(1, 6, 1, None)); // shared edge of two six-membered rings
        assert_eq!(
            find_selected_matches(&m, &[0]).unwrap().matches,
            vec![vec![1], vec![6]]
        );
        m.atoms[0].query = Some(".R1".into());
        assert_eq!(
            find_selected_matches(&m, &[0]).unwrap().matches,
            vec![
                vec![2],
                vec![3],
                vec![4],
                vec![5],
                vec![7],
                vec![8],
                vec![9],
                vec![10]
            ]
        );
        m.atoms[0].query = Some(".R0".into());
        assert!(find_selected_matches(&m, &[0]).unwrap().matches.is_empty());
    }

    #[test]
    fn evaluates_ring_bond_count_and_exact_minimum_ring_size() {
        let mut m = Molecule::default();
        m.atoms.extend([
            atom("C", Some(".r3;.rb2")),
            atom("C", None),
            atom("C", None),
            atom("C", None),
        ]);
        m.bonds.extend([
            bond(1, 2, 1, None),
            bond(2, 3, 1, None),
            bond(3, 1, 1, None),
        ]);
        let report = find_selected_matches(&m, &[0]).unwrap();
        assert_eq!(report.matches, vec![vec![1], vec![2], vec![3]]);
    }

    #[test]
    fn bond_order_query_matches_only_allowed_orders() {
        let mut m = Molecule::default();
        m.atoms.extend((0..6).map(|_| atom("C", None)));
        m.atoms[0] = atom("C", None);
        m.atoms[1] = atom("C", None);
        m.bonds.extend([
            bond(0, 1, 1, Some("single-or-double")),
            bond(2, 3, 1, None),
            bond(4, 5, 3, None),
        ]);
        assert_eq!(
            find_selected_matches(&m, &[0, 1]).unwrap().matches,
            vec![vec![2, 3]]
        );
        m.bonds[0].query = Some("double-or-aromatic".into());
        assert!(find_selected_matches(&m, &[0, 1])
            .unwrap()
            .matches
            .is_empty());
    }

    #[test]
    fn ring_topology_bond_query_rejects_acyclic_edges() {
        let mut m = Molecule::default();
        m.atoms.extend((0..7).map(|_| atom("C", None)));
        m.bonds.extend([
            bond(0, 1, 1, Some("topology-ring")),
            bond(2, 3, 1, None),
            bond(3, 4, 1, None),
            bond(4, 2, 1, None),
            bond(5, 6, 1, None),
        ]);
        assert_eq!(
            find_selected_matches(&m, &[0, 1]).unwrap().matches,
            vec![vec![2, 3], vec![2, 4], vec![3, 4]]
        );
    }

    #[test]
    fn chain_topology_bond_query_rejects_ring_edges() {
        let mut m = Molecule::default();
        m.atoms.extend((0..7).map(|_| atom("C", None)));
        m.bonds.extend([
            bond(0, 1, 1, Some("topology-chain")),
            bond(2, 3, 1, None),
            bond(3, 4, 1, None),
            bond(4, 2, 1, None),
            bond(5, 6, 1, None),
        ]);
        assert_eq!(
            find_selected_matches(&m, &[0, 1]).unwrap().matches,
            vec![vec![5, 6]]
        );
    }

    #[test]
    fn aromatic_query_uses_aromatic_bond_evidence_even_without_atom_flag() {
        let mut m = Molecule::default();
        m.atoms.push(atom("C", Some(".a")));
        m.atoms.extend((0..6).map(|_| atom("C", None)));
        for i in 0..6 {
            let mut aromatic = bond(i + 1, (i + 1) % 6 + 1, 4, None);
            aromatic.style = BondStyle::Aromatic;
            m.bonds.push(aromatic);
        }
        assert_eq!(
            find_selected_matches(&m, &[0]).unwrap().matches,
            (1..=6).map(|index| vec![index]).collect::<Vec<_>>()
        );
    }

    #[test]
    fn evaluates_non_ring_atom_query_without_conflating_ring_system_count() {
        let mut m = Molecule::default();
        m.atoms
            .extend([atom("C", Some(".R0")), atom("C", None), atom("C", None)]);
        m.bonds.push(bond(1, 2, 1, None));
        assert_eq!(
            find_selected_matches(&m, &[0]).unwrap().matches,
            vec![vec![1], vec![2]]
        );
    }

    #[test]
    fn evaluates_exact_substitution_counts_from_heavy_and_isotopic_hydrogen_neighbors() {
        let mut m = Molecule::default();
        m.atoms.extend([
            atom("C", Some(".s3")), // Query atom
            atom("C", None),
            atom("C", None),
            atom("C", None), // Target central carbon
            atom("O", None),
            atom("N", None),
            atom("H", None), // Deuterium counts as a substituent
        ]);
        m.atoms[6].isotope = Some(2);
        m.bonds.extend([
            bond(3, 4, 1, None),
            bond(3, 5, 1, None),
            bond(3, 6, 1, None),
        ]);
        assert_eq!(
            find_selected_matches(&m, &[0]).unwrap().matches,
            vec![vec![3]]
        );
    }

    #[test]
    fn substitution_as_drawn_uses_selected_query_neighbors() {
        let mut m = Molecule::default();
        m.atoms.extend([
            atom("C", Some(".s*")),
            atom("C", None),
            atom("C", None),
            atom("C", None),
        ]);
        m.bonds.extend([bond(0, 1, 1, None), bond(2, 3, 1, None)]);
        assert_eq!(
            find_selected_matches(&m, &[0, 1]).unwrap().matches,
            vec![vec![2, 3]]
        );
    }

    #[test]
    fn refuses_a_partial_selection_that_cuts_through_a_query_molecule() {
        let mut m = Molecule::default();
        m.atoms
            .extend([atom("C", None), atom("C", None), atom("O", None)]);
        m.bonds.extend([bond(0, 1, 1, None), bond(1, 2, 1, None)]);
        assert!(find_selected_matches(&m, &[0, 1])
            .unwrap_err()
            .contains("complete query molecule"));
    }

    #[test]
    fn rejects_disconnected_query_selection_and_does_not_match_query_component_itself() {
        let mut m = Molecule::default();
        m.atoms.extend((0..3).map(|_| atom("C", None)));
        m.bonds.push(bond(0, 1, 1, None));
        assert!(find_selected_matches(&m, &[0, 2]).is_err());
        assert!(find_selected_matches(&m, &[0, 1])
            .unwrap()
            .matches
            .is_empty());
    }
}
