use crate::{chemistry, smarts, Atom, Bond, BondStyle, Molecule, BOND};
use macroquad::prelude::{vec2, vec3, Vec2};
use std::collections::{BTreeMap, BTreeSet};

const MAX_CML_BYTES: usize = 10 * 1024 * 1024;
const MAX_XML_NODES: usize = 100_000;
const MAX_XML_DEPTH: usize = 64;

#[derive(Debug, Default)]
struct Node {
    name: String,
    attributes: BTreeMap<String, String>,
    children: Vec<Node>,
    text: String,
}

impl Node {
    fn attribute(&self, name: &str) -> Option<&str> {
        self.attributes.get(name).map(String::as_str)
    }

    fn children_named<'a>(&'a self, name: &'a str) -> impl Iterator<Item = &'a Node> {
        self.children.iter().filter(move |child| child.name == name)
    }

    fn first_descendant(&self, name: &str) -> Option<&Node> {
        self.children.iter().find_map(|child| {
            if child.name == name {
                Some(child)
            } else {
                child.first_descendant(name)
            }
        })
    }
}

fn local_name(name: &str) -> &str {
    name.rsplit(':').next().unwrap_or(name)
}

fn decode_entities(value: &str) -> Result<String, String> {
    let mut output = String::with_capacity(value.len());
    let mut remaining = value;
    while let Some(index) = remaining.find('&') {
        output.push_str(&remaining[..index]);
        remaining = &remaining[index + 1..];
        let end = remaining
            .find(';')
            .ok_or("CML XML contains an unterminated entity reference")?;
        let entity = &remaining[..end];
        let decoded = match entity {
            "amp" => '&',
            "lt" => '<',
            "gt" => '>',
            "quot" => '"',
            "apos" => '\'',
            _ if entity.starts_with("#x") => {
                let code = u32::from_str_radix(&entity[2..], 16)
                    .map_err(|_| "CML XML contains an invalid hexadecimal character reference")?;
                char::from_u32(code)
                    .ok_or("CML XML contains an invalid Unicode character reference")?
            }
            _ if entity.starts_with('#') => {
                let code = entity[1..]
                    .parse::<u32>()
                    .map_err(|_| "CML XML contains an invalid decimal character reference")?;
                char::from_u32(code)
                    .ok_or("CML XML contains an invalid Unicode character reference")?
            }
            _ => return Err("CML XML contains an unsupported entity reference".into()),
        };
        output.push(decoded);
        remaining = &remaining[end + 1..];
    }
    output.push_str(remaining);
    Ok(output)
}

fn append_node(node: Node, stack: &mut [Node], roots: &mut Vec<Node>) {
    if let Some(parent) = stack.last_mut() {
        parent.children.push(node);
    } else {
        roots.push(node);
    }
}

fn skip_until(input: &str, position: &mut usize, end: &str) -> Result<(), String> {
    let rest = &input[*position..];
    let length = rest
        .find(end)
        .ok_or("CML XML contains an unterminated declaration")?;
    *position += length + end.len();
    Ok(())
}

fn read_xml(input: &str) -> Result<Node, String> {
    if input.len() > MAX_CML_BYTES {
        return Err("CML file exceeds the 10 MiB import limit".into());
    }
    let mut position = 0usize;
    let mut stack: Vec<Node> = Vec::new();
    let mut roots = Vec::new();
    let mut node_count = 0usize;
    while position < input.len() {
        let rest = &input[position..];
        let Some(open) = rest.find('<') else {
            let text = decode_entities(rest)?;
            if !text.trim().is_empty() && stack.is_empty() {
                return Err("CML XML has text outside its root element".into());
            }
            if let Some(node) = stack.last_mut() {
                node.text.push_str(&text);
            }
            break;
        };
        if open > 0 {
            let text = decode_entities(&rest[..open])?;
            if let Some(node) = stack.last_mut() {
                node.text.push_str(&text);
            } else if !text.trim().is_empty() {
                return Err("CML XML has text outside its root element".into());
            }
            position += open;
        }
        if input[position..].starts_with("<!--") {
            position += 4;
            skip_until(input, &mut position, "-->")?;
            continue;
        }
        if input[position..].starts_with("<?") {
            position += 2;
            skip_until(input, &mut position, "?>")?;
            continue;
        }
        if input[position..].starts_with("<![CDATA[") {
            position += 9;
            let rest = &input[position..];
            let end = rest
                .find("]]>")
                .ok_or("CML XML contains an unterminated CDATA section")?;
            if let Some(node) = stack.last_mut() {
                node.text.push_str(&rest[..end]);
            } else if !rest[..end].trim().is_empty() {
                return Err("CML XML has text outside its root element".into());
            }
            position += end + 3;
            continue;
        }
        if input[position..].starts_with("<!") {
            return Err("CML XML declarations such as DOCTYPE are not supported".into());
        }
        position += 1;
        if input[position..].starts_with('/') {
            position += 1;
            let end = input[position..]
                .find('>')
                .ok_or("CML XML has an unterminated closing tag")?
                + position;
            let closing = local_name(input[position..end].trim());
            let node = stack
                .pop()
                .ok_or("CML XML contains a closing tag without an opening tag")?;
            if node.name != closing {
                return Err(format!(
                    "CML XML closing tag </{closing}> does not match <{}>",
                    node.name
                ));
            }
            append_node(node, &mut stack, &mut roots);
            position = end + 1;
            continue;
        }

        let mut cursor = position;
        let bytes = input.as_bytes();
        while cursor < input.len()
            && !bytes[cursor].is_ascii_whitespace()
            && !matches!(bytes[cursor], b'/' | b'>')
        {
            cursor += 1;
        }
        if cursor == position {
            return Err("CML XML contains an empty element name".into());
        }
        let name = local_name(&input[position..cursor]).to_string();
        let mut node = Node {
            name,
            ..Node::default()
        };
        let mut self_closing = false;
        loop {
            while cursor < input.len() && bytes[cursor].is_ascii_whitespace() {
                cursor += 1;
            }
            if cursor >= input.len() {
                return Err("CML XML contains an unterminated start tag".into());
            }
            match bytes[cursor] {
                b'>' => {
                    cursor += 1;
                    break;
                }
                b'/' if bytes.get(cursor + 1) == Some(&b'>') => {
                    cursor += 2;
                    self_closing = true;
                    break;
                }
                _ => {}
            }
            let attr_start = cursor;
            while cursor < input.len()
                && !bytes[cursor].is_ascii_whitespace()
                && !matches!(bytes[cursor], b'=' | b'/' | b'>')
            {
                cursor += 1;
            }
            if cursor == attr_start {
                return Err("CML XML contains an invalid attribute name".into());
            }
            let raw_name = &input[attr_start..cursor];
            let name = local_name(raw_name).to_string();
            while cursor < input.len() && bytes[cursor].is_ascii_whitespace() {
                cursor += 1;
            }
            if bytes.get(cursor) != Some(&b'=') {
                return Err(format!("CML XML attribute {raw_name} has no value"));
            }
            cursor += 1;
            while cursor < input.len() && bytes[cursor].is_ascii_whitespace() {
                cursor += 1;
            }
            let quote = *bytes
                .get(cursor)
                .ok_or("CML XML has an incomplete attribute value")?;
            if quote != b'\'' && quote != b'"' {
                return Err("CML XML attribute values must be quoted".into());
            }
            cursor += 1;
            let value_start = cursor;
            while cursor < input.len() && bytes[cursor] != quote {
                cursor += 1;
            }
            if cursor >= input.len() {
                return Err("CML XML has an unterminated attribute value".into());
            }
            let value = decode_entities(&input[value_start..cursor])?;
            cursor += 1;
            if !raw_name.starts_with("xmlns")
                && node.attributes.insert(name.clone(), value).is_some()
            {
                return Err(format!("CML XML repeats attribute {name}"));
            }
        }
        node_count += 1;
        if node_count > MAX_XML_NODES {
            return Err("CML XML exceeds the 100,000-element import limit".into());
        }
        if self_closing {
            append_node(node, &mut stack, &mut roots);
        } else {
            if stack.len() >= MAX_XML_DEPTH {
                return Err("CML XML exceeds the 64-level nesting limit".into());
            }
            stack.push(node);
        }
        position = cursor;
    }
    if !stack.is_empty() {
        return Err("CML XML has an unclosed element".into());
    }
    if roots.len() != 1 {
        return Err("CML XML must contain exactly one root element".into());
    }
    Ok(roots.remove(0))
}

fn parse_number(node: &Node, key: &str) -> Result<Option<f32>, String> {
    node.attribute(key)
        .map(|value| {
            value
                .parse::<f32>()
                .ok()
                .filter(|number| number.is_finite())
                .ok_or_else(|| format!("CML attribute {key} must be a finite number"))
        })
        .transpose()
}

fn parse_small_integer<T>(node: &Node, key: &str) -> Result<Option<T>, String>
where
    T: std::str::FromStr,
{
    node.attribute(key)
        .map(|value| {
            value
                .parse::<T>()
                .map_err(|_| format!("CML attribute {key} has an invalid integer value"))
        })
        .transpose()
}

fn find_molecules<'a>(node: &'a Node, out: &mut Vec<&'a Node>) {
    if node.name == "molecule" {
        out.push(node);
    } else {
        for child in &node.children {
            find_molecules(child, out);
        }
    }
}

pub fn parse(input: &str, center: Vec2) -> Result<Molecule, String> {
    let root = read_xml(input)?;
    if root.name != "cml" && root.name != "molecule" {
        return Err("CML root must be <cml> or <molecule>".into());
    }
    let mut molecules = Vec::new();
    find_molecules(&root, &mut molecules);
    if molecules.len() != 1 {
        return Err(if molecules.is_empty() {
            "CML file contains no molecule".into()
        } else {
            "CML import currently accepts one molecule per file; split the moleculeList first"
                .into()
        });
    }
    let molecule = molecules[0];
    if molecule.children.iter().any(|child| {
        !matches!(
            child.name.as_str(),
            "name" | "atomArray" | "bondArray" | "metadataList"
        )
    }) {
        return Err("CML molecule has annotations outside the supported CML subset".into());
    }
    let atom_array = molecule.first_descendant("atomArray");
    let bond_array = molecule.first_descendant("bondArray");
    let atom_nodes: Vec<&Node> = atom_array
        .map(|array| array.children_named("atom").collect())
        .unwrap_or_default();
    if atom_nodes.is_empty() {
        return Err("CML molecule has no atoms in its atomArray".into());
    }
    if atom_nodes.len() > 100_000 {
        return Err("CML molecule exceeds the 100,000-atom import limit".into());
    }
    let mut result = Molecule {
        name: molecule
            .attribute("title")
            .map(str::to_string)
            .or_else(|| {
                molecule
                    .children_named("name")
                    .next()
                    .map(|name| name.text.trim().to_string())
            })
            .unwrap_or_default(),
        ..Molecule::default()
    };
    let mut atom_ids = BTreeMap::new();
    let mut raw_coordinates = Vec::with_capacity(atom_nodes.len());
    for node in atom_nodes {
        let id = node
            .attribute("id")
            .ok_or("CML atom is missing its required id attribute")?;
        if id.trim().is_empty() || atom_ids.contains_key(id) {
            return Err(format!("CML atom id {id:?} is empty or duplicated"));
        }
        if !node.children.is_empty() {
            return Err(format!(
                "CML atom {id} has child annotations outside the supported CML subset"
            ));
        }
        let element = node
            .attribute("elementType")
            .ok_or_else(|| format!("CML atom {id} has no elementType"))?;
        if !chemistry::valid_element(element) {
            return Err(format!(
                "CML atom {id} uses unsupported element {element:?}"
            ));
        }
        if node.attribute("atomParity").is_some() || node.attribute("stereo").is_some() {
            return Err(format!(
                "CML atom {id} has stereochemistry not supported by this importer"
            ));
        }
        if node.attribute("dictRef").is_some() || node.attribute("aromatic").is_some() {
            return Err(format!(
                "CML atom {id} has query or aromatic annotations outside the supported CML subset"
            ));
        }
        let charge = parse_small_integer::<i8>(node, "formalCharge")?.unwrap_or(0);
        if !(-9..=9).contains(&charge) {
            return Err(format!("CML atom {id} formalCharge is outside -9..9"));
        }
        let isotope = parse_small_integer::<u16>(node, "isotopeNumber")?;
        if isotope == Some(0) {
            return Err(format!("CML atom {id} isotopeNumber must be positive"));
        }
        let hydrogen_count = parse_small_integer::<u8>(node, "hydrogenCount")?;
        if hydrogen_count.is_some_and(|count| count > 5) {
            return Err(format!(
                "CML atom {id} hydrogenCount above 5 is not supported"
            ));
        }
        let x2 = parse_number(node, "x2")?;
        let y2 = parse_number(node, "y2")?;
        if x2.is_some() != y2.is_some() {
            return Err(format!(
                "CML atom {id} must provide both x2 and y2 coordinates"
            ));
        }
        let x3 = parse_number(node, "x3")?;
        let y3 = parse_number(node, "y3")?;
        let z3 = parse_number(node, "z3")?;
        let coord3d = match (x3, y3, z3) {
            (Some(x), Some(y), Some(z)) => Some(vec3(x, y, z)),
            (None, None, None) => None,
            _ => return Err(format!("CML atom {id} must provide x3, y3 and z3 together")),
        };
        let position = match (x2, y2) {
            (Some(x), Some(y)) => vec2(center.x + x * (BOND / 1.5), center.y - y * (BOND / 1.5)),
            _ => match (x3, y3) {
                (Some(x), Some(y)) => {
                    vec2(center.x + x * (BOND / 1.5), center.y - y * (BOND / 1.5))
                }
                _ => Vec2::ZERO,
            },
        };
        let mut atom = Atom::new(position, element);
        atom.charge = charge;
        atom.isotope = isotope;
        atom.implicit_h_override = hydrogen_count;
        atom.coord3d = coord3d;
        atom.z = z3;
        atom_ids.insert(id.to_string(), result.atoms.len());
        raw_coordinates.push(x2.is_some() || x3.is_some());
        result.atoms.push(atom);
    }

    let bond_nodes: Vec<&Node> = bond_array
        .map(|array| array.children_named("bond").collect())
        .unwrap_or_default();
    if bond_nodes.len() > 200_000 {
        return Err("CML molecule exceeds the 200,000-bond import limit".into());
    }
    let mut seen_bonds = BTreeSet::new();
    for node in bond_nodes {
        if !node.children.is_empty() {
            return Err("CML bond has child annotations outside the supported CML subset".into());
        }
        let refs: Vec<_> = node
            .attribute("atomRefs2")
            .ok_or("CML bond is missing atomRefs2")?
            .split_whitespace()
            .collect();
        if refs.len() != 2 {
            return Err("CML bond atomRefs2 must contain exactly two atom ids".into());
        }
        let a = *atom_ids
            .get(refs[0])
            .ok_or_else(|| format!("CML bond refers to unknown atom {}", refs[0]))?;
        let b = *atom_ids
            .get(refs[1])
            .ok_or_else(|| format!("CML bond refers to unknown atom {}", refs[1]))?;
        if a == b {
            return Err("CML bond cannot connect an atom to itself".into());
        }
        let key = (a.min(b), a.max(b));
        if !seen_bonds.insert(key) {
            return Err("CML molecule contains a duplicate bond".into());
        }
        let order_text = node.attribute("order").unwrap_or("1");
        let order = match order_text {
            "1" => 1,
            "2" => 2,
            "3" => 3,
            "A" | "1.5" => 4,
            _ => return Err(format!("CML bond order {order_text:?} is not supported")),
        };
        let stereo = node.attribute("stereo").map(str::to_string);
        if stereo
            .as_deref()
            .is_some_and(|s| !matches!(s, "E" | "Z" | "cis" | "trans"))
        {
            return Err("CML bond stereo must be E, Z, cis or trans".into());
        }
        if stereo.is_some() && order != 2 {
            return Err("CML bond stereo can only be assigned to a double bond".into());
        }
        result.bonds.push(Bond {
            a,
            b,
            order,
            style: if order == 4 {
                BondStyle::Aromatic
            } else {
                BondStyle::Normal
            },
            stereo,
            query: None,
        });
    }
    for bond in result.bonds.iter().filter(|bond| bond.order == 4) {
        result.atoms[bond.a].aromatic = true;
        result.atoms[bond.b].aromatic = true;
    }
    if raw_coordinates.iter().any(|has| *has) && !raw_coordinates.iter().all(|has| *has) {
        return Err(
            "CML molecule has coordinates for only some atoms; refusing to discard them".into(),
        );
    }
    if raw_coordinates.iter().all(|has| !*has) {
        smarts::layout_imported_molecule(&mut result);
        for atom in &mut result.atoms {
            atom.pos += center;
        }
    }
    for metadata in molecule
        .first_descendant("metadataList")
        .into_iter()
        .flat_map(|list| list.children_named("metadata"))
    {
        if let (Some(name), Some(value)) =
            (metadata.attribute("name"), metadata.attribute("content"))
        {
            result
                .properties
                .push((name.to_string(), value.to_string()));
        }
    }
    Ok(result)
}

fn xml_escape(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

pub fn write(molecule: &Molecule) -> String {
    let scale = BOND / 1.5;
    let mut output = String::from(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<cml xmlns=\"http://www.xml-cml.org/schema\">\n",
    );
    output.push_str(&format!(
        "  <molecule id=\"m1\" title=\"{}\">\n",
        xml_escape(&molecule.name)
    ));
    output.push_str("    <atomArray>\n");
    for (index, atom) in molecule.atoms.iter().enumerate() {
        output.push_str(&format!(
            "      <atom id=\"a{}\" elementType=\"{}\" x2=\"{:.5}\" y2=\"{:.5}\"",
            index + 1,
            xml_escape(&atom.element),
            (atom.pos.x / scale),
            (-atom.pos.y / scale)
        ));
        if atom.charge != 0 {
            output.push_str(&format!(" formalCharge=\"{}\"", atom.charge));
        }
        if let Some(isotope) = atom.isotope {
            output.push_str(&format!(" isotopeNumber=\"{isotope}\""));
        }
        if let Some(hydrogens) = atom.implicit_h_override {
            output.push_str(&format!(" hydrogenCount=\"{hydrogens}\""));
        }
        if let Some(point) = atom.coord3d {
            output.push_str(&format!(
                " x3=\"{:.5}\" y3=\"{:.5}\" z3=\"{:.5}\"",
                point.x, point.y, point.z
            ));
        }
        output.push_str("/>\n");
    }
    output.push_str("    </atomArray>\n    <bondArray>\n");
    for (index, bond) in molecule.bonds.iter().enumerate() {
        let order = match bond.order {
            1 => "1",
            2 => "2",
            3 => "3",
            4 => "A",
            _ => "1",
        };
        output.push_str(&format!(
            "      <bond id=\"b{}\" atomRefs2=\"a{} a{}\" order=\"{}\"",
            index + 1,
            bond.a + 1,
            bond.b + 1,
            order
        ));
        if let Some(stereo) = &bond.stereo {
            output.push_str(&format!(" stereo=\"{}\"", xml_escape(stereo)));
        }
        output.push_str("/>\n");
    }
    output.push_str("    </bondArray>\n");
    if !molecule.properties.is_empty() {
        output.push_str("    <metadataList>");
        for (name, value) in &molecule.properties {
            output.push_str(&format!(
                "<metadata name=\"{}\" content=\"{}\"/>",
                xml_escape(name),
                xml_escape(value)
            ));
        }
        output.push_str("</metadataList>\n");
    }
    output.push_str("  </molecule>\n</cml>\n");
    output
}

pub fn export_losses(molecule: &Molecule, has_graphics: bool) -> Vec<&'static str> {
    let mut losses = Vec::new();
    if has_graphics {
        losses.push("graphics, reaction arrows and text");
    }
    if molecule.atoms.iter().enumerate().any(|(index, atom)| {
        !chemistry::valid_element(&atom.element)
            || !atom.pos.x.is_finite()
            || !atom.pos.y.is_finite()
            || atom
                .coord3d
                .is_some_and(|p| !p.x.is_finite() || !p.y.is_finite() || !p.z.is_finite())
            || atom.alias.is_some()
            || atom.r_group.is_some()
            || atom.query.is_some()
            || atom.radical != 0
            || atom.valence_override.is_some()
            || atom.stereo.is_some()
            || atom.atom_map.is_some()
            || atom.explicit_h > 0
            || atom.implicit_h_override.is_some_and(|count| count > 5)
            || (atom.z.is_some() && atom.coord3d.is_none())
            || (atom.aromatic
                && !molecule
                    .bonds
                    .iter()
                    .any(|bond| bond.order == 4 && (bond.a == index || bond.b == index)))
    }) {
        losses.push("atom labels, queries, radicals, stereochemistry, maps or valence overrides");
    }
    let mut seen_bonds = BTreeSet::new();
    if molecule.bonds.iter().any(|bond| {
        let key = (bond.a.min(bond.b), bond.a.max(bond.b));
        bond.a >= molecule.atoms.len()
            || bond.b >= molecule.atoms.len()
            || bond.a == bond.b
            || !seen_bonds.insert(key)
            || !matches!(bond.order, 1..=4)
            || bond.query.is_some()
            || !matches!(bond.style, BondStyle::Normal | BondStyle::Aromatic)
            || bond
                .stereo
                .as_deref()
                .is_some_and(|s| !matches!(s, "E" | "Z" | "cis" | "trans"))
    }) {
        losses.push("bond query, display styles or unsupported stereo labels");
    }
    if molecule.atoms.iter().any(|atom| {
        atom.reaction_role.is_some() || atom.reaction_step != 1 || atom.coefficient != 1
    }) {
        losses.push("reaction roles, steps or stoichiometric coefficients");
    }
    losses
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cml_roundtrip_preserves_supported_atom_bond_and_metadata_fields() {
        let mut molecule = Molecule {
            name: "ethanol & sample".into(),
            ..Molecule::default()
        };
        let mut carbon = Atom::new(vec2(0., 0.), "C");
        carbon.implicit_h_override = Some(3);
        molecule.atoms.push(carbon);
        let mut oxygen = Atom::new(vec2(BOND * 2., -BOND), "O");
        oxygen.charge = -1;
        oxygen.isotope = Some(18);
        oxygen.coord3d = Some(vec3(1.0, 2.0, 3.0));
        oxygen.z = Some(3.0);
        molecule.atoms.push(oxygen);
        molecule.bonds.push(Bond {
            a: 0,
            b: 1,
            order: 1,
            style: BondStyle::Normal,
            stereo: None,
            query: None,
        });
        molecule
            .properties
            .push(("SOURCE".into(), "lab <A> & notes".into()));

        let loaded = parse(&write(&molecule), vec2(100., 200.)).unwrap();
        assert_eq!(loaded.name, molecule.name);
        assert_eq!(loaded.atoms.len(), 2);
        assert_eq!(loaded.bonds, molecule.bonds);
        assert_eq!(loaded.atoms[0].implicit_h_override, Some(3));
        assert_eq!(loaded.atoms[1].charge, -1);
        assert_eq!(loaded.atoms[1].isotope, Some(18));
        assert_eq!(loaded.atoms[1].coord3d, Some(vec3(1.0, 2.0, 3.0)));
        assert_eq!(loaded.properties, molecule.properties);
    }

    #[test]
    fn cml_parser_handles_namespaces_aromaticity_and_graph_layout() {
        let input = r#"<?xml version="1.0"?><cml xmlns="http://www.xml-cml.org/schema"><molecule title="benzene"><atomArray><atom id="a1" elementType="C"/><atom id="a2" elementType="C"/><atom id="a3" elementType="C"/><atom id="a4" elementType="C"/><atom id="a5" elementType="C"/><atom id="a6" elementType="C"/></atomArray><bondArray><bond atomRefs2="a1 a2" order="A"/><bond atomRefs2="a2 a3" order="A"/><bond atomRefs2="a3 a4" order="A"/><bond atomRefs2="a4 a5" order="A"/><bond atomRefs2="a5 a6" order="A"/><bond atomRefs2="a6 a1" order="A"/></bondArray></molecule></cml>"#;
        let molecule = parse(input, vec2(300., 300.)).unwrap();
        assert_eq!(molecule.name, "benzene");
        assert_eq!(molecule.atoms.len(), 6);
        assert_eq!(molecule.bonds.len(), 6);
        assert!(molecule.bonds.iter().all(|bond| bond.order == 4));
        assert!(molecule.atoms.iter().all(|atom| atom.aromatic));
        assert!(molecule.atoms.iter().all(|atom| atom.pos.is_finite()));
        let partial = r#"<molecule><atomArray><atom id="a1" elementType="C" x2="0" y2="0"/><atom id="a2" elementType="C"/></atomArray><bondArray><bond atomRefs2="a1 a2"/></bondArray></molecule>"#;
        assert!(parse(partial, vec2(0., 0.))
            .unwrap_err()
            .contains("only some atoms"));
    }

    #[test]
    fn cml_parser_rejects_bad_references_and_unsupported_stereo_without_partial_graph() {
        let bad_reference = r#"<molecule><atomArray><atom id="a1" elementType="C"/></atomArray><bondArray><bond atomRefs2="a1 missing"/></bondArray></molecule>"#;
        assert!(parse(bad_reference, Vec2::ZERO)
            .unwrap_err()
            .contains("unknown atom"));
        let unsupported = r#"<molecule><atomArray><atom id="a1" elementType="C" atomParity="1"/></atomArray></molecule>"#;
        assert!(parse(unsupported, Vec2::ZERO)
            .unwrap_err()
            .contains("stereochemistry"));
        let child_stereo = r#"<molecule><atomArray><atom id="a1" elementType="C"/><atom id="a2" elementType="C"/></atomArray><bondArray><bond atomRefs2="a1 a2"><bondStereo>W</bondStereo></bond></bondArray></molecule>"#;
        assert!(parse(child_stereo, Vec2::ZERO)
            .unwrap_err()
            .contains("child annotations"));
    }

    #[test]
    fn cml_xml_rejects_malformed_markup_and_doctype() {
        assert!(read_xml("<cml><molecule></cml>").is_err());
        assert!(read_xml("<!DOCTYPE cml><cml/>").is_err());
        assert_eq!(decode_entities("A &amp; B &#x3C; C").unwrap(), "A & B < C");
    }

    #[test]
    fn cml_export_reports_fields_outside_its_supported_interchange_subset() {
        let mut molecule = Molecule::default();
        let mut atom = Atom::new(Vec2::ZERO, "C");
        atom.atom_map = Some(12);
        atom.query = Some(".a".into());
        molecule.atoms.push(atom);
        assert_eq!(
            export_losses(&molecule, false),
            vec!["atom labels, queries, radicals, stereochemistry, maps or valence overrides"]
        );
        assert_eq!(
            export_losses(&Molecule::default(), true),
            vec!["graphics, reaction arrows and text"]
        );
        molecule.bonds.push(Bond {
            a: 0,
            b: 0,
            order: 7,
            style: BondStyle::Normal,
            stereo: None,
            query: None,
        });
        assert!(export_losses(&molecule, false)
            .iter()
            .any(|loss| loss.contains("bond query")));
    }
}
