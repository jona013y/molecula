//! Optional direct PubChem lookup through the official PUG REST API.
//! Results are database records; this module does not calculate properties.

use std::process::Command;

#[cfg(test)]
use std::collections::VecDeque;

const CID_ENDPOINT: &str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastidentity/smiles/cids/TXT?identity_type=same_stereo_isotope&MaxRecords=10";
const PROPERTY_NAMES: &str = "MolecularFormula,MolecularWeight,InChIKey,IUPACName,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount";

fn curl(args: &[&str]) -> Result<String, String> {
    let output = Command::new("curl")
        .args([
            "--fail",
            "--silent",
            "--show-error",
            "--connect-timeout",
            "8",
            "--max-time",
            "30",
            "--max-filesize",
            "1048576",
            "--user-agent",
            "Molecula/1.0.12 (https://github.com/jona013y/molecula)",
        ])
        .args(args)
        .output()
        .map_err(|error| format!("Could not start curl: {error}"))?;
    if !output.status.success() {
        let detail = String::from_utf8_lossy(&output.stderr);
        return Err(format!(
            "PubChem request failed: {}",
            detail.trim().chars().take(400).collect::<String>()
        ));
    }
    String::from_utf8(output.stdout).map_err(|_| "PubChem returned non-UTF-8 data".into())
}

trait PubChemClient {
    fn exact_identity_cids(&mut self, smiles: &str) -> Result<String, String>;
    fn properties_csv(&mut self, cids: &str) -> Result<String, String>;
}

struct CurlPubChemClient;

impl PubChemClient for CurlPubChemClient {
    fn exact_identity_cids(&mut self, smiles: &str) -> Result<String, String> {
        let form = format!("smiles={smiles}");
        curl(&["--request", "POST", "--data-urlencode", &form, CID_ENDPOINT])
    }

    fn properties_csv(&mut self, cids: &str) -> Result<String, String> {
        let url = format!(
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}/property/{PROPERTY_NAMES}/CSV"
        );
        curl(&[&url])
    }
}

fn cids_from_text(body: &str) -> Result<Vec<u64>, String> {
    let mut cids = Vec::new();
    for line in body.lines().map(str::trim).filter(|line| !line.is_empty()) {
        let cid = line
            .parse::<u64>()
            .map_err(|_| "PubChem returned an invalid compound identifier".to_string())?;
        if cid == 0 {
            return Err("PubChem returned an invalid zero compound identifier".into());
        }
        if !cids.contains(&cid) {
            cids.push(cid);
        }
        if cids.len() > 10 {
            return Err("PubChem returned more results than the requested limit".into());
        }
    }
    Ok(cids)
}

/// Parse RFC 4180-style CSV, including quoted commas, escaped quotes and newlines.
fn parse_csv(input: &str) -> Result<Vec<Vec<String>>, String> {
    let mut rows = Vec::new();
    let mut row = Vec::new();
    let mut field = String::new();
    let mut chars = input.chars().peekable();
    let mut quoted = false;
    let mut closed_quote = false;

    while let Some(ch) = chars.next() {
        if quoted {
            match ch {
                '"' if chars.peek() == Some(&'"') => {
                    chars.next();
                    field.push('"');
                }
                '"' => {
                    quoted = false;
                    closed_quote = true;
                }
                other => field.push(other),
            }
            continue;
        }
        if closed_quote && !matches!(ch, ',' | '\r' | '\n') {
            return Err("PubChem returned malformed CSV after a quoted field".into());
        }
        match ch {
            '"' if field.is_empty() && !closed_quote => quoted = true,
            ',' => {
                row.push(std::mem::take(&mut field));
                closed_quote = false;
            }
            '\r' | '\n' => {
                if ch == '\r' && chars.peek() == Some(&'\n') {
                    chars.next();
                }
                row.push(std::mem::take(&mut field));
                rows.push(std::mem::take(&mut row));
                closed_quote = false;
            }
            other => field.push(other),
        }
    }
    if quoted {
        return Err("PubChem returned an unterminated CSV field".into());
    }
    if !field.is_empty() || !row.is_empty() || closed_quote {
        row.push(field);
        rows.push(row);
    }
    Ok(rows)
}

fn format_properties(csv: &str) -> Result<String, String> {
    let rows = parse_csv(csv)?;
    let Some(header) = rows.first() else {
        return Err("PubChem returned an empty property table".into());
    };
    let required = [
        "CID",
        "MolecularFormula",
        "MolecularWeight",
        "InChIKey",
        "IUPACName",
        "XLogP",
        "TPSA",
        "HBondDonorCount",
        "HBondAcceptorCount",
    ];
    let positions: Vec<usize> = required
        .iter()
        .map(|name| {
            header
                .iter()
                .position(|item| item == name)
                .ok_or_else(|| format!("PubChem omitted the requested `{name}` column"))
        })
        .collect::<Result<_, _>>()?;
    let mut report = String::from(
        "PubChem database records (direct PUG REST lookup; not calculated by Molecula)\n\n",
    );
    let mut count = 0usize;
    for row in rows.iter().skip(1) {
        if row.iter().all(String::is_empty) {
            continue;
        }
        if row.len() != header.len() {
            return Err("PubChem returned a malformed property row".into());
        }
        let get = |column: usize| row[positions[column]].trim();
        let cid = get(0)
            .parse::<u64>()
            .map_err(|_| "PubChem returned an invalid compound identifier".to_string())?;
        if cid == 0 {
            return Err("PubChem returned an invalid zero compound identifier".into());
        }
        count += 1;
        report.push_str(&format!(
            "CID {cid} · https://pubchem.ncbi.nlm.nih.gov/compound/{cid}\n"
        ));
        for (label, column, unit) in [
            ("Formula", 1, ""),
            ("Molecular weight", 2, "g/mol"),
            ("InChIKey", 3, ""),
            ("PubChem IUPAC name", 4, ""),
            ("PubChem XLogP", 5, ""),
            ("PubChem TPSA", 6, "Å²"),
            ("H-bond donors", 7, ""),
            ("H-bond acceptors", 8, ""),
        ] {
            let value = get(column);
            report.push_str(&format!(
                "{label}: {}{}\n",
                if value.is_empty() {
                    "Not listed"
                } else {
                    value
                },
                if !unit.is_empty() && !value.is_empty() {
                    format!(" {unit}")
                } else {
                    String::new()
                }
            ));
        }
        report.push('\n');
    }
    if count == 0 {
        return Err("PubChem found no exact structure match (stereochemistry and isotopes included). The drawing was not changed.".into());
    }
    if count == 10 {
        report.push_str("Showing up to 10 exact matches.\n");
    }
    Ok(report)
}

/// Query PubChem for an exact stereochemical and isotopic structure identity.
fn lookup_with(client: &mut impl PubChemClient, smiles: &str) -> Result<String, String> {
    if smiles.is_empty() || smiles.len() > 16_384 {
        return Err("The exported structure string is empty or too long for PubChem lookup".into());
    }
    let response = client.exact_identity_cids(smiles)?;
    let cids = cids_from_text(&response)?;
    if cids.is_empty() {
        return Err("PubChem found no exact structure match (stereochemistry and isotopes included). The drawing was not changed.".into());
    }
    let cid_list = cids
        .iter()
        .map(u64::to_string)
        .collect::<Vec<_>>()
        .join(",");
    let csv = client.properties_csv(&cid_list)?;
    format_properties(&csv)
}

pub fn lookup(smiles: &str) -> Result<String, String> {
    lookup_with(&mut CurlPubChemClient, smiles)
}

#[cfg(test)]
mod tests {
    use super::*;

    struct FixturePubChemClient {
        responses: VecDeque<Result<String, String>>,
        requests: Vec<String>,
    }

    impl PubChemClient for FixturePubChemClient {
        fn exact_identity_cids(&mut self, smiles: &str) -> Result<String, String> {
            self.requests
                .push(format!("POST {CID_ENDPOINT}\nsmiles={smiles}"));
            self.responses.pop_front().expect("identity fixture")
        }

        fn properties_csv(&mut self, cids: &str) -> Result<String, String> {
            self.requests.push(format!(
                "GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cids}/property/{PROPERTY_NAMES}/CSV"
            ));
            self.responses.pop_front().expect("property fixture")
        }
    }

    #[test]
    fn parses_and_deduplicates_bounded_pubchem_cid_text() {
        assert_eq!(
            cids_from_text("2244\r\n2244\n1983\n").unwrap(),
            [2244, 1983]
        );
        assert!(cids_from_text("0").is_err());
        assert!(cids_from_text("not-a-cid").is_err());
        assert!(cids_from_text(&(1..=11).map(|id| format!("{id}\n")).collect::<String>()).is_err());
    }

    #[test]
    fn csv_parser_preserves_quoted_commas_quotes_and_multiline_values() {
        assert_eq!(
            parse_csv("CID,Name,Note\r\n2244,\"aspirin, acetyl\",\"said \"\"yes\"\"\nnext\"\r\n")
                .unwrap(),
            vec![
                vec!["CID", "Name", "Note"],
                vec!["2244", "aspirin, acetyl", "said \"yes\"\nnext"]
            ]
        );
        assert!(parse_csv("CID,Name\n1,\"unfinished").is_err());
        assert!(parse_csv("CID,Name\n1,\"ok\"oops").is_err());
    }

    #[test]
    fn formats_pubchem_rows_with_provenance_and_units() {
        let csv = "CID,MolecularFormula,MolecularWeight,InChIKey,IUPACName,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount\n2244,C9H8O4,180.16,BSYNRYMUTXBXSQ-UHFFFAOYSA-N,2-acetyloxybenzoic acid,1.2,63.6,1,4\n";
        let report = format_properties(csv).unwrap();
        assert!(report.contains("not calculated by Molecula"));
        assert!(report.contains("CID 2244"));
        assert!(report.contains("PubChem IUPAC name: 2-acetyloxybenzoic acid"));
        assert!(report.contains("PubChem TPSA: 63.6 Å²"));
        assert!(format_properties("CID,MolecularFormula\n1,C\n").is_err());
    }

    #[test]
    fn lookup_uses_exact_identity_query_before_fetching_properties() {
        let smiles = "[13CH3]CO";
        let csv = "CID,MolecularFormula,MolecularWeight,InChIKey,IUPACName,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount\n123,C2H6O,47.07,TEST-KEY,example,0.4,20.2,1,1\n";
        let mut client = FixturePubChemClient {
            responses: VecDeque::from([Ok("123\n".into()), Ok(csv.into())]),
            requests: Vec::new(),
        };

        let report = lookup_with(&mut client, smiles).unwrap();

        assert!(CID_ENDPOINT.contains("identity_type=same_stereo_isotope"));
        assert!(report.contains("CID 123"));
        assert_eq!(client.requests.len(), 2);
        assert_eq!(
            client.requests[0],
            format!("POST {CID_ENDPOINT}\nsmiles={smiles}")
        );
        assert_eq!(
            client.requests[1],
            format!(
                "GET https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/123/property/{PROPERTY_NAMES}/CSV"
            )
        );
    }

    #[test]
    fn exact_identity_miss_does_not_request_properties() {
        let mut client = FixturePubChemClient {
            responses: VecDeque::from([Ok(String::new())]),
            requests: Vec::new(),
        };
        let error = lookup_with(&mut client, "CCO").unwrap_err();
        assert!(error.contains("no exact structure match"));
        assert_eq!(client.requests.len(), 1);
        assert!(client.requests[0].contains("same_stereo_isotope"));
    }

    #[test]
    fn invalid_inputs_are_refused_before_contacting_pubchem() {
        let mut client = FixturePubChemClient {
            responses: VecDeque::new(),
            requests: Vec::new(),
        };
        assert!(lookup_with(&mut client, "").is_err());
        assert!(lookup_with(&mut client, &"C".repeat(16_385)).is_err());
        assert!(client.requests.is_empty());
    }
}
