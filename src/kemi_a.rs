//! Unit-aware, idealized calculations for common Chemistry A exercises.
//! Every result includes the equation, substituted values, and units.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BasicCalculation {
    AmountFromMass,
    MassFromAmount,
    Concentration,
    Dilution,
    PhFromHydronium,
    HydroniumFromPh,
    PohFromHydroxide,
    HydroxideFromPoh,
    PhFromPoh,
    PohFromPh,
}

impl BasicCalculation {
    pub const ALL: [Self; 10] = [
        Self::AmountFromMass,
        Self::MassFromAmount,
        Self::Concentration,
        Self::Dilution,
        Self::PhFromHydronium,
        Self::HydroniumFromPh,
        Self::PohFromHydroxide,
        Self::HydroxideFromPoh,
        Self::PhFromPoh,
        Self::PohFromPh,
    ];

    pub fn title(self, danish: bool) -> &'static str {
        match (self, danish) {
            (Self::AmountFromMass, true) => "Stofmængde fra masse",
            (Self::MassFromAmount, true) => "Masse fra stofmængde",
            (Self::Concentration, true) => "Koncentration",
            (Self::Dilution, true) => "Fortynding",
            (Self::PhFromHydronium, true) => "pH fra [H₃O⁺]",
            (Self::HydroniumFromPh, true) => "[H₃O⁺] fra pH",
            (Self::PohFromHydroxide, true) => "pOH fra [OH⁻]",
            (Self::HydroxideFromPoh, true) => "[OH⁻] fra pOH",
            (Self::PhFromPoh, true) => "pH fra pOH",
            (Self::PohFromPh, true) => "pOH fra pH",
            (Self::AmountFromMass, false) => "Amount of substance from mass",
            (Self::MassFromAmount, false) => "Mass from amount of substance",
            (Self::Concentration, false) => "Concentration",
            (Self::Dilution, false) => "Dilution",
            (Self::PhFromHydronium, false) => "pH from [H₃O⁺]",
            (Self::HydroniumFromPh, false) => "[H₃O⁺] from pH",
            (Self::PohFromHydroxide, false) => "pOH from [OH⁻]",
            (Self::HydroxideFromPoh, false) => "[OH⁻] from pOH",
            (Self::PhFromPoh, false) => "pH from pOH",
            (Self::PohFromPh, false) => "pOH from pH",
        }
    }

    pub fn prompt(self, danish: bool) -> &'static str {
        match (self, danish) {
            (Self::AmountFromMass, true) => "Masse i g; molarmasse i g/mol",
            (Self::MassFromAmount, true) => "Stofmængde i mol; molarmasse i g/mol",
            (Self::Concentration, true) => "Stofmængde i mol; volumen i L",
            (Self::Dilution, true) => {
                "Startkoncentration i mol/L; overført volumen i L; slutvolumen i L"
            }
            (Self::PhFromHydronium, true) => "[H₃O⁺] i mol/L",
            (Self::HydroniumFromPh, true) => "pH",
            (Self::PohFromHydroxide, true) => "[OH⁻] i mol/L",
            (Self::HydroxideFromPoh, true) => "pOH",
            (Self::PhFromPoh, true) => "pOH ved 25 °C",
            (Self::PohFromPh, true) => "pH ved 25 °C",
            (Self::AmountFromMass, false) => "Mass in g; molar mass in g/mol",
            (Self::MassFromAmount, false) => "Amount in mol; molar mass in g/mol",
            (Self::Concentration, false) => "Amount in mol; volume in L",
            (Self::Dilution, false) => {
                "Stock concentration in mol/L; transferred volume in L; final volume in L"
            }
            (Self::PhFromHydronium, false) => "[H₃O⁺] in mol/L",
            (Self::HydroniumFromPh, false) => "pH",
            (Self::PohFromHydroxide, false) => "[OH⁻] in mol/L",
            (Self::HydroxideFromPoh, false) => "pOH",
            (Self::PhFromPoh, false) => "pOH at 25 °C",
            (Self::PohFromPh, false) => "pH at 25 °C",
        }
    }

    pub fn example(self, danish: bool) -> &'static str {
        match (self, danish) {
            (Self::AmountFromMass, true) => "Eksempel: 5,00; 58,44",
            (Self::MassFromAmount, true) => "Eksempel: 0,100; 58,44",
            (Self::Concentration, true) => "Eksempel: 0,0250; 0,250",
            (Self::Dilution, true) => "Eksempel: 1,00; 0,0100; 0,100",
            (Self::PhFromHydronium, true) => "Eksempel: 0,0010",
            (Self::HydroniumFromPh, true) => "Eksempel: 3,00",
            (Self::PohFromHydroxide, true) => "Eksempel: 0,010",
            (Self::HydroxideFromPoh, true) => "Eksempel: 2,00",
            (Self::PhFromPoh, true) => "Eksempel: 4,00",
            (Self::PohFromPh, true) => "Eksempel: 10,00",
            (Self::AmountFromMass, false) => "Example: 5.00; 58.44",
            (Self::MassFromAmount, false) => "Example: 0.100; 58.44",
            (Self::Concentration, false) => "Example: 0.0250; 0.250",
            (Self::Dilution, false) => "Example: 1.00; 0.0100; 0.100",
            (Self::PhFromHydronium, false) => "Example: 0.0010",
            (Self::HydroniumFromPh, false) => "Example: 3.00",
            (Self::PohFromHydroxide, false) => "Example: 0.010",
            (Self::HydroxideFromPoh, false) => "Example: 2.00",
            (Self::PhFromPoh, false) => "Example: 4.00",
            (Self::PohFromPh, false) => "Example: 10.00",
        }
    }
}

fn parse_values(input: &str, count: usize) -> Result<Vec<f64>, String> {
    let values: Vec<f64> = input
        .split(';')
        .map(str::trim)
        .map(|field| {
            field
                .replace(',', ".")
                .parse::<f64>()
                .map_err(|_| format!("Invalid number: {field}"))
        })
        .collect::<Result<_, _>>()?;
    if values.len() != count {
        return Err(format!(
            "Enter exactly {count} value(s), separated by semicolons"
        ));
    }
    if values.iter().any(|v| !v.is_finite()) {
        return Err("Values must be finite numbers".into());
    }
    Ok(values)
}

fn require_positive(value: f64, label: &str) -> Result<(), String> {
    if value > 0.0 {
        Ok(())
    } else {
        Err(format!("{label} must be greater than zero"))
    }
}

fn require_nonnegative(value: f64, label: &str) -> Result<(), String> {
    if value >= 0.0 {
        Ok(())
    } else {
        Err(format!("{label} cannot be negative"))
    }
}

fn report(title: &str, equation: &str, substitution: &str, result: &str, note: &str) -> String {
    format!(
        "{title}\n\nEquation: {equation}\nSubstitution: {substitution}\nResult: {result}\n\n{note}"
    )
}

/// Format a calculated value to a consistent number of significant figures.
pub fn format_significant(value: f64, figures: u8) -> String {
    let figures = figures.clamp(2, 6) as usize;
    if value == 0.0 {
        return format!("{:.0}", value);
    }
    let exponent = value.abs().log10().floor() as i32;
    if !(-3..figures as i32).contains(&exponent) {
        format!("{:.*e}", figures - 1, value)
    } else {
        let decimals = (figures as i32 - 1 - exponent).max(0) as usize;
        format!("{value:.decimals$}")
    }
}

pub fn format_significant_locale(value: f64, figures: u8, danish: bool) -> String {
    let formatted = format_significant(value, figures);
    if danish {
        formatted.replace('.', ",")
    } else {
        formatted
    }
}

/// pH, pOH and pKa are logarithmic quantities; report decimal places rather
/// than applying ordinary significant-figure formatting to the log value.
pub fn format_log_quantity(value: f64, figures: u8, danish: bool) -> String {
    let decimals = figures.clamp(2, 6) as usize;
    let formatted = format!("{value:.decimals$}");
    if danish {
        formatted.replace('.', ",")
    } else {
        formatted
    }
}

#[cfg(test)]
fn calculate(tool: BasicCalculation, input: &str, danish: bool) -> Result<String, String> {
    calculate_with_significant_figures(tool, input, danish, 4)
}

pub fn calculate_with_significant_figures(
    tool: BasicCalculation,
    input: &str,
    danish: bool,
    figures: u8,
) -> Result<String, String> {
    let da = danish;
    let f = |value| format_significant_locale(value, figures, da);
    let log = |value| format_log_quantity(value, figures, da);
    let result = match tool {
        BasicCalculation::AmountFromMass => {
            let v = parse_values(input, 2)?;
            require_nonnegative(v[0], if da { "Masse" } else { "Mass" })?;
            require_positive(v[1], if da { "Molarmasse" } else { "Molar mass" })?;
            let n = v[0] / v[1];
            Ok(report(
                tool.title(da),
                "n = m / M",
                &format!("{} g / {} g mol⁻¹", f(v[0]), f(v[1])),
                &format!("{} mol", f(n)),
                if da {
                    "Kontrollér, at massen og molarmassen gælder samme stof."
                } else {
                    "Check that the mass and molar mass refer to the same substance."
                },
            ))
        }
        BasicCalculation::MassFromAmount => {
            let v = parse_values(input, 2)?;
            require_nonnegative(v[0], if da { "Stofmængde" } else { "Amount" })?;
            require_positive(v[1], if da { "Molarmasse" } else { "Molar mass" })?;
            let m = v[0] * v[1];
            Ok(report(
                tool.title(da),
                "m = n · M",
                &format!("{} mol · {} g mol⁻¹", f(v[0]), f(v[1])),
                &format!("{} g", f(m)),
                if da {
                    "Brug molarmassen for stoffets korrekte formel."
                } else {
                    "Use the molar mass for the correct chemical formula."
                },
            ))
        }
        BasicCalculation::Concentration => {
            let v = parse_values(input, 2)?;
            require_nonnegative(v[0], if da { "Stofmængde" } else { "Amount" })?;
            require_positive(v[1], if da { "Volumen" } else { "Volume" })?;
            let c = v[0] / v[1];
            Ok(report(
                tool.title(da),
                "c = n / V",
                &format!("{} mol / {} L", f(v[0]), f(v[1])),
                &format!("{} mol/L", f(c)),
                if da {
                    "Volumen skal angives i liter."
                } else {
                    "Enter volume in litres."
                },
            ))
        }
        BasicCalculation::Dilution => {
            let v = parse_values(input, 3)?;
            require_nonnegative(
                v[0],
                if da {
                    "Startkoncentration"
                } else {
                    "Stock concentration"
                },
            )?;
            require_nonnegative(
                v[1],
                if da {
                    "Overført volumen"
                } else {
                    "Transferred volume"
                },
            )?;
            require_positive(v[2], if da { "Slutvolumen" } else { "Final volume" })?;
            if v[2] < v[1] {
                return Err(if da {
                    "Slutvolumen kan ikke være mindre end det overførte volumen".into()
                } else {
                    "Final volume cannot be smaller than the transferred volume".into()
                });
            }
            let c2 = v[0] * v[1] / v[2];
            Ok(report(
                tool.title(da),
                "c₁V₁ = c₂V₂",
                &format!("c₂ = ({} mol/L · {} L) / {} L", f(v[0]), f(v[1]), f(v[2])),
                &format!("{} mol/L", f(c2)),
                if da {
                    "Antager, at stofmængden af opløst stof bevares ved fortynding."
                } else {
                    "Assumes the amount of dissolved substance is conserved during dilution."
                },
            ))
        }
        BasicCalculation::PhFromHydronium => {
            let v = parse_values(input, 1)?;
            require_positive(v[0], "[H₃O⁺]")?;
            let ph = -v[0].log10();
            Ok(report(
                tool.title(da),
                "pH = −log₁₀([H₃O⁺])",
                &format!("−log₁₀({} mol/L)", f(v[0])),
                &format!("pH = {}", log(ph)),
                if da {
                    "Bruger koncentration som tilnærmelse til aktivitet."
                } else {
                    "Uses concentration as an approximation to activity."
                },
            ))
        }
        BasicCalculation::HydroniumFromPh => {
            let v = parse_values(input, 1)?;
            if !(-5.0..=25.0).contains(&v[0]) {
                return Err("pH must be between -5 and 25".into());
            }
            let c = 10_f64.powf(-v[0]);
            Ok(report(
                tool.title(da),
                "[H₃O⁺] = 10⁻ᵖᴴ mol/L",
                &format!("10⁻{}", log(v[0])),
                &format!("{} mol/L", f(c)),
                if da {
                    "Koncentration beregnet fra pH; pH afhænger i praksis af aktivitet."
                } else {
                    "Concentration calculated from pH; pH is strictly related to activity."
                },
            ))
        }
        BasicCalculation::PohFromHydroxide => {
            let v = parse_values(input, 1)?;
            require_positive(v[0], "[OH⁻]")?;
            let poh = -v[0].log10();
            Ok(report(
                tool.title(da),
                "pOH = −log₁₀([OH⁻])",
                &format!("−log₁₀({} mol/L)", f(v[0])),
                &format!("pOH = {}", log(poh)),
                if da {
                    "Bruger koncentration som tilnærmelse til aktivitet."
                } else {
                    "Uses concentration as an approximation to activity."
                },
            ))
        }
        BasicCalculation::HydroxideFromPoh => {
            let v = parse_values(input, 1)?;
            if !(-5.0..=25.0).contains(&v[0]) {
                return Err("pOH must be between -5 and 25".into());
            }
            let c = 10_f64.powf(-v[0]);
            Ok(report(
                tool.title(da),
                "[OH⁻] = 10⁻ᵖᴼᴴ mol/L",
                &format!("10⁻{}", log(v[0])),
                &format!("{} mol/L", f(c)),
                if da {
                    "Koncentration beregnet fra pOH."
                } else {
                    "Concentration calculated from pOH."
                },
            ))
        }
        BasicCalculation::PhFromPoh | BasicCalculation::PohFromPh => {
            let v = parse_values(input, 1)?;
            let label = if tool == BasicCalculation::PhFromPoh {
                "pOH"
            } else {
                "pH"
            };
            if !(-5.0..=25.0).contains(&v[0]) {
                return Err(format!("{label} must be between -5 and 25"));
            }
            let output = 14.0 - v[0];
            let result_label = if tool == BasicCalculation::PhFromPoh {
                "pH"
            } else {
                "pOH"
            };
            Ok(report(
                tool.title(da),
                "pH + pOH = 14.00",
                &format!("14 − {}", log(v[0])),
                &format!("{result_label} = {}", log(output)),
                if da {
                    "Gælder for vand ved 25 °C med Kᵥ = 1,0×10⁻¹⁴."
                } else {
                    "For water at 25 °C with Kᵥ = 1.0×10⁻¹⁴."
                },
            ))
        }
    };
    result.map(|mut text| {
        text.push_str(&format!("\n\nDisplayed values use {figures} significant figures (pH/pOH use {figures} decimal places); calculations use unrounded inputs."));
        text
    })
}

/// Evaluate limiting reagent and theoretical/actual yield for a balanced
/// equation. Reactant coefficients and amounts are aligned in input order.
pub fn stoichiometric_yield(input: &str, danish: bool, figures: u8) -> Result<String, String> {
    let f = |value| format_significant_locale(value, figures, danish);
    let fields: Vec<&str> = input.split(';').map(str::trim).collect();
    if !(4..=5).contains(&fields.len()) {
        return Err(if danish {
            "Indtast koefficienter; stofmængder i mol; produktkoefficient; produktets molarmasse i g/mol; eventuelt faktisk produktmasse i g".into()
        } else {
            "Enter coefficients; reactant amounts in mol; product coefficient; product molar mass in g/mol; optional actual product mass in g".into()
        });
    }
    let parse_list = |text: &str| -> Result<Vec<f64>, String> {
        text.split_whitespace()
            .map(|token| {
                token
                    .replace(',', ".")
                    .parse::<f64>()
                    .map_err(|_| format!("Invalid number: {token}"))
            })
            .collect()
    };
    let coefficients = parse_list(fields[0])?;
    let amounts = parse_list(fields[1])?;
    if coefficients.len() < 2 || coefficients.len() != amounts.len() {
        return Err(if danish {
            "Angiv mindst to koefficienter og lige så mange stofmængder.".into()
        } else {
            "Enter at least two coefficients and the same number of reactant amounts.".into()
        });
    }
    if coefficients
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0 || value.fract() != 0.0)
        || amounts
            .iter()
            .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(if danish {
            "Koefficienter skal være positive heltal, og stofmængder positive endelige tal.".into()
        } else {
            "Coefficients must be positive integers and amounts must be finite positive numbers."
                .into()
        });
    }
    let parse_number = |field: &str| {
        field
            .replace(',', ".")
            .parse::<f64>()
            .map_err(|_| format!("Invalid number: {field}"))
    };
    let product_coefficient = parse_number(fields[2])?;
    let molar_mass = parse_number(fields[3])?;
    if !product_coefficient.is_finite()
        || product_coefficient <= 0.0
        || product_coefficient.fract() != 0.0
    {
        return Err("Product coefficient must be a positive integer".into());
    }
    if !molar_mass.is_finite() {
        return Err("Molar mass must be finite".into());
    }
    require_positive(molar_mass, if danish { "Molarmasse" } else { "Molar mass" })?;
    let limiting_index = amounts
        .iter()
        .zip(&coefficients)
        .enumerate()
        .min_by(|(_, (a, ca)), (_, (b, cb))| (*a / *ca).total_cmp(&(*b / *cb)))
        .map(|(index, _)| index)
        .ok_or("No reactants were entered")?;
    let extent = amounts[limiting_index] / coefficients[limiting_index];
    let product_moles = extent * product_coefficient;
    let theoretical_mass = product_moles * molar_mass;
    let mut output = format!(
        "{}\n\n{}\n{}\n{}\n{}\n{}\n{}",
        if danish {
            "Støkiometrisk udbytte"
        } else {
            "Stoichiometric yield"
        },
        if danish {
            "Begrænsende reaktant"
        } else {
            "Limiting reactant"
        },
        if danish {
            format!(
                "Reaktant {} (n/ν = {} mol / {} = {} mol reaktionsomfang)",
                limiting_index + 1,
                f(amounts[limiting_index]),
                f(coefficients[limiting_index]),
                f(extent)
            )
        } else {
            format!(
                "Reactant {} (n/ν = {} mol / {} = {} mol reaction extent)",
                limiting_index + 1,
                f(amounts[limiting_index]),
                f(coefficients[limiting_index]),
                f(extent)
            )
        },
        if danish {
            format!("Teoretisk produktmængde = {} mol", f(product_moles))
        } else {
            format!("Theoretical product amount = {} mol", f(product_moles))
        },
        if danish {
            format!("Teoretisk produktmasse = {} g", f(theoretical_mass))
        } else {
            format!("Theoretical product mass = {} g", f(theoretical_mass))
        },
        if fields.len() == 5 {
            let actual = parse_number(fields[4])?;
            if !actual.is_finite() {
                return Err("Actual mass must be finite".into());
            }
            require_nonnegative(
                actual,
                if danish {
                    "Faktisk masse"
                } else {
                    "Actual mass"
                },
            )?;
            let percent = actual / theoretical_mass * 100.0;
            if danish {
                format!(
                    "Faktisk masse = {} g\nProcentudbytte = {} %",
                    f(actual),
                    f(percent)
                )
            } else {
                format!(
                    "Actual mass = {} g\nPercent yield = {} %",
                    f(actual),
                    f(percent)
                )
            }
        } else {
            if danish {
                "Faktisk udbytte: ikke angivet".into()
            } else {
                "Actual yield: not entered".into()
            }
        },
        if danish {
            "Forudsætning: den indtastede ligning er afstemt, og koefficienterne følger reaktanterne i samme rækkefølge. Renheds- og sideproduktkorrektioner er ikke medregnet."
        } else {
            "Assumption: the entered equation is balanced, and coefficients match reactants in the same order. Purity and side-product corrections are not included."
        }
    );
    output.push_str(&if danish {
        format!("\n\nViste værdier er afrundet til {figures} betydende cifre; beregninger bruger uafrundede input.")
    } else {
        format!("\n\nDisplayed values rounded to {figures} significant figures; calculations use unrounded inputs.")
    });
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stoichiometry_calculators_show_equation_substitution_and_units() {
        let r = calculate(BasicCalculation::AmountFromMass, "5,844;58,44", true).unwrap();
        assert!(r.contains("n = m / M") && r.contains("0,1000 mol"));
        let r = calculate(BasicCalculation::MassFromAmount, "0.1;58.44", false).unwrap();
        assert!(r.contains("5.844 g"));
        let r = calculate(BasicCalculation::Concentration, "0.025;0.25", true).unwrap();
        assert!(r.contains("0,1000 mol/L"));
    }

    #[test]
    fn dilution_and_acid_base_conversions_match_reference_values() {
        let r = calculate(BasicCalculation::Dilution, "1;0.01;0.1", false).unwrap();
        assert!(r.contains("0.1000 mol/L"));
        let r = calculate(BasicCalculation::PhFromHydronium, "0.001", false).unwrap();
        assert!(r.contains("pH = 3.0000"));
        let r = calculate(BasicCalculation::HydroniumFromPh, "3", false).unwrap();
        assert!(r.contains("0.001000 mol/L"));
        let r = calculate(BasicCalculation::PohFromHydroxide, "0.01", false).unwrap();
        assert!(r.contains("pOH = 2.0000"));
        let r = calculate(BasicCalculation::HydroxideFromPoh, "2", false).unwrap();
        assert!(r.contains("0.01000 mol/L"));
        assert!(calculate(BasicCalculation::PhFromPoh, "4", false)
            .unwrap()
            .contains("pH = 10.0000"));
        assert!(calculate(BasicCalculation::PohFromPh, "10", false)
            .unwrap()
            .contains("pOH = 4.0000"));
    }

    #[test]
    fn invalid_units_ranges_and_volume_relations_are_rejected() {
        assert!(calculate(BasicCalculation::Concentration, "1;0", true).is_err());
        assert!(calculate(BasicCalculation::HydroniumFromPh, "NaN", false).is_err());
        assert!(calculate(BasicCalculation::Dilution, "1;0.2;0.1", false).is_err());
        assert!(calculate(BasicCalculation::AmountFromMass, "-1;58", false).is_err());
        assert!(calculate(BasicCalculation::PhFromPoh, "30", false).is_err());
        assert!(calculate(BasicCalculation::Concentration, "1;2;3", false).is_err());
    }

    #[test]
    fn significant_figure_formatter_handles_small_and_large_values() {
        assert_eq!(format_significant(0.012345, 3), "0.0123");
        assert_eq!(format_significant(12_345.0, 4), "1.234e4");
        assert_eq!(format_significant(0.000_0123, 3), "1.23e-5");
        assert_eq!(format_significant(0.1, 4), "0.1000");
        assert_eq!(format_log_quantity(3.0, 4, false), "3.0000");
        assert_eq!(format_log_quantity(5.061, 4, true), "5,0610");
    }

    #[test]
    fn limiting_reactant_and_product_yield_follow_balanced_coefficients() {
        let report = stoichiometric_yield("2 1; 0.30 0.10; 2; 18.015; 3.20", false, 4).unwrap();
        assert!(report.contains("Reactant 2"));
        assert!(report.contains("0.2000 mol"));
        assert!(report.contains("3.603 g"));
        assert!(report.contains("Percent yield = "), "{report}");
        assert!(stoichiometric_yield("1 2; 0.3; 1; 10", false, 4).is_err());
        assert!(stoichiometric_yield("1 2; 0.3 0.1; 1; inf", false, 4).is_err());
    }
}
