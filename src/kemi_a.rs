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

pub fn calculate(tool: BasicCalculation, input: &str, danish: bool) -> Result<String, String> {
    let da = danish;
    match tool {
        BasicCalculation::AmountFromMass => {
            let v = parse_values(input, 2)?;
            require_nonnegative(v[0], if da { "Masse" } else { "Mass" })?;
            require_positive(v[1], if da { "Molarmasse" } else { "Molar mass" })?;
            let n = v[0] / v[1];
            Ok(report(
                tool.title(da),
                "n = m / M",
                &format!("{:.6} g / {:.6} g mol⁻¹", v[0], v[1]),
                &format!("{n:.6} mol"),
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
                &format!("{:.6} mol · {:.6} g mol⁻¹", v[0], v[1]),
                &format!("{m:.6} g"),
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
                &format!("{:.6} mol / {:.6} L", v[0], v[1]),
                &format!("{c:.6} mol/L"),
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
                &format!("c₂ = ({:.6} mol/L · {:.6} L) / {:.6} L", v[0], v[1], v[2]),
                &format!("{c2:.6} mol/L"),
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
                &format!("−log₁₀({:.6} mol/L)", v[0]),
                &format!("pH = {ph:.4}"),
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
                &format!("10⁻{:.4}", v[0]),
                &format!("{c:.6e} mol/L"),
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
                &format!("−log₁₀({:.6} mol/L)", v[0]),
                &format!("pOH = {poh:.4}"),
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
                &format!("10⁻{:.4}", v[0]),
                &format!("{c:.6e} mol/L"),
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
                &format!("14.00 − {:.4}", v[0]),
                &format!("{result_label} = {output:.4}"),
                if da {
                    "Gælder for vand ved 25 °C med Kᵥ = 1,0×10⁻¹⁴."
                } else {
                    "For water at 25 °C with Kᵥ = 1.0×10⁻¹⁴."
                },
            ))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stoichiometry_calculators_show_equation_substitution_and_units() {
        let r = calculate(BasicCalculation::AmountFromMass, "5,844;58,44", true).unwrap();
        assert!(r.contains("n = m / M") && r.contains("0.100000 mol"));
        let r = calculate(BasicCalculation::MassFromAmount, "0.1;58.44", false).unwrap();
        assert!(r.contains("5.844000 g"));
        let r = calculate(BasicCalculation::Concentration, "0.025;0.25", true).unwrap();
        assert!(r.contains("0.100000 mol/L"));
    }

    #[test]
    fn dilution_and_acid_base_conversions_match_reference_values() {
        let r = calculate(BasicCalculation::Dilution, "1;0.01;0.1", false).unwrap();
        assert!(r.contains("0.100000 mol/L"));
        let r = calculate(BasicCalculation::PhFromHydronium, "0.001", false).unwrap();
        assert!(r.contains("pH = 3.0000"));
        let r = calculate(BasicCalculation::HydroniumFromPh, "3", false).unwrap();
        assert!(r.contains("1.000000e-3 mol/L"));
        let r = calculate(BasicCalculation::PohFromHydroxide, "0.01", false).unwrap();
        assert!(r.contains("pOH = 2.0000"));
        let r = calculate(BasicCalculation::HydroxideFromPoh, "2", false).unwrap();
        assert!(r.contains("1.000000e-2 mol/L"));
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
}
