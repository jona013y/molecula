//! Aqueous acid/base calculations for the Danish Chemistry A curriculum.
//!
//! Equilibria use ideal concentrations and Kw = 1e-14 at 25 °C. pKa values
//! are supplied by the learner from a cited source or exercise; this module
//! never predicts or invents a pKa.

const KW: f64 = 1e-14;
const MIN_PH: f64 = -2.0;
const MAX_PH: f64 = 16.0;

fn positive_finite(value: f64, label: &str) -> Result<(), String> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(format!("{label} must be a finite positive number"))
    }
}

fn nonnegative_finite(value: f64, label: &str) -> Result<(), String> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(format!("{label} must be finite and non-negative"))
    }
}

fn validate_pka(pka: f64) -> Result<(), String> {
    if pka.is_finite() && (-20.0..=30.0).contains(&pka) {
        Ok(())
    } else {
        Err("pKa must be finite and between -20 and 30".into())
    }
}

fn validate_pkas(pkas: &[f64]) -> Result<(), String> {
    if pkas.is_empty() || pkas.len() > 8 {
        return Err("Enter between one and eight successive pKa values".into());
    }
    for (index, &pka) in pkas.iter().enumerate() {
        validate_pka(pka)?;
        if index > 0 && pka <= pkas[index - 1] {
            return Err("Polyprotic pKa values must be in strictly increasing order".into());
        }
    }
    Ok(())
}

/// Calculate pH from a strong acid and strong base after mixing, given their
/// final analytical concentrations in the mixed solution.
pub fn strong_mixture_ph(acid_molar: f64, base_molar: f64) -> Result<f64, String> {
    nonnegative_finite(acid_molar, "Strong acid concentration")?;
    nonnegative_finite(base_molar, "Strong base concentration")?;
    let excess = acid_molar - base_molar;
    let discriminant_half = (excess * 0.5).hypot(KW.sqrt());
    let h = if excess >= 0.0 {
        excess * 0.5 + discriminant_half
    } else {
        KW / (-excess * 0.5 + discriminant_half)
    };
    Ok(-h.log10())
}

/// Henderson–Hasselbalch pH for a conjugate acid/base buffer.
pub fn buffer_ph(pka: f64, acid_molar: f64, base_molar: f64) -> Result<f64, String> {
    validate_pka(pka)?;
    positive_finite(acid_molar, "Conjugate acid concentration")?;
    positive_finite(base_molar, "Conjugate base concentration")?;
    Ok(pka + base_molar.log10() - acid_molar.log10())
}

/// Derive pKa from a measured buffer pH and the conjugate base/acid ratio.
pub fn pka_from_buffer(ph: f64, acid_molar: f64, base_molar: f64) -> Result<f64, String> {
    if !ph.is_finite() || !(-5.0..=25.0).contains(&ph) {
        return Err("Measured pH must be finite and between -5 and 25".into());
    }
    positive_finite(acid_molar, "Conjugate acid concentration")?;
    positive_finite(base_molar, "Conjugate base concentration")?;
    Ok(ph - base_molar.log10() + acid_molar.log10())
}

fn fractions_at_ph(ph: f64, pkas: &[f64]) -> Vec<f64> {
    let log10 = std::f64::consts::LN_10;
    let mut weights = Vec::with_capacity(pkas.len() + 1);
    let mut cumulative_pka = 0.0;
    for degree in 0..=pkas.len() {
        if degree > 0 {
            cumulative_pka += pkas[degree - 1];
        }
        weights.push(
            (f64::from(degree as u32) * ph * log10 - cumulative_pka * log10).clamp(-745.0, 700.0),
        );
    }
    let maximum = weights.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let total: f64 = weights.iter().map(|weight| (weight - maximum).exp()).sum();
    weights
        .into_iter()
        .map(|weight| (weight - maximum).exp() / total)
        .collect()
}

/// Fractions of H_nA, H_(n-1)A⁻, …, Aⁿ⁻ at a requested pH.
pub fn polyprotic_fractions(ph: f64, pkas: &[f64]) -> Result<Vec<f64>, String> {
    validate_pkas(pkas)?;
    if !ph.is_finite() || !(-5.0..=25.0).contains(&ph) {
        return Err("pH must be finite and between -5 and 25".into());
    }
    Ok(fractions_at_ph(ph, pkas))
}

fn equilibrium_ph(mut residual: impl FnMut(f64) -> f64) -> Result<f64, String> {
    let mut low = MIN_PH;
    let mut high = MAX_PH;
    let low_value = residual(low);
    let high_value = residual(high);
    if !low_value.is_finite() || !high_value.is_finite() || low_value < 0.0 || high_value > 0.0 {
        return Err(
            "Equilibrium could not be bracketed in the supported pH range (-2 to 16)".into(),
        );
    }
    for _ in 0..160 {
        let mid = (low + high) / 2.0;
        if residual(mid) > 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    Ok((low + high) / 2.0)
}

/// Exact ideal-solution pH for a weak monoprotic acid, including water
/// autoionization. pKa is the learner-supplied aqueous value at 25 °C.
pub fn weak_acid_ph(pka: f64, total_molar: f64) -> Result<f64, String> {
    validate_pka(pka)?;
    positive_finite(total_molar, "Weak acid concentration")?;
    let ka = 10_f64.powf(-pka);
    equilibrium_ph(|ph| {
        let h = 10_f64.powf(-ph);
        h - KW / h - total_molar * ka / (ka + h)
    })
}

/// Exact ideal-solution pH for a weak base, entered with the pKa of its
/// conjugate acid (BH⁺). Includes water autoionization.
pub fn weak_base_ph(conjugate_acid_pka: f64, total_molar: f64) -> Result<f64, String> {
    validate_pka(conjugate_acid_pka)?;
    positive_finite(total_molar, "Weak base concentration")?;
    let ka = 10_f64.powf(-conjugate_acid_pka);
    equilibrium_ph(|ph| {
        let h = 10_f64.powf(-ph);
        h - KW / h + total_molar * h / (h + ka)
    })
}

fn validate_titration_inputs(
    pkas: &[f64],
    acid_volume_l: f64,
    acid_molar: f64,
    base_molar: f64,
    base_volume_l: f64,
) -> Result<(), String> {
    validate_pkas(pkas)?;
    positive_finite(acid_volume_l, "Initial acid volume")?;
    positive_finite(acid_molar, "Initial acid concentration")?;
    positive_finite(base_molar, "Strong base concentration")?;
    nonnegative_finite(base_volume_l, "Added base volume")?;
    Ok(())
}

/// pH after titrating a polyprotic acid with a strong base. Volumes are in L
/// and concentrations in mol/L. Charge balance is solved numerically.
pub fn polyprotic_titration_ph(
    pkas: &[f64],
    acid_volume_l: f64,
    acid_molar: f64,
    base_molar: f64,
    base_volume_l: f64,
) -> Result<f64, String> {
    validate_titration_inputs(pkas, acid_volume_l, acid_molar, base_molar, base_volume_l)?;
    let total_volume = acid_volume_l + base_volume_l;
    if !total_volume.is_finite() || total_volume <= 0.0 {
        return Err("Total titration volume is invalid".into());
    }
    let acid_total = acid_volume_l * acid_molar / total_volume;
    let sodium = base_volume_l * base_molar / total_volume;
    equilibrium_ph(|ph| {
        let h = 10_f64.powf(-ph);
        let hydroxide = KW / h;
        let mean_charge: f64 = fractions_at_ph(ph, pkas)
            .iter()
            .enumerate()
            .map(|(charge, fraction)| charge as f64 * fraction)
            .sum();
        h + sodium - hydroxide - acid_total * mean_charge
    })
}

/// A table suitable for a Bjerrum diagram or titration-curve plot.
pub fn bjerrum_points(
    pkas: &[f64],
    ph_min: f64,
    ph_max: f64,
    step: f64,
) -> Result<Vec<(f64, Vec<f64>)>, String> {
    validate_pkas(pkas)?;
    if !ph_min.is_finite()
        || !ph_max.is_finite()
        || !step.is_finite()
        || ph_min < -5.0
        || ph_max > 25.0
        || ph_min >= ph_max
        || step <= 0.0
    {
        return Err("Use a finite range within pH -5 to 25 and a positive step".into());
    }
    let count = ((ph_max - ph_min) / step).ceil() as usize + 1;
    if count > 501 {
        return Err("Bjerrum table is limited to 501 pH points; increase the step".into());
    }
    Ok((0..count)
        .map(|index| {
            let ph = (ph_min + index as f64 * step).min(ph_max);
            (ph, fractions_at_ph(ph, pkas))
        })
        .collect())
}

/// Generate a pH curve for titration of H_nA by a strong base.
pub fn polyprotic_titration_curve(
    pkas: &[f64],
    acid_volume_l: f64,
    acid_molar: f64,
    base_molar: f64,
    max_equivalents: f64,
    points: usize,
) -> Result<Vec<(f64, f64)>, String> {
    validate_pkas(pkas)?;
    positive_finite(acid_volume_l, "Initial acid volume")?;
    positive_finite(acid_molar, "Initial acid concentration")?;
    positive_finite(base_molar, "Strong base concentration")?;
    if !max_equivalents.is_finite()
        || max_equivalents <= 0.0
        || max_equivalents > pkas.len() as f64 + 2.0
    {
        return Err("Maximum titration extent must be positive and no more than two equivalents beyond complete deprotonation".into());
    }
    if !(2..=501).contains(&points) {
        return Err("A titration curve must contain 2 to 501 points".into());
    }
    let acid_moles = acid_volume_l * acid_molar;
    (0..points)
        .map(|index| {
            let equivalents = max_equivalents * index as f64 / (points - 1) as f64;
            let base_volume = equivalents * acid_moles / base_molar;
            polyprotic_titration_ph(pkas, acid_volume_l, acid_molar, base_molar, base_volume)
                .map(|ph| (base_volume, ph))
        })
        .collect()
}

/// Analyze measured monoprotic acid titration data. Inputs are added titrant
/// volume (mL) and measured pH. Equivalence is estimated by the largest
/// adjacent-point slope; pKa is interpolated at half that volume.
pub fn analyze_monoprotic_titration_data(data: &[(f64, f64)]) -> Result<(f64, f64, f64), String> {
    if data.len() < 5 || data.len() > 2000 {
        return Err("Enter 5 to 2000 measured volume,pH pairs".into());
    }
    for (i, &(volume, ph)) in data.iter().enumerate() {
        if !volume.is_finite() || volume < 0.0 || !ph.is_finite() || !(-2.0..=16.0).contains(&ph) {
            return Err(
                "Volumes must be finite and non-negative; measured pH must be between -2 and 16"
                    .into(),
            );
        }
        if i > 0 && volume <= data[i - 1].0 {
            return Err("Titrant volumes must be strictly increasing".into());
        }
    }
    let (index, slope) = data
        .windows(2)
        .enumerate()
        .map(|(i, pair)| (i, (pair[1].1 - pair[0].1) / (pair[1].0 - pair[0].0)))
        .max_by(|a, b| a.1.total_cmp(&b.1))
        .ok_or("No adjacent data points")?;
    if slope <= 0.0 {
        return Err("Measured pH must increase overall for a strong-base titration".into());
    }
    let (v0, _) = data[index];
    let (v1, _) = data[index + 1];
    let eq = (v0 + v1) / 2.0;
    let half = eq / 2.0;
    let bracket = data
        .windows(2)
        .find(|pair| pair[0].0 <= half && half <= pair[1].0)
        .ok_or("Measurements must bracket the half-equivalence volume to estimate pKa")?;
    let fraction = (half - bracket[0].0) / (bracket[1].0 - bracket[0].0);
    let pka = bracket[0].1 + fraction * (bracket[1].1 - bracket[0].1);
    Ok((eq, pka, slope))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strong_mixture_uses_water_balance_and_neutralization() {
        assert!((strong_mixture_ph(0.1, 0.0).unwrap() - 1.0).abs() < 1e-10);
        assert!((strong_mixture_ph(0.0, 0.1).unwrap() - 13.0).abs() < 1e-10);
        assert!((strong_mixture_ph(0.1, 0.1).unwrap() - 7.0).abs() < 1e-10);
        assert!((strong_mixture_ph(0.0, 0.0).unwrap() - 7.0).abs() < 1e-10);
    }

    #[test]
    fn buffer_and_experimental_pka_are_inverse_henderson_hasselbalch_calculations() {
        let ph = buffer_ph(4.76, 0.1, 0.2).unwrap();
        assert!((ph - 5.061_029_995_7).abs() < 1e-9);
        assert!((pka_from_buffer(ph, 0.1, 0.2).unwrap() - 4.76).abs() < 1e-12);
        assert_eq!(buffer_ph(4.76, 0.1, 0.1).unwrap(), 4.76);
    }

    #[test]
    fn weak_acid_and_base_ph_solve_equilibrium_with_water() {
        let acid_ph = weak_acid_ph(4.76, 0.1).unwrap();
        assert!((acid_ph - 2.88).abs() < 0.02);
        let base_ph = weak_base_ph(9.25, 0.1).unwrap();
        assert!((base_ph - 11.13).abs() < 0.03);
    }

    #[test]
    fn polyprotic_fractions_sum_to_one_and_are_equal_at_the_first_pka() {
        let fractions = polyprotic_fractions(2.15, &[2.15, 7.20, 12.35]).unwrap();
        assert!((fractions.iter().sum::<f64>() - 1.0).abs() < 1e-12);
        assert!((fractions[0] - fractions[1]).abs() < 1e-12);
    }

    #[test]
    fn buffer_half_equivalence_and_equivalence_are_solved_by_charge_balance() {
        let half = polyprotic_titration_ph(&[4.76], 0.025, 0.1, 0.1, 0.0125).unwrap();
        let equivalent = polyprotic_titration_ph(&[4.76], 0.025, 0.1, 0.1, 0.025).unwrap();
        assert!((half - 4.76).abs() < 0.03);
        assert!(equivalent > 7.0);
    }

    #[test]
    fn bjerrum_and_titration_curves_have_bounded_finite_points() {
        let bjerrum = bjerrum_points(&[4.76, 10.33], 0.0, 14.0, 1.0).unwrap();
        assert_eq!(bjerrum.len(), 15);
        assert!(bjerrum
            .iter()
            .all(|(_, fractions)| { (fractions.iter().sum::<f64>() - 1.0).abs() < 1e-12 }));
        let curve = polyprotic_titration_curve(&[4.76], 0.025, 0.1, 0.1, 1.5, 31).unwrap();
        assert_eq!(curve.len(), 31);
        assert!(curve
            .iter()
            .all(|(volume, ph)| volume.is_finite() && ph.is_finite()));
    }

    #[test]
    fn measured_titration_analysis_finds_equivalence_and_half_equivalence_pka() {
        let data = [
            (0.0, 2.1),
            (2.0, 2.4),
            (4.0, 2.8),
            (6.0, 3.2),
            (8.0, 4.1),
            (9.0, 5.0),
            (10.0, 8.8),
            (11.0, 10.8),
            (12.0, 11.5),
        ];
        let (equivalence, pka, slope) = analyze_monoprotic_titration_data(&data).unwrap();
        assert_eq!(equivalence, 9.5);
        assert!((pka - 2.95).abs() < 1e-12);
        assert!((slope - 3.8).abs() < 1e-12);
        assert!(analyze_monoprotic_titration_data(&data[..4]).is_err());
        assert!(analyze_monoprotic_titration_data(&[
            (0.0, 2.0),
            (1.0, 3.0),
            (1.0, 4.0),
            (2.0, 5.0),
            (3.0, 6.0)
        ])
        .is_err());
    }

    #[test]
    fn invalid_concentrations_and_unsorted_pkas_are_rejected() {
        assert!(weak_acid_ph(4.76, 0.0).is_err());
        assert!(buffer_ph(4.76, 0.1, 0.0).is_err());
        assert!(polyprotic_fractions(7.0, &[7.0, 4.0]).is_err());
        assert!(strong_mixture_ph(f64::NAN, 0.1).is_err());
    }
}
