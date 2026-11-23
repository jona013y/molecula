//! Curated, offline study content and bounded intent matching for the local tutor.
//! This module deliberately does not generate chemistry claims from free text.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Topic {
    GettingStarted,
    Buffer,
    AcidBase,
    Titration,
    AmountConcentration,
    Reaction,
    Naming,
    Editor,
}

impl Topic {
    pub const ALL: [Self; 8] = [
        Self::GettingStarted,
        Self::Buffer,
        Self::AcidBase,
        Self::Titration,
        Self::AmountConcentration,
        Self::Reaction,
        Self::Naming,
        Self::Editor,
    ];

    pub fn title(self, da: bool) -> &'static str {
        match (self, da) {
            (Self::GettingStarted, true) => "Kom godt i gang",
            (Self::Buffer, true) => "Buffere og pKa",
            (Self::AcidBase, true) => "Syrer, baser og pH",
            (Self::Titration, true) => "Titrering",
            (Self::AmountConcentration, true) => "Stofmængde og koncentration",
            (Self::Reaction, true) => "Reaktioner og afstemning",
            (Self::Naming, true) => "Navngivning",
            (Self::Editor, true) => "Molecula-editoren",
            (Self::GettingStarted, false) => "Getting started",
            (Self::Buffer, false) => "Buffers and pKa",
            (Self::AcidBase, false) => "Acids, bases and pH",
            (Self::Titration, false) => "Titration",
            (Self::AmountConcentration, false) => "Amount and concentration",
            (Self::Reaction, false) => "Reactions and balancing",
            (Self::Naming, false) => "Chemical naming",
            (Self::Editor, false) => "The Molecula editor",
        }
    }

    pub fn article(self, da: bool) -> &'static str {
        match (self, da) {
            (Self::GettingStarted, true) => "Vælg et emne til venstre, eller spørg den lokale kemiassistent.\n\nBeregningerne viser formel, indsatte værdier, enhed og antagelser. Du kan bruge decimal-komma eller punktum.\n\nHjælp og assistent kører lokalt. Der sendes ikke spørgsmål eller tegninger til en tjeneste.",
            (Self::GettingStarted, false) => "Choose a topic on the left, or ask the local chemistry helper.\n\nCalculations show the equation, substituted values, units and assumptions. Decimal commas and points are accepted.\n\nThe help and tutor run locally. Questions and drawings are not sent to a service.",
            (Self::Buffer, true) => "En buffer indeholder en svag syre HA og dens korresponderende base A⁻.\n\npH = pKa + log₁₀([A⁻]/[HA])\n\nEksempel: pKa 4,76, [HA] = 0,10 mol/L og [A⁻] = 0,20 mol/L giver pH = 4,76 + log₁₀(2) ≈ 5,06.\n\nNår [A⁻] = [HA], er pH = pKa. Brug pKa fra en måling eller en kilde; Molecula udleder ikke pKa fra tegningen.",
            (Self::Buffer, false) => "A buffer contains a weak acid HA and its conjugate base A⁻.\n\npH = pKa + log₁₀([A⁻]/[HA])\n\nExample: pKa 4.76, [HA] = 0.10 mol/L and [A⁻] = 0.20 mol/L gives pH = 4.76 + log₁₀(2) ≈ 5.06.\n\nWhen [A⁻] = [HA], pH = pKa. Use pKa from a measurement or source; Molecula does not infer pKa from a drawing.",
            (Self::AcidBase, true) => "pH = −log₁₀([H₃O⁺]) og pOH = −log₁₀([OH⁻]).\n\nI vand ved 25 °C gælder pH + pOH = 14,00, fordi Kᵥ = 1,0×10⁻¹⁴.\n\nEn stærk syre antages fuldstændigt protolyseret i disse skoleberegninger. En svag monoprot syre beregnes ud fra ligevægten HA ⇌ H₃O⁺ + A⁻. Ved fortyndede, ideelle opløsninger kan aktivitet tilnærmes med koncentration.\n\nKemi A-værktøjet antager ideel vandig opløsning ved 25 °C; det modellerer ikke ionstyrke eller aktivitetskoefficienter.",
            (Self::AcidBase, false) => "pH = −log₁₀([H₃O⁺]) and pOH = −log₁₀([OH⁻]).\n\nFor water at 25 °C, pH + pOH = 14.00 because Kᵥ = 1.0×10⁻¹⁴.\n\nA strong acid is assumed to dissociate completely in these classroom calculations. A weak monoprotic acid is calculated from HA ⇌ H₃O⁺ + A⁻. In dilute ideal solutions, activity is approximated by concentration.\n\nThe Chemistry A tools assume an ideal aqueous solution at 25 °C; they do not model ionic strength or activity coefficients.",
            (Self::Titration, true) => "Ved en syre-base-titrering tilsætter man en opløsning med kendt koncentration til analytten.\n\nVed ækvivalenspunktet er stofmængderne af syre og base reageret i det støkiometriske forhold fra reaktionsskemaet. For en monoprot syre og en stærk base er n(OH⁻) = n(HA) ved ækvivalens.\n\nFør ækvivalenspunktet kan en buffer ofte beskrive opløsningen; ved ækvivalens afhænger pH af produkternes syre-base-egenskaber. Efter ækvivalens bestemmer overskuddet af titrator pH.\n\nMoleculas titreringskurve for polyprot syre bruger ideal ligevægt ved 25 °C og kræver pKa-værdier.",
            (Self::Titration, false) => "In an acid-base titration, a solution of known concentration is added to the analyte.\n\nAt the equivalence point, the amounts of acid and base have reacted in the stoichiometric ratio from the reaction equation. For a monoprotic acid and strong base, n(OH⁻) = n(HA) at equivalence.\n\nBefore equivalence, the solution may behave as a buffer; at equivalence, product acid-base properties determine pH. After equivalence, excess titrant determines pH.\n\nMolecula's polyprotic titration curve uses ideal equilibrium at 25 °C and requires pKa values.",
            (Self::AmountConcentration, true) => "Stofmængde og koncentration forbindes med n = m/M og c = n/V.\n\nn måles i mol, m i g, M i g/mol, V i L og c i mol/L.\n\nVed fortynding bevares stofmængden af opløst stof: c₁V₁ = c₂V₂. Brug samme volumen-enhed på begge sider.\n\nVælg Beregninger → Kemi A: stofmængde, fortynding eller koncentration for at få trin og enheder vist.",
            (Self::AmountConcentration, false) => "Amount and concentration are related by n = m/M and c = n/V.\n\nn is in mol, m in g, M in g/mol, V in L and c in mol/L.\n\nDuring dilution, the amount of dissolved substance is conserved: c₁V₁ = c₂V₂. Use the same volume unit on both sides.\n\nChoose Calculations → Chemistry A: amount, dilution or concentration to see the worked steps and units.",
            (Self::Reaction, true) => "En reaktions ligning afstemmes ved at bevare hvert grundstofs antal atomer og den samlede ladning. Koefficienter ændres; kemiske indeks i formler ændres ikke.\n\nI editoren: tegn reaktanter og produkter som adskilte komponenter, sæt deres rolle i inspektøren, og tilføj reaktionspil. Vælg Beregninger → Afstem reaktion eller Reaktionsanalyse for at kontrollere grafens formler, molarmasser og atom-/ladningsbalance.\n\nBrug derefter Begrænsende reaktant og udbytte til at indtaste de afstemte reaktantkoefficienter og de kendte stofmængder. Mængder kan ikke udledes af tegningen, så de indtastes særskilt.\n\nHvis der findes flere løsninger eller utilstrækkelige oplysninger, afviser Molecula automatisk ændringen i stedet for at gætte.",
            (Self::Reaction, false) => "A reaction equation is balanced by conserving the number of atoms of every element and total charge. Change coefficients; do not change subscripts in formulas.\n\nIn the editor, draw reactants and products as separate components, assign their role in the inspector, and add a reaction arrow. Choose Calculations → Auto-balance reaction or Reaction analysis to check graph formulas, molar masses, and atom/charge balance.\n\nThen use Limiting reactant & yield to enter the balanced reactant coefficients and known amounts. Quantities cannot be inferred from a drawing, so enter them separately.\n\nIf multiple solutions or insufficient information exist, Molecula refuses to guess.",
            (Self::Naming, true) => "Moleculas navngivningslektion viser de trin, som den indbyggede navngiver faktisk understøtter, og fremhæver atomerne.\n\nNavngiveren er et undervisningsudsnit, ikke en komplet IUPAC-motor. Hvis en struktur ikke understøttes, viser Molecula årsagen i stedet for at præsentere et usikkert navn som korrekt.\n\nTegn en støttet struktur, og vælg Hjælp → Navngivningslektion.",
            (Self::Naming, false) => "Molecula's naming lesson shows the steps supported by its built-in namer and highlights the atoms.\n\nThe namer is a teaching subset, not a complete IUPAC engine. If a structure is unsupported, Molecula reports why instead of presenting an uncertain name as correct.\n\nDraw a supported structure and choose Help → Naming lesson.",
            (Self::Editor, true) => "Vigtige genveje: 1/2/3 vælger enkelt-, dobbelt- og tripelbinding; 4 vælger aromatisk binding; Esc skifter til markering; Delete sletter markeringen.\n\nTræk med markeringsværktøjet for at flytte eller markere. Klik midt på en eksisterende binding med bindingstasten 2 for at gøre den dobbelt.\n\nHelp → Shortcuts viser flere genveje. Kommandopaletten åbnes med Cmd/Ctrl+Shift+P.",
            (Self::Editor, false) => "Key shortcuts: 1/2/3 select single, double and triple bonds; 4 selects aromatic bond; Esc selects the selection tool; Delete removes the selection.\n\nDrag with the selection tool to move or select. Click the midpoint of an existing bond with bond key 2 to make it double.\n\nHelp → Shortcuts shows more. Open the command palette with Cmd/Ctrl+Shift+P.",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Reply {
    pub topic: Option<Topic>,
    pub title: String,
    pub body: String,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Exercise {
    pub title: &'static str,
    pub prompt: &'static str,
    pub hint_one: &'static str,
    pub hint_two: &'static str,
    pub worked_steps: &'static str,
}

pub const CHEMISTRY_A_EXERCISE_COUNT: usize = 7;

pub fn exercise(index: usize, da: bool) -> Option<Exercise> {
    Some(match (index, da) {
        (0, true) => Exercise {
            title: "Buffer-pH",
            prompt: "En buffer har pKa = 4,76, [HA] = 0,10 mol/L og [A⁻] = 0,20 mol/L. Beregn pH.",
            hint_one: "Identificér HA som syren og A⁻ som den korresponderende base.",
            hint_two: "Brug pH = pKa + log₁₀([A⁻]/[HA]).",
            worked_steps: "Forholdet [A⁻]/[HA] = 0,20/0,10 = 2. Derfor pH = 4,76 + log₁₀(2) ≈ 5,06.",
        },
        (1, true) => Exercise {
            title: "Svag syre",
            prompt: "En monoprot svag syre har pKa = 4,76 og startkoncentration 0,10 mol/L. Find pH med ligevægtsmodellen.",
            hint_one: "Skriv protolysen HA + H₂O ⇌ H₃O⁺ + A⁻.",
            hint_two: "Ka = 10⁻ᵖᴷᵃ. Ligevægten og vands autoprotolyse løses samlet.",
            worked_steps: "Kemi A-modellen løser ladningsbalancen for den monoprotte ligevægt og medregner Kw = 1,0×10⁻¹⁴ ved 25 °C. Se det numeriske ligevægtsresultat nedenfor.",
        },
        (2, true) => Exercise {
            title: "Fortynding",
            prompt: "Du overfører 10,0 mL af en 1,00 mol/L stamopløsning til en målekolbe og fylder til 100,0 mL. Find slutkoncentrationen i mol/L.",
            hint_one: "Stofmængden af opløst stof ændres ikke ved fortynding.",
            hint_two: "Brug c₁V₁ = c₂V₂ med samme volumen-enhed på begge sider.",
            worked_steps: "c₂ = c₁V₁/V₂ = 1,00 mol/L × 10,0 mL / 100,0 mL = 0,100 mol/L.",
        },
        (3, true) => Exercise {
            title: "Titrering før ækvivalenspunktet",
            prompt: "25,0 mL af en 0,100 mol/L monoprot syre med pKa 4,76 titreres med 0,100 mol/L stærk base. Der er tilsat 12,5 mL base. Beregn pH.",
            hint_one: "Beregn først stofmængden af HA og den tilsatte OH⁻.",
            hint_two: "Før ækvivalenspunktet er der både HA og A⁻; deres forhold bestemmer buffer-pH.",
            worked_steps: "Der er tilsat halvdelen af den nødvendige base til ækvivalens. I den ideelle monoprotte model er [A⁻] = [HA], så pH = pKa = 4,76.",
        },
        (0, false) => Exercise {
            title: "Buffer pH",
            prompt: "A buffer has pKa = 4.76, [HA] = 0.10 mol/L and [A⁻] = 0.20 mol/L. Calculate pH.",
            hint_one: "Identify HA as the acid and A⁻ as the conjugate base.",
            hint_two: "Use pH = pKa + log₁₀([A⁻]/[HA]).",
            worked_steps: "The ratio [A⁻]/[HA] = 0.20/0.10 = 2. Therefore pH = 4.76 + log₁₀(2) ≈ 5.06.",
        },
        (1, false) => Exercise {
            title: "Weak acid",
            prompt: "A monoprotic weak acid has pKa = 4.76 and initial concentration 0.10 mol/L. Find pH with the equilibrium model.",
            hint_one: "Write the protolysis as HA + H₂O ⇌ H₃O⁺ + A⁻.",
            hint_two: "Ka = 10⁻ᵖᴷᵃ. Solve the equilibrium together with water autoionization.",
            worked_steps: "The Chemistry A model solves charge balance for the monoprotic equilibrium and includes Kw = 1.0×10⁻¹⁴ at 25 °C. See the calculated equilibrium result below.",
        },
        (2, false) => Exercise {
            title: "Dilution",
            prompt: "Transfer 10.0 mL of a 1.00 mol/L stock solution to a volumetric flask and fill to 100.0 mL. Find the final concentration in mol/L.",
            hint_one: "The amount of dissolved substance is conserved during dilution.",
            hint_two: "Use c₁V₁ = c₂V₂ with the same volume unit on both sides.",
            worked_steps: "c₂ = c₁V₁/V₂ = 1.00 mol/L × 10.0 mL / 100.0 mL = 0.100 mol/L.",
        },
        (3, false) => Exercise {
            title: "Titration before equivalence",
            prompt: "Titrate 25.0 mL of 0.100 mol/L monoprotic acid, pKa 4.76, with 0.100 mol/L strong base. 12.5 mL base has been added. Calculate pH.",
            hint_one: "First calculate the amount of HA and added OH⁻.",
            hint_two: "Before equivalence, both HA and A⁻ are present; their ratio sets buffer pH.",
            worked_steps: "Half the base needed for equivalence has been added. In the ideal monoprotic model [A⁻] = [HA], so pH = pKa = 4.76.",
        },
        (4, true) => Exercise { title: "Begrænsende reaktant", prompt: "Reaktionen 2 H₂ + O₂ → 2 H₂O. Du har 0,30 mol H₂ og 0,10 mol O₂. Hvor meget H₂O kan højst dannes i mol?", hint_one: "Sammenlign n/ν for hver reaktant.", hint_two: "Reaktionsomfanget er min(0,30/2; 0,10/1). Gang med produktkoefficienten 2.", worked_steps: "O₂ begrænser reaktionen: ξ = 0,10 mol. Derfor dannes højst 2ξ = 0,20 mol H₂O." },
        (5, true) => Exercise { title: "Procentudbytte", prompt: "Det teoretiske udbytte er 0,200 mol H₂O (M = 18,015 g/mol). Der isoleres 3,20 g. Beregn procentudbyttet.", hint_one: "Beregn først den teoretiske masse fra n · M.", hint_two: "Procentudbytte = faktisk masse / teoretisk masse · 100 %.", worked_steps: "m(teoretisk) = 0,200 · 18,015 = 3,603 g. Udbytte = 3,20/3,603 · 100 % ≈ 88,8 %." },
        (6, true) => Exercise { title: "Betydende cifre", prompt: "Afrund 0,012345 g til tre betydende cifre. Indtast talværdien.", hint_one: "De tre første betydende cifre er 1, 2 og 3.", hint_two: "Det næste ciffer er 4, så det tredje ciffer forbliver uændret.", worked_steps: "0,012345 g afrundet til tre betydende cifre er 0,0123 g." },
        (4, false) => Exercise { title: "Limiting reactant", prompt: "For 2 H₂ + O₂ → 2 H₂O, you have 0.30 mol H₂ and 0.10 mol O₂. What maximum amount of H₂O forms, in mol?", hint_one: "Compare n/ν for each reactant.", hint_two: "The reaction extent is min(0.30/2; 0.10/1). Multiply by the product coefficient 2.", worked_steps: "O₂ is limiting: ξ = 0.10 mol. Therefore at most 2ξ = 0.20 mol H₂O forms." },
        (5, false) => Exercise { title: "Percent yield", prompt: "The theoretical yield is 0.200 mol H₂O (M = 18.015 g/mol). You isolate 3.20 g. Calculate percent yield.", hint_one: "First calculate theoretical mass from n · M.", hint_two: "Percent yield = actual mass / theoretical mass · 100%.", worked_steps: "m(theoretical) = 0.200 · 18.015 = 3.603 g. Yield = 3.20/3.603 · 100% ≈ 88.8%." },
        (6, false) => Exercise { title: "Significant figures", prompt: "Round 0.012345 g to three significant figures. Enter the numerical value.", hint_one: "The first three significant digits are 1, 2 and 3.", hint_two: "The next digit is 4, so the third digit stays unchanged.", worked_steps: "0.012345 g rounded to three significant figures is 0.0123 g." },
        _ => return None,
    })
}

fn exercise_expected(index: usize) -> Option<f64> {
    Some(match index {
        0 => crate::acid_base::buffer_ph(4.76, 0.10, 0.20).ok()?,
        1 => crate::acid_base::weak_acid_ph(4.76, 0.10).ok()?,
        2 => 1.00 * 0.0100 / 0.1000,
        3 => {
            crate::acid_base::polyprotic_titration_ph(&[4.76], 0.0250, 0.100, 0.100, 0.0125).ok()?
        }
        4 => 0.200,
        5 => 3.20 / (0.200 * 18.015) * 100.0,
        6 => 0.0123,
        _ => return None,
    })
}

pub fn exercise_hint(index: usize, level: u8, da: bool) -> Option<&'static str> {
    let exercise = exercise(index, da)?;
    Some(if level <= 1 {
        exercise.hint_one
    } else {
        exercise.hint_two
    })
}

pub fn search_topic(query: &str, da: bool) -> Option<Topic> {
    let query = query.trim().to_lowercase();
    if query.is_empty() {
        return None;
    }
    let terms: Vec<&str> = query
        .split(|ch: char| !ch.is_alphanumeric())
        .filter(|term| term.chars().count() >= 2)
        .collect();
    let mut best = None;
    let mut best_score = 0;
    for topic in Topic::ALL {
        let title = topic.title(da).to_lowercase();
        let article = topic.article(da).to_lowercase();
        let score = terms
            .iter()
            .map(|term| {
                (if title.contains(term) { 4 } else { 0 })
                    + (if article.contains(term) { 1 } else { 0 })
            })
            .sum();
        if score > best_score {
            best = Some(topic);
            best_score = score;
        }
    }
    best
}

pub fn check_exercise(index: usize, input: &str, da: bool) -> Option<(bool, String)> {
    let exercise = exercise(index, da)?;
    let expected = exercise_expected(index)?;
    let answer = match input.trim().replace(',', ".").parse::<f64>() {
        Ok(value) => value,
        Err(_) => {
            return Some((
                false,
                if da {
                    "Skriv et tal med decimal-komma eller punktum.".into()
                } else {
                    "Enter a number using a decimal comma or point.".into()
                },
            ));
        }
    };
    if !answer.is_finite() {
        return Some((
            false,
            if da {
                "Skriv et endeligt tal."
            } else {
                "Enter a finite number."
            }
            .into(),
        ));
    }
    let tolerance = if index == 0 || index == 1 || index == 3 {
        0.03
    } else if index == 5 {
        0.15
    } else if index == 6 {
        0.000_01
    } else {
        0.002
    };
    if (answer - expected).abs() <= tolerance {
        Some((
            true,
            format!(
                "{}\n\n{}\n\n{}: {expected:.4}",
                if da {
                    "Korrekt — svaret stemmer med den indbyggede beregningsmodel."
                } else {
                    "Correct — the answer matches the built-in calculation model."
                },
                exercise.worked_steps,
                if da {
                    "Beregnet svar"
                } else {
                    "Calculated answer"
                }
            ),
        ))
    } else {
        Some((
            false,
            format!(
                "{}\n\n{}\n\n{}: {expected:.4}",
                if da {
                    "Ikke helt. Kontrollér formel, forhold og enheder."
                } else {
                    "Not quite. Check the equation, ratio and units."
                },
                exercise.worked_steps,
                if da {
                    "Beregnet svar"
                } else {
                    "Calculated answer"
                }
            ),
        ))
    }
}

fn match_topic(question: &str) -> Option<Topic> {
    let q = question.to_lowercase();
    let rules: &[(Topic, &[&str])] = &[
        (
            Topic::Buffer,
            &[
                "buffer",
                "pka",
                "henderson",
                "[ha]",
                "a⁻",
                "korresponderende base",
                "conjugate base",
                "bufferkapacitet",
            ],
        ),
        (
            Topic::Titration,
            &[
                "titrer",
                "titration",
                "ækvivalens",
                "equivalence",
                "titrator",
            ],
        ),
        (
            Topic::AmountConcentration,
            &[
                "stofmængde",
                "molarmasse",
                "molarmass",
                "fortynd",
                "dilution",
                "koncentration",
                "concentration",
                "mol",
                "molaritet",
                "m =",
                "c =",
                "n =",
            ],
        ),
        (
            Topic::AcidBase,
            &[
                "syre",
                "base",
                "acid",
                "ph",
                "poh",
                "h3o",
                "oh-",
                "hydroxid",
                "oxonium",
                "ka",
                "kw",
                "ligevægt",
                "equilibrium",
            ],
        ),
        (
            Topic::Reaction,
            &[
                "reaktion",
                "reaction",
                "afstem",
                "balance",
                "koefficient",
                "coefficient",
            ],
        ),
        (
            Topic::Naming,
            &[
                "navn",
                "navngiv",
                "iupac",
                "naming",
                "aldehyd",
                "aldehyde",
                "heptan",
                "funktionel gruppe",
            ],
        ),
        (
            Topic::Editor,
            &[
                "binding",
                "dobbeltbinding",
                "double bond",
                "genvej",
                "shortcut",
                "tegn",
                "tegner",
                "draw",
                "editor",
                "molecula",
                "markering",
                "select",
            ],
        ),
    ];
    let mut best = None;
    let mut best_score = 0;
    for (topic, keywords) in rules {
        let score = keywords
            .iter()
            .filter(|keyword| {
                let keyword = **keyword;
                if keyword.chars().all(char::is_alphanumeric) && keyword.chars().count() <= 3 {
                    q.split(|ch: char| !ch.is_alphanumeric())
                        .any(|word| word == keyword)
                } else {
                    q.contains(keyword)
                }
            })
            .count();
        if score > best_score {
            best = Some(*topic);
            best_score = score;
        }
    }
    best.or_else(|| {
        if q.trim().is_empty() {
            None
        } else {
            Some(Topic::GettingStarted)
        }
    })
}

pub fn answer(question: &str, danish: bool) -> Reply {
    let q = question.to_lowercase();
    if (q.contains("pka") || q.contains("pka-værdi"))
        && [
            "tegnet",
            "tegning",
            "struktur",
            "molekyle",
            "drawn",
            "structure",
            "molecule",
        ]
        .iter()
        .any(|word| q.contains(word))
    {
        return Reply {
            topic: None,
            title: if danish { "pKa fra en tegnet struktur" } else { "pKa from a drawn structure" }.into(),
            body: if danish { "Molecula beregner ikke pKa ud fra den tegnede struktur. Brug en pKa-værdi fra undervisningsmateriale, en måling eller en pålidelig kilde. Kemi A-værktøjet kan beregne buffer-pH fra pKa og koncentrationer eller udlede pKa fra indtastede måledata med Henderson–Hasselbalch." } else { "Molecula does not calculate pKa from the drawn structure. Use a pKa value from class material, a measurement, or a reliable source. The Chemistry A tools can calculate buffer pH from pKa and concentrations, or derive pKa from measured values with Henderson–Hasselbalch." }.into(),
        };
    }
    let Some(topic) = match_topic(question) else {
        return Reply {
            topic: None,
            title: if danish { "Skriv et spørgsmål" } else { "Enter a question" }.into(),
            body: if danish { "Spørg fx om buffer-pH, pKa, titrering, stofmængde, reaktionsafstemning eller en editorfunktion." } else { "For example, ask about buffer pH, pKa, titration, amount of substance, reaction balancing, or an editor feature." }.into(),
        };
    };
    if topic == Topic::GettingStarted {
        return Reply {
            topic: Some(topic),
            title: if danish { "Det emne kan jeg ikke besvare sikkert endnu" } else { "I cannot answer that topic reliably yet" }.into(),
            body: if danish { "Jeg svarer kun ud fra det indbyggede, gennemgåede fagindhold. Prøv et af emnerne i Hjælpecenteret, eller spørg mere konkret om pH, buffere, titrering, stofmængde, reaktioner, navngivning eller editoren. Jeg sender ikke spørgsmålet til internettet." } else { "I only answer from the reviewed content built into the app. Try a topic in the Learning Center, or ask specifically about pH, buffers, titration, amount of substance, reactions, naming, or the editor. Your question is not sent to the internet." }.into(),
        };
    }
    Reply {
        topic: Some(topic),
        title: topic.title(danish).into(),
        body: topic.article(danish).into(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn assistant_matches_common_danish_and_english_questions() {
        assert_eq!(
            answer("Hvordan beregner jeg buffer pH?", true).topic,
            Some(Topic::Buffer)
        );
        assert_eq!(
            answer("How does titration equivalence work?", false).topic,
            Some(Topic::Titration)
        );
        assert_eq!(
            answer("hvordan finder jeg stofmængde fra masse", true).topic,
            Some(Topic::AmountConcentration)
        );
        assert_eq!(
            answer("How do I change a double bond?", false).topic,
            Some(Topic::Editor)
        );
        assert_eq!(
            answer("reaktionsafstemning", true).topic,
            Some(Topic::Reaction)
        );
    }

    #[test]
    fn assistant_refuses_to_invent_answers_outside_curated_topics() {
        let da = answer("Hvad bliver pKa for molekylet jeg har tegnet?", true);
        assert_eq!(da.topic, None);
        assert!(da.body.contains("beregner ikke pKa"));
        let en = answer("Will this drug cure my illness?", false);
        assert_eq!(en.topic, Some(Topic::GettingStarted));
        assert!(en.title.contains("cannot answer"));
    }

    #[test]
    fn curated_buffer_article_contains_verified_equation_and_boundary() {
        let r = answer("explain pKa and buffer", true);
        assert!(r.body.contains("pH = pKa + log₁₀([A⁻]/[HA])"));
        assert!(r.body.contains("udleder ikke pKa"));
    }

    #[test]
    fn chemistry_a_exercises_check_against_the_calculation_engines() {
        let (ok, _) = check_exercise(0, "5,061", true).unwrap();
        assert!(ok);
        let (ok, _) = check_exercise(0, "7.2", false).unwrap();
        assert!(!ok);
        let (ok, _) =
            check_exercise(1, &format!("{:.3}", exercise_expected(1).unwrap()), false).unwrap();
        assert!(ok);
        let (ok, _) = check_exercise(2, "0,100", true).unwrap();
        assert!(ok);
        assert!(check_exercise(4, "0,2", true).unwrap().0);
        assert!(check_exercise(5, "88,8", true).unwrap().0);
        assert!(check_exercise(6, "0,0123", true).unwrap().0);
        let (ok, feedback) = check_exercise(3, "4.76", false).unwrap();
        assert!(ok, "{feedback}");
        assert!(check_exercise(2, "not a number", false)
            .unwrap()
            .1
            .contains("Enter a number"));
        assert_eq!(exercise(CHEMISTRY_A_EXERCISE_COUNT, true), None);
    }

    #[test]
    fn learning_center_search_matches_curated_danish_and_english_articles() {
        assert_eq!(search_topic("buffer pKa", true), Some(Topic::Buffer));
        assert_eq!(
            search_topic("titration equivalence", false),
            Some(Topic::Titration)
        );
        assert_eq!(search_topic("molekyle-editor", true), Some(Topic::Editor));
        assert_eq!(search_topic("", false), None);
        assert_eq!(search_topic("xyzqv", false), None);
    }
}
