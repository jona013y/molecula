# Changelog

## v1.0.15

- Add a Chemistry A limiting-reactant and yield calculator for theoretical product amount/mass and optional percent yield from actual mass.
- Add configurable 3–5 significant-figure display precision; basic and acid/base reports now use the selected precision with decimal-comma Danish display.
- Expand the interactive Chemistry A practice set from four to seven exercises, adding limiting reagent, percent yield and significant figures.
- Link the drawn reaction analysis guide to yield calculations while requiring known amounts to be entered separately.
- Keep the existing pKa, polyprotic distribution/Bjerrum charts and interactive titration curves, which were already present in v1.0.14.
- Verify calculation equations, input validation, settings persistence, exercise grading and the full local test suite.

## v1.0.14

- Add a searchable Danish/English learning center with eight reviewed, bundled Chemistry A and editor guides.
- Add an offline local chemistry helper with bounded intent matching and an explicit fallback for unsupported questions; it does not use a generative model, network service, or structure-based pKa prediction.
- Add Chemistry A calculators for amount, mass, concentration, dilution, and pH/pOH conversions; show equations, substitutions, units, and validated answers.
- Add four Chemistry A practice problems with progressive hints, checked answers, and locally persisted progress.
- Verify the full Rust suite (231 tests), formatting, Clippy with warnings denied, and optimized release build.

## v1.0.13

- Add a Chemistry A acid/base workbench: Henderson–Hasselbalch buffer pH and measured pKa, weak monoprotic acid/base equilibrium pH, strong acid/base mixtures, polyprotic speciation, Bjerrum diagrams and titration calculations.
- Plot Bjerrum distributions and polyprotic titration curves in the app, with CSV copy for plotted points. Localize reports and input examples in Danish; state the ideal aqueous 25 °C assumptions and require user-supplied pKa values instead of inventing structure-based predictions.
- Verify acid/base equations, edge cases, Danish outputs and plotted data; the complete Rust test suite passes.

## Unreleased

- Extend the optional RDKit descriptor report with deterministic-conformer solvent-accessible surface area using rdFreeSASA Shrake-Rupley; include the probe, force field, radii and conformer method, and report when the engine cannot calculate it.
- Match selected substructure queries against explicitly stored R/S atom labels and E/Z/cis/trans double-bond labels. The matcher still refuses stereo groups and wedge/hash stereochemistry.
- Add a bounded terminal nitrile naming path, including unsaturation in the parent chain (ethanenitrile, prop-2-enenitrile); count the nitrile carbon in the parent and refuse mixed nitrile/principal-functional-group structures.

- Extend optional AtomMapper reaction mapping from one component pair to all explicitly selected reactant/product components in one reaction step, sent together as disconnected reaction sides; refuse agents, partial components, mixed roles and cross-step selections, and apply the verified mapping in one undoable edit.

- Extend bounded cyclic naming to cover hydrocarbon dienes and up to three noncumulated ring double bonds with a suffix alcohol or ketone, including principal-group numbering and locant selection. Verify cyclohexa-1,3-diene, cyclohexa-2,4-dien-1-ol and cyclohexa-2,4-dien-1-one.

- Preserve RDKit-imported E/Z alkene assignments in the molecular graph and carry them through native-document save/reopen and chemistry-engine SMILES export; explicit E/Z annotations now override conflicting 2D layout during export. Verify cis/trans alkenes, layout conflicts and tetrahedral chirality end to end.

- Add optional RDKit Gasteiger-Marsili partial-charge calculation with explicit per-atom and implicit-hydrogen charge output, method, iteration count, and net formal-charge check; run in background without altering the molecule.
- Add a validated, undoable custom reaction stoichiometry coefficient (1–99) for the selected connected component; localize the menu and entry dialog in Danish.
- Add optional automatic reaction atom mapping with AtomMapper + RDKit in Complete, Matching and Changing modes for an explicitly selected reactant/product component pair; run in background, validate results and apply as one undoable edit; reject duplicate correspondences and preserve existing maps when a mode yields no assignments. Serialize background chemistry jobs so results cannot be overwritten by another calculation.
- Document the optional mapping engine requirements, algorithm limitations and single-component-pair scope; localize its commands and quick-start guidance in Danish.
- Localize the quick-start guide in Danish or English according to the saved interface setting.
- Explain the one-click workflow for changing an existing bond to a double bond and undoing the edit.
- Keep bond-order queries consistent when a query bond is changed in place; preserve independent ring/chain constraints and make query-bond tools edit the existing edge.

## v1.0.12

- Add **Edit → Copy drawing as image** to copy a 1920×1080 PNG through native macOS, Windows and X11/Wayland clipboard backends; surface unavailable-tool errors.
- Render normal, aromatic, filled-wedge, hashed-wedge, wavy, coordinate and query-any bonds correctly in SVG, PNG and PDF exports; verify exported geometry for every style.

- Write V2000 atom-map fields at the correct columns, emit a valid V3000 header and omit empty V3000 bond blocks. Normalize RDKit explicit map-zero values so SMILES exports omit `:0` while preserving assigned maps; add real-RDKit V2000/V3000 round-trip regressions.

- Exclude the stale prebuilt app executable from the public source archive and remap Rust build-home paths before packaging to avoid publishing local build paths.

- Add bounded local CML import/export for one molecule, preserving supported graph fields, coordinates and simple metadata. Reject unsupported CML annotations and warn about export losses; CML no longer depends on optional RDKit.

- Add direct PubChem structure lookup through the official NCBI PUG REST API, separate from opening a browser search.
- Request PubChem same-stereochemistry/isotope identity matching for locally representable SMILES, then show database properties with the PubChem CID and source link. Refuse local stereochemical annotations that the SMILES writer cannot preserve.
- Report PubChem network failures as PubChem errors, and keep failed background jobs from showing a completed status.
- Add depth-ordered 3D spacefill and monochrome rendering modes, expand common-element CPK colors, avoid drawing atoms with missing 3D coordinates, and add a camera reset plus a persisted atom-symbol toggle.
- Keep the drawing unchanged and refuse structures that cannot be represented faithfully by Molecula's local SMILES writer.
- Add a direct “Make double” action to the selected-bond inspector. Select a bond and set it to a double bond without redrawing it; the edit remains undoable.
- Persist the canvas zoom preference between sessions and reuse it for new documents; clamp invalid stored values and preserve other settings.
- Name supported neutral primary amine substituents once as `amino`; refuse amides and uncovered heteroatom groups instead of returning a misleading partial name.

## v1.0.11

- Extend the bounded teaching namer to recognize common branched saturated alkyl substituents, including propan-2-yl and 2-methylpropan-1-yl.
- Correct saturated aldehyde parent formation so a heptane-derived aldehyde is named heptanal.
- Keep unsupported functionalized branches as explicit refusals instead of generating incomplete names.

## v1.0.10

- Let the selected single, double, triple or styled-bond tool edit an existing bond with one click on its midpoint; atom clicks still start the usual bond drawing gesture.
- Keep the in-place edit undoable and select the changed bond for immediate inspection.

## v1.0.9

- Prevent stale E/Z annotations when a double bond is changed to an aromatic or single bond.
- Keep wedge, hash and wavy stereobond display marks on single bonds; converting a double bond to that style is undoable.
- Report imported bond stereochemistry/display-style combinations that conflict with their bond order.

## v1.0.8

- Make reaction analysis check isotope counts separately from element totals and formal charge.
- Require isotope-labelled atoms to be conserved by automatic stoichiometric balancing; report isotope differences explicitly.
- Add regression coverage for isotope loss and isotope-preserving reactions.

## v1.0.7

- Add a per-component reaction table to the analysis report, listing role, coefficient, formula, formal charge and whether the role was explicitly assigned or inferred from geometry.
- Cover explicit/inferred role labels and coefficient/formula reporting with regression tests.

## v1.0.6

- Include per-step reaction conditions in the calculation report, showing the text above and below each arrow.
- Keep conditions isolated to their own step and omit the conditions section when both arrow labels are blank.

## v1.0.5

- Add atomic automatic balancing for every reaction step. The editor solves each step before applying any coefficient changes, reports the result by step, and records the complete operation as one undoable edit.
- Add regression coverage for multi-step balancing, undo, and failure without partial edits.

## v1.0.4

- Add searchable, graph-backed templates for naphthalene, indole, quinoline, anthracene, phenanthrene, decalin and norbornane.
- Confirm detected data loss before writing SMILES, SMARTS, image and exchange-format exports.
- Add a bounded local SMILES import/export path for supported structures without requiring RDKit.
- Refresh the Mac application package and improve source-based Mac/Windows CI packaging.

Molecula remains an independent editor with a limited, documented chemistry scope; this release does not claim feature parity with MarvinSketch or completeness for research-grade cheminformatics.
