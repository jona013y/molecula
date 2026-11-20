# Molecula — roadmap for approval

All packages A–G were approved by Jonathan (“add all”). They are being implemented incrementally. A checked item means implemented to the limits in `STATUS.md`; it does not mean MarvinSketch feature parity.

## Previous approved work

- [x] Editor interaction: zoom, pan, multi-select, move/rotate/duplicate/delete, undo/redo.
- [x] Reaction model: manual component roles, step numbers, coefficients and balance report.
- [x] Auto-balance single-step reactions when the positive integer solution is unique; include formal charge and reject unresolved/ambiguous cases.
- [x] Teaching name explanation with supported-scope warnings and highlighted atoms.
- [x] Extend the bounded teaching namer with carboxylic-acid priority, saturated/unsaturated parent names, and simple cycloalkanes/cycloalkenes with straight alkyl/halo substituents, cyclic ketones and cyclic alcohols, including ring dienes and one to three noncumulated ring double bonds alongside the suffix group; test priority, numbering and prefixes.
- [x] Name neutral straight-chain alkyl alkanoate esters with acyl-parent and alkoxy-chain identification; refuse branched, multiply esterified, or otherwise unsupported cases and test both accepted names and refusal.
- [x] Add bounded naming for terminal nitriles, counting the nitrile carbon in the parent and ranking the nitrile suffix above aldehyde, ketone and alcohol; refuse mixed principal groups until their prefix forms are supported.
- [x] Rank equal-length acyclic parent-chain orientations by principal-group locant, multiple-bond locants, then substituent locants; include locant comparison in the learner explanation.
- [x] Name supported saturated acyclic diols and diketones with multiple principal suffix locants, requiring a parent chain that contains every suffix group.
- [x] Name common saturated branched alkyl substituents such as propan-2-yl and 2-methylpropan-1-yl; preserve explicit refusal for functionalized branches. Correct saturated aldehyde suffix construction and cover the requested heptanal naming workflow.
- [x] Name neutral primary amine substituents with the amino prefix, and refuse amides or unsupported heteroatom groups rather than mislabeling them as aldehydes or simple prefixes.
- [x] Guided naming/reaction balancing lessons with real graph-based feedback.
- [x] Local graph topology in the properties report; optional RDKit descriptor/candidate-analysis bridge.
- [x] Document flow: tabs, autosave/recovery, recent files, tab-local undo.
- [x] Optional RDKit bridge and MOL/SDF/SMILES paths, with limits described in `STATUS.md`.

## Approved feature packages — in progress

### A. Safe document workflow and everyday polish — core implemented

- Close a tab with a save/discard/cancel choice; ask before exiting with unsaved tabs.
- Add a recovery chooser when autosaved content is found, keeping originals untouched.
- Add dark/light theme, saved zoom/appearance preferences, show/hide implicit hydrogens, and a shortcut reference that reflects the active platform.
- [x] Persist the current canvas zoom as the default for new documents and later sessions; clamp corrupt or out-of-range stored values and preserve other settings.
- Show a compact zoom percentage and tool/selection-specific status.

**Acceptance:** scripted tests cover all close/save/discard branches and recovery never overwrites source files; UI check confirms theme and shortcut state persists after restart.

**Dependencies:** no chemistry engine.

### B. Reaction studio — partial

- [x] Add an inspector for assigning a connected component as reactant, product or agent and editing its step and coefficient, including custom integer coefficients from 1 to 99 with validation and undo.
- [x] Map all explicitly selected reactant/product components in a single step as disconnected reaction sides; leave agents out, require complete components and one shared explicit reaction step, and apply atom maps as one undoable operation.
- [x] Add atom-map editing and a visible reaction-change report that highlights formed, broken and order-changed bonds on the canvas; `M+`/`M=` assigns collision-free map numbers to selected atoms or all atoms when no selection exists.
- [x] Solve and apply unique single-step stoichiometric coefficients with undo support.
- [x] Navigate reaction steps beyond the first three and restrict reaction analysis/highlighting to the displayed step.
- [x] Include each reaction arrow's above/below conditions in its corresponding per-step analysis report; test multi-step condition isolation and blank-condition behavior.
- [x] Add an explicit reaction component role/stoichiometry balance table, including explicit/inferred role source and each component's coefficient, formula and formal charge.
- [x] Preserve isotope counts during reaction analysis and automatic coefficient balancing.
- [x] Add optional local automatic atom mapping with AtomMapper + RDKit for all fully selected, explicitly assigned reactant/product components in one reaction step, presented as disconnected reaction sides. Offer Complete, Matching and Changing modes; run off the UI thread, apply maps as one undoable edit, serialize background chemistry jobs, and leave the graph unchanged on unavailable, incomplete or malformed engine results.

**Acceptance:** multi-step documents survive save/reopen; role and coefficient edits undo correctly; reports never silently infer atom correspondences.

**Dependencies:** manual tools require none; automatic multi-component-side mapping requires Python 3.11+, RDKit and AtomMapper (which also installs NetworkX and SciPy). AtomMapper returns a shared core across the disconnected sides; `complete` means search completion, not proof of the globally best mapping or a separately scored one-to-one component assignment.

### C. Structure editor and template library — partial

- [x] Add searchable built-in and user-template tabs; save selected connected fragments as named native `.mols` templates, reload them between sessions, insert them as real atom/bond graphs, rename while preserving graph data, reject duplicate names and confirm deletion.
- [x] Add a pointer-following translucent preview for built-in and custom templates; rotate by 60° with `R` before placement and allow repeated placement of saved user templates.
- [x] Attach a template to exactly one selected atom with `A`; preview an adjacent anchor atom, create a real single-bond graph connection, and undo the whole insertion as one command. Multiple named attachment sites remain future work.
- Add broader categorized built-in templates with common functional groups, fused/bridged systems and biomolecular building blocks.
- [x] Add searchable naphthalene, indole, quinoline, anthracene and phenanthrene aromatic fused-ring graphs; add decalin and bicyclo[2.2.1]heptane (norbornane) saturated fused/bridged graphs. Verify formula, connectivity, cycle rank, valence warnings, aromaticity and single-step undo.
- [x] Add per-atom valence override editing, persist it in native documents, use it for implicit-H estimates and preserve it in supported MOL V2000/V3000 exports.
- [x] Add undoable bond order and display-style controls plus a manual E/Z toggle for selected double bonds.
- [x] Let an active bond tool change an existing bond on a midpoint click, while clicks on atoms continue to start a new bond; make the in-place edit selectable and undoable.
- [x] Add a direct selected-bond inspector action to make the clicked bond double without redrawing it; keep the graph edit undoable and localize the control in Danish.
- [x] Parse and combine supported atom-query shortcuts (aromaticity, H counts, connectivity, degree, ring data, valence, substitutions and element include/exclude lists), rejecting malformed or contradictory same-field predicates.
- [x] Add a bounded local substructure matcher for selected whole query molecules against other disconnected target components; evaluate supported element/list, aromaticity, unsaturation, H-count, degree/connectivity, valence, bond-order, substitution count/as-drawn, incident ring-bond count, exact minimum-ring-size, deterministic minimum-cycle-basis ring counts, non-ring-atom predicates and exact stored R/S and E/Z/cis/trans labels. Reject unsupported stereo groups and partial query selections.
- [x] Add a selected-atom query inspector for aromatic, aliphatic, unsaturated, any-atom and heteroatom predicates, with a clear-query action and undo support.
- [x] Add direct selected-atom cycling for supported singlet, doublet and triplet radical states, including undo/redo.
- [x] Add typed selected-atom isotope and alias editing with a complete-token buffer, alias rendering on carbon, clear actions and undo/redo; preserve isotope data through supported MOL export/import.
- [x] Add direct atom-inspector cycling for automatic or explicit implicit-hydrogen counts H0–H5, with undo/redo and MOL V2000/V3000 persistence.
- [x] Apply explicitly requested RDKit CIP R/S assignments as undoable atom annotations; show labels on the canvas and refuse malformed/duplicate worker output.
- Still needed: broader stereochemical matching beyond exact stored CIP and double-bond labels (including wedge/hash configuration, stereogroups and geometric verification), and toolkit-specific ring semantics. RDKit assignment remains optional and requires end-to-end verification where installed.
- Add fragment alignment and keyboard-driven bond editing with an explicit preview.
- Add aromatic/Kekulé conversion only when a chemistry engine can validate the resulting graph.
- [x] Add formula-checked aromatic pyrrole, furan, thiophene and correctly protonated imidazole templates.
- [x] Add formula-checked graph templates for all 20 standard proteinogenic amino acids; keep their alpha-carbon stereochemistry unspecified instead of implying L/D configuration.

**Acceptance:** template insertion creates chemically connected graph data; invalid chemistry remains editable but is visibly flagged; every edit is undoable.

**Dependencies:** core template/editor work needs none; aromaticity perception may require optional RDKit.

### D. 3D structures and conformers — partial; optional engine

- Generate 3D coordinates for supported molecules through optional RDKit; show the force field and optimization status.
- [x] Add optional solvent-accessible surface area to the RDKit report using a deterministic ETKDGv3 conformer, MMFF94/UFF minimization, RDKit van der Waals radii, Shrake-Rupley, and a 1.4 Å probe; report the method and an unavailable reason when the build lacks rdFreeSASA or a usable force field.
- [x] Add a rotatable 3D view with wireframe, sticks, ball-and-stick, spacefill and monochrome styles; depth-order atoms, use a common-element CPK palette and common-element van der Waals radii (fallback for other elements), skip/report atoms without 3D coordinates, add a reset-view control, persist the atom-symbol visibility toggle, and keep camera movement separate from the 2D drawing.
- Keep 3D generation separate from 2D stereochemical identity; failures must preserve the source structure.

**Acceptance:** test geometry generation and failure handling with/without RDKit; camera rotation cannot change the stored 2D graph or stereochemistry.

**Dependencies:** separate Python + RDKit installation; package licensing and installation instructions need review.

### E. File interoperability and export — partial

- Add V3000 MOL/SDF support for structures larger than V2000 limits, retaining supported atom/bond data.
- Add RXN import/export for reaction graphs, with a before-export loss report.
- [x] Export drawing documents to SVG, PNG and PDF with rendered bond display styles (including aromatic, wedge, hash, wavy and coordinate bonds).
- [x] Add a 1920×1080 PNG copy-to-clipboard command using the native macOS, Windows and X11/Wayland clipboard paths, with explicit platform errors.
- [x] Preserve and validate MOL V2000/V3000 singlet, doublet and triplet radical records; round-trip tests cover valid multi-atom records and invalid radical types.
- [x] Preserve MOL V2000 hydrogen-count overrides (H0–H4), V3000 `HCOUNT` query values (H0–H5), and the exact-count `SUBST`, `UNSAT` and `RBCNT` atom-query fields; warn before export when a format or unsupported value would drop data.
- [x] Preserve V3000 bond query types 5–8 and ring/chain `TOPO=1/2`; exercise their round trips and ring-aware substructure matching.
- [x] Add a selected-bond inspector control to cycle through unconstrained, ring-only and chain-only query topology; preserve other bond predicates and support undo/redo.
- [x] Clear E/Z annotations when edited bond order becomes non-double; constrain wedge/hash/wavy display marks to single bonds and warn about conflicting imported stereo/order data.
- [x] Preserve V3000 wildcard/heteroatom selectors and included/excluded element lists; evaluate them in the local substructure matcher and reject malformed lists.
- [x] Preserve MOL V2000 parity 3 and V3000 atom CFG=3 as explicitly unspecified/either stereochemistry through round trips.
- [x] Export a bounded, local SMARTS subset to clipboard or `.smarts`, using graph traversal for branches/rings and refusing unsupported predicates or annotations rather than dropping them.
- [x] Import the supported SMARTS subset from `.smarts`; reject unsupported recursive/stereo syntax and malformed patterns instead of guessing, with graph-aware chain/ring/fragment layout and parser/exporter round-trip tests.
- [x] Add bounded local CML import/export for one molecule, including supported graph data, coordinates and simple metadata; refuse unsupported import annotations and warn before lossy export.
- [x] Add a local bounded SMILES reader/writer for ordinary and aromatic atoms, branches, rings, fragments, isotope/charge/H/map brackets; route `.smi`/`.smiles` and clipboard through it, with optional RDKit fallback and explicit refusal of unsupported stereo/query features.
- [x] Restore aromatic atom flags from imported type-4 aromatic bonds and avoid MOL/SDF loss warnings for supported atom stereo, wedge/hash/wavy, aromatic bonds, maps, isotopes, charges and radicals; keep warnings for unsupported query/alias/atom extensions.
- [x] Add editable molecule data fields (`KEY=VALUE`, `-KEY` removal and molecule `NAME`) with undo and native/SDF persistence checks.
- [x] Add `Cmd/Ctrl+L` copy as SMILES through the local writer when representable, with optional RDKit fallback and platform clipboard; refuse unsupported alias/query/R-group and nonstandard bond data.
- Add robust malformed-file messages and round-trip tests for every shipped format.

**Acceptance:** round trips preserve every supported field; unsupported fields are listed before export; no conversion reports success after dropping data silently.

**Dependencies:** image/PDF export depends on the selected rendering backend; RXN and V3000 need format-specific verification.

### F. Teaching mode — partial

- [x] Add a searchable Danish/English learning center with bundled, reviewed guides and a context menu entry.
- [x] Add an offline local chemistry helper that answers a bounded set of reviewed topics and clearly declines unsupported or structure-based pKa questions; no network, account, API key or generative model is used.
- [x] Add Kemi A amount/mass/concentration/dilution and pH/pOH calculators with equations, substitutions, units and validation; retain the existing pKa and equilibrium workbench.
- [x] Retain pKa from measured buffer data, polyprotic speciation/Bjerrum diagrams and interactive charge-balance titration curves; reports require supplied pKa values and state ideal-solution conditions.
- [x] Add limiting-reactant, theoretical-yield and percent-yield calculations from a balanced equation's entered coefficients and measured amounts; amounts are not inferred from the drawing.
- [x] Add 3–5 significant-figure display precision to basic and acid/base calculation reports, with persisted preference, decimal-comma Danish output, and unrounded internal calculations.
- [x] Add seven interactive numerical Kemi A exercises with progressive hints, answer grading through shared calculation routines, and locally persisted solved progress.
- [x] Link reaction-analysis guidance to the yield calculator while explicitly asking learners to enter quantities separately from the drawn reaction.
- Add a guided, interactive lesson flow for drawing common structures and reactions.
- Explain naming steps on the canvas and let the learner reveal functional-group priority, parent chain, numbering and prefixes one at a time.
- Add reaction balancing practice with immediate, evidence-based feedback; do not grade unsupported structures.

**Acceptance:** lessons have deterministic expected outcomes, are keyboard/mouse accessible, and unsupported chemistry produces a clear “not evaluated” result.

**Dependencies:** no external service for the initial lessons.

**Current implementation boundary:** the helper uses a small keyword/intent matcher over bundled reviewed content. It is not a generative AI, has no access to the drawn structure and does not predict structure-based pKa or infer quantities from a drawing. The guides and seven exercises are a starter set, not a complete Chemistry A curriculum.

### G. Chemistry lookup and property reports — partial

- [x] Extend local formula/mass reporting with CIAAW 2024 standard weights (using the midpoint for interval values), explicit isotope masses, and exact-mass correction for ionic electron gain/loss; keep unsupported elements unavailable.
- [x] Add optional RDKit descriptor and Gasteiger-Marsili partial-charge reports with method, units, engine version, explicit/implicit-H charge separation and total-charge validation. Runtime smoke checks use a temporary RDKit environment; the app remains independent of it.
- [x] Add local classical heavy-atom graph descriptors (Wiener index, diameter, first/second Zagreb indices, and Randić connectivity), with explicit graph definition and refusal for disconnected/query/nonstandard structures.
- [x] Provide bounded, documented tautomer and unspecified-stereoisomer candidate enumeration without changing the source graph; optional RDKit required and runtime verification is pending.
- [x] Keep browser search distinct from a separate direct PubChem PUG REST lookup; request exact structure identity including stereo/isotopes, fetch linked property records, and never modify the drawing.
- Do not promise pKa, NMR, logD, solubility or IUPAC names until a real suitable engine/service and its license/data terms are verified.

**Acceptance:** each displayed result names its source/method; unavailable properties stay unavailable; network lookup can be disabled and never requires an embedded API key.

**Dependencies:** RDKit for local descriptors; optional internet for direct database lookup.

## Current boundary

All packages are authorized, but remain incomplete until their acceptance criteria are met. Current implemented scope and remaining limitations are tracked in `STATUS.md`. Continue from this codebase; do not treat this roadmap as a claim that every item is complete.
