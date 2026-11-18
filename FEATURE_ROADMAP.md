# Molecula — roadmap for approval

All packages A–G were approved by Jonathan (“add all”). They are being implemented incrementally. A checked item means implemented to the limits in `STATUS.md`; it does not mean MarvinSketch feature parity.

## Previous approved work

- [x] Editor interaction: zoom, pan, multi-select, move/rotate/duplicate/delete, undo/redo.
- [x] Reaction model: manual component roles, step numbers, coefficients and balance report.
- [x] Auto-balance single-step reactions when the positive integer solution is unique; include formal charge and reject unresolved/ambiguous cases.
- [x] Teaching name explanation with supported-scope warnings and highlighted atoms.
- [x] Extend the bounded teaching namer with carboxylic-acid priority and saturated/unsaturated parent names; test acid numbering and prefixes.
- [x] Guided naming/reaction balancing lessons with real graph-based feedback.
- [x] Local graph topology in the properties report; optional RDKit descriptor/candidate-analysis bridge.
- [x] Document flow: tabs, autosave/recovery, recent files, tab-local undo.
- [x] Optional RDKit bridge and MOL/SDF/SMILES paths, with limits described in `STATUS.md`.

## Approved feature packages — in progress

### A. Safe document workflow and everyday polish — core implemented

- Close a tab with a save/discard/cancel choice; ask before exiting with unsaved tabs.
- Add a recovery chooser when autosaved content is found, keeping originals untouched.
- Add dark/light theme, saved zoom/appearance preferences, and a shortcut reference that reflects the active platform.
- Show a compact zoom percentage and tool/selection-specific status.

**Acceptance:** scripted tests cover all close/save/discard branches and recovery never overwrites source files; UI check confirms theme and shortcut state persists after restart.

**Dependencies:** no chemistry engine.

### B. Reaction studio — partial

- Add an inspector for assigning an entire connected component as reactant, product or agent, and editing its step and coefficient.
- Add atom-map editing and a visible reaction-change report that highlights formed, broken and order-changed bonds on the canvas.
- [x] Solve and apply unique single-step stoichiometric coefficients with undo support.
- Add reaction-step navigation with per-step conditions, plus an explicit role/stoichiometry balance table.
- Keep automatic atom mapping optional and unavailable unless a verified local engine is installed.

**Acceptance:** multi-step documents survive save/reopen; role and coefficient edits undo correctly; reports never silently infer atom correspondences.

**Dependencies:** manual tools require none; automatic mapping requires a separately verified engine.

### C. Structure editor and template library — partial

- Add a searchable, categorized template browser with common functional groups, aromatic/heterocyclic rings and amino acids.
- [x] Add per-atom valence override editing, persist it in native documents, use it for implicit-H estimates and preserve it in supported MOL V2000/V3000 exports.
- [x] Add undoable bond order and display-style controls plus a manual E/Z toggle for selected double bonds.
- [x] Parse and combine supported atom-query shortcuts (aromaticity, H counts, connectivity, degree, ring data, valence, substitutions and element include/exclude lists), rejecting malformed or contradictory same-field predicates.
- [x] Add a bounded local substructure matcher for selected whole query molecules against other disconnected target components; evaluate supported element/list, aromaticity, unsaturation, H-count, degree/connectivity, valence, bond-order, substitution count/as-drawn, incident ring-bond count, exact minimum-ring-size, deterministic minimum-cycle-basis ring counts, and non-ring-atom predicates. Refuse stereochemical predicates explicitly; reject partial query selections.
- Still needed: query-atom inspector and broader matching (including toolkit-specific ring semantics and stereochemistry), plus stereobond configuration and geometric stereochemistry validation.
- Add fragment alignment/rotation and keyboard-driven bond editing with an explicit preview.
- Add aromatic/Kekulé conversion only when a chemistry engine can validate the resulting graph.
- [x] Add formula-checked aromatic pyrrole, furan, thiophene and correctly protonated imidazole templates.

**Acceptance:** template insertion creates chemically connected graph data; invalid chemistry remains editable but is visibly flagged; every edit is undoable.

**Dependencies:** core template/editor work needs none; aromaticity perception may require optional RDKit.

### D. 3D structures and conformers — partial; optional engine

- Generate 3D coordinates for supported molecules through optional RDKit; show the force field and optimization status.
- Add a rotatable 3D view with wireframe, sticks and ball-and-stick styles, zoom/reset, and a control to return to the original 2D drawing.
- Keep 3D generation separate from 2D stereochemical identity; failures must preserve the source structure.

**Acceptance:** test geometry generation and failure handling with/without RDKit; camera rotation cannot change the stored 2D graph or stereochemistry.

**Dependencies:** separate Python + RDKit installation; package licensing and installation instructions need review.

### E. File interoperability and export — partial

- Add V3000 MOL/SDF support for structures larger than V2000 limits, retaining supported atom/bond data.
- Add RXN import/export for reaction graphs, with a before-export loss report.
- Add PNG and PDF drawing export and clipboard copy as SVG/image where the platform supports it.
- [x] Preserve and validate MOL V2000/V3000 singlet, doublet and triplet radical records; round-trip tests cover valid multi-atom records and invalid radical types.
- [x] Preserve MOL V2000 hydrogen-count overrides (H0–H4), V3000 `HCOUNT` query values (H0–H5), and the exact-count `SUBST`, `UNSAT` and `RBCNT` atom-query fields; warn before export when a format or unsupported value would drop data.
- [x] Preserve V3000 bond query types 5–8 and ring-only `TOPO=1`; exercise their round trips and ring-aware substructure matching.
- [x] Preserve V3000 wildcard/heteroatom selectors and included/excluded element lists; evaluate them in the local substructure matcher and reject malformed lists.
- [x] Preserve MOL V2000 parity 3 and V3000 atom CFG=3 as explicitly unspecified/either stereochemistry through round trips.
- [x] Restore aromatic atom flags from imported type-4 aromatic bonds and avoid MOL/SDF loss warnings for supported atom stereo, wedge/hash/wavy, aromatic bonds, maps, isotopes, charges and radicals; keep warnings for unsupported query/alias/atom extensions.
- Add robust malformed-file messages and round-trip tests for every shipped format.

**Acceptance:** round trips preserve every supported field; unsupported fields are listed before export; no conversion reports success after dropping data silently.

**Dependencies:** image/PDF export depends on the selected rendering backend; RXN and V3000 need format-specific verification.

### F. Teaching mode — partial

- Add a guided, interactive lesson flow for drawing common structures and reactions.
- Explain naming steps on the canvas and let the learner reveal functional-group priority, parent chain, numbering and prefixes one at a time.
- Add reaction balancing practice with immediate, evidence-based feedback; do not grade unsupported structures.

**Acceptance:** lessons have deterministic expected outcomes, are keyboard/mouse accessible, and unsupported chemistry produces a clear “not evaluated” result.

**Dependencies:** no external service for the initial lessons.

### G. Chemistry lookup and property reports — partial

- [x] Extend local formula/mass reporting with CIAAW 2024 standard weights (using the midpoint for interval values), explicit isotope masses, and exact-mass correction for ionic electron gain/loss; keep unsupported elements unavailable.
- [x] Add a local descriptor report for structures accepted by the optional RDKit bridge, including method, units and engine version; this remains unverified at runtime without RDKit installed.
- [x] Provide bounded, documented tautomer and unspecified-stereoisomer candidate enumeration without changing the source graph; optional RDKit required and runtime verification is pending.
- Offer PubChem lookup as a clearly labeled browser search first; add direct API retrieval only as a separate opt-in integration.
- Do not promise pKa, NMR, logD, solubility or IUPAC names until a real suitable engine/service and its license/data terms are verified.

**Acceptance:** each displayed result names its source/method; unavailable properties stay unavailable; network lookup can be disabled and never requires an embedded API key.

**Dependencies:** RDKit for local descriptors; optional internet for direct database lookup.

## Current boundary

All packages are authorized, but remain incomplete until their acceptance criteria are met. Current implemented scope and remaining limitations are tracked in `STATUS.md`. Continue from this codebase; do not treat this roadmap as a claim that every item is complete.
