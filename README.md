# Molecula

Molecula is an independent Rust desktop editor for chemical structures and reaction diagrams. It includes a molecular graph editor, selection and transforms, zoom/pan, reaction roles and stoichiometry, a guided naming subset, local formula/mass analysis, multi-document tabs, crash recovery, MOL V2000/V3000, SDF, RXN, SVG, PNG and PDF support, optional RDKit operations, a rotatable 3D viewer, and a searchable command palette (`Cmd/Ctrl+Shift+P`). PubChem lookup opens a clearly labeled browser search. It is not a MarvinSketch clone or a complete cheminformatics package. See [STATUS.md](STATUS.md) for verified scope and limits.

The editor and visual assets are original Molecula work; no MarvinSketch source code, icons, or screenshots are included. Local analysis reports formula, supported mass values, valence warnings and graph topology (connected components, cycle rank and average degree); methods and boundaries are stated in the report. Average masses use [CIAAW standard atomic weights](https://ciaaw.org/atomic-weights.htm) (2024 table; midpoint for interval values); isotopically labelled atoms use listed isotope masses. Exact ionic masses account for electron gain/loss. Elements without a standard weight or a listed isotope mass are reported as unavailable, not guessed. Choose **Calculations → Auto-balance reaction** to solve uniquely determined positive integer coefficients while preserving element counts and formal charge. Ambiguous or unresolved systems are refused, and the edit is undoable. Third-party components and their licenses are listed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). This is a practical source and dependency review, not legal advice or a trademark clearance opinion.

## Download and first-run security warnings

Download only from the [official releases](https://github.com/jona013y/molecula/releases). Unsigned builds may show publisher-verification warnings because they are not notarized by Apple or signed with a Windows publisher certificate. That warning alone does not establish whether a file is safe; verify that it came from this repository and compare its SHA-256 if a checksum is published.

- **Mac:** Open Molecula once. If the warning appears, click **Done**. Go to **System Settings → Privacy & Security**, scroll down, click **Open Anyway** next to Molecula, then confirm **Open**. Apple describes this one-app exception in its [Gatekeeper instructions](https://support.apple.com/en-us/102445).
- **Windows:** On **Windows protected your PC**, click **More info → Run anyway**. Microsoft explains why new or unsigned downloads can show this [SmartScreen warning](https://learn.microsoft.com/windows/apps/package-and-deploy/smartscreen-reputation).

These steps apply only to the publisher-verification warning. If macOS says the app will damage your Mac, or Windows Defender reports a specific threat, cancel and do not override that alert.

## Build and run

Install Rust, then run:

```sh
cargo test --locked --offline
cargo run --release --locked --offline
```

The Windows installer script requires Inno Setup. The macOS package script creates an arm64 package on Apple Silicon; universal builds require both `aarch64-apple-darwin` and `x86_64-apple-darwin` Rust targets. macOS bundles are checked for their executable, `PkgInfo`, valid ICNS Dock icon, and CPU architecture before packaging.

The optional RDKit bridge requires a separate Python environment with RDKit installed. Its descriptor report includes RDKit logP/TPSA, H-bond counts, ring and topology measures, Labute ASA and related properties, with method and unit labels. Candidate reports enumerate up to 100 tautomers or up to 64 unique isomers at unspecified stereo centers without changing the source document. The bridge also provides SMILES conversion, 3D conformer generation, aromatize/Kekulé operations and stereochemistry reporting. RDKit-backed import is wired for MOL2, PDB and InChI; export is wired for CML, PDB, InChI and InChIKey. These conversions are unavailable until RDKit is installed. PDB export also needs 3D coordinates, and MOL2 relies on compatible atom typing. They do not provide every property or name in the full specification. See `STATUS.md` for exact coverage and limitations.

## Files and document recovery

Use **File → New/Open/Save** to work with documents. Native `.mols` preserves Molecula document data. MOL/SDF import and export supports V2000 and the implemented V3000 atom/bond subset, including charge, isotope, atom maps, supported local stereobond data and singlet/doublet/triplet radical records; RXN supports basic V2000-style reactant/product/agent components. Exchange formats can omit editor graphics, stereochemical/query detail or stoichiometry; review the loss notice before export. `.smi`/`.smiles` requires optional RDKit, as do `.mol2`, `.pdb`, and `.inchi` imports, and `.cml`, `.pdb`, `.inchi` and `.inchikey` exports. A row of tabs appears below the toolbar. Unsaved tabs are marked with a dot. Dirty workspaces autosave every 20 seconds to the application-support/config directory in a separate `workspace.autosave` file and prompt for recovery at startup. **File → Recent files** lists recently used paths.

To install RDKit separately, use a Python environment supported by RDKit and set `MOLECULA_PYTHON` to that environment's interpreter. The app remains usable without it; dependent commands show an explicit error. MOL2 parsing depends on atom typing recognized by RDKit. PDB is best treated as a coordinate-oriented exchange format, not as a lossless archive for editor metadata.

## Repository layout

- `Molecula-source.zip` is the source archive used by the installer workflow. Verify it with `shasum -a 256 -c Molecula-source.zip.sha256` on macOS/Linux or `certutil -hashfile Molecula-source.zip SHA256` on Windows.
- `Molecula-source.zip.sha256` is the published SHA-256 checksum for the archive.
- Installer workflows build from the checked-out source on the selected commit, so release builds use the exact published code.
- `THIRD_PARTY_NOTICES.md` documents dependency and bundled-font license notices.
- `src/` contains the Rust editor and chemistry model.
- `molecula.app/Contents/Resources/Molecula.icns` contains the original Molecula Dock icon.
- `scripts/verify_macos_bundle.sh` checks app bundle metadata and the icon during packaging.
- `SIGNING.md` describes optional publisher signing; no certificates or credentials belong in this repository.
