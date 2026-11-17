# Molecula

Molecula is an independent Rust desktop editor for chemical structures and reaction diagrams. It includes a molecular graph model, 2D editor, reaction arrows, local formula/mass analysis, MOL V2000 import/export, editable `.mols` documents, and SVG export. The feature-by-feature implementation status and limitations are in [STATUS.md](STATUS.md).

## Build and run

Install Rust, then run:

```sh
cargo test --locked --offline
cargo run --release --locked --offline
```

The Windows installer script requires Inno Setup. The macOS package script creates an arm64 package on Apple Silicon; universal builds require both `aarch64-apple-darwin` and `x86_64-apple-darwin` Rust targets. macOS bundles are checked for their executable, `PkgInfo`, valid ICNS Dock icon, and CPU architecture before packaging.

The optional RDKit descriptor bridge requires a separate Python environment with RDKit installed. See `STATUS.md` for the exact descriptor coverage and unsupported chemistry.

## Repository layout

- `Molecula-source.zip` is the source archive used by the installer workflow.
- `src/` contains the Rust editor and chemistry model.
- `molecula.app/Contents/Resources/Molecula.icns` contains the original Molecula Dock icon.
- `scripts/verify_macos_bundle.sh` checks app bundle metadata and the icon during packaging.
- `SIGNING.md` describes optional publisher signing; no certificates or credentials belong in this repository.
