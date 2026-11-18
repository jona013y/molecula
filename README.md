# Molecula

Molecula is an independent Rust desktop editor for chemical structures and reaction diagrams. It includes a molecular graph model, 2D editor, reaction arrows, local formula/mass analysis, MOL V2000 import/export, editable `.mols` documents, SVG export, and a searchable command palette (`Cmd/Ctrl+Shift+P`). The feature-by-feature implementation status and limitations are in [STATUS.md](STATUS.md).

The editor and visual assets are original Molecula work; no MarvinSketch source code, icons, or screenshots are included. Third-party components and their licenses are listed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). This is a practical source and dependency review, not legal advice or a trademark clearance opinion.

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

The optional RDKit descriptor bridge requires a separate Python environment with RDKit installed. See `STATUS.md` for the exact descriptor coverage and unsupported chemistry.

## Repository layout

- `Molecula-source.zip` is the source archive used by the installer workflow.
- `THIRD_PARTY_NOTICES.md` documents dependency and bundled-font license notices.
- `src/` contains the Rust editor and chemistry model.
- `molecula.app/Contents/Resources/Molecula.icns` contains the original Molecula Dock icon.
- `scripts/verify_macos_bundle.sh` checks app bundle metadata and the icon during packaging.
- `SIGNING.md` describes optional publisher signing; no certificates or credentials belong in this repository.
