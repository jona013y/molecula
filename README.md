# Molecula

Molecula is an independent Rust desktop editor for chemical structures and reaction diagrams. It includes a molecular graph model, 2D editor, reaction arrows, local formula/mass analysis, MOL V2000 import/export, editable `.mols` documents, and SVG export. The feature-by-feature implementation status and limitations are in [STATUS.md](STATUS.md).

## Download and first-run security warnings

Download Molecula only from the [official v1.0.2 GitHub release](https://github.com/jona013y/molecula/releases/tag/v1.0.2). The Mac and Windows files are not signed by commercial publisher certificates. On macOS the app is ad-hoc signed and not notarized by Apple; on Windows the installer has no Authenticode signature. This is why Gatekeeper may say that Apple cannot check the app, and SmartScreen may say Windows protected your PC or identify an unknown publisher. Those warnings describe missing publisher verification or file reputation; they do not prove that a file is safe or malicious.

Before opening either installer, verify its SHA-256 against the value shown here and on the release page:

| Download | SHA-256 |
| --- | --- |
| `Molecula-macOS-universal.dmg` | `99134a51f4464b11ecef41c7021a84f23c4b5eaeb6b1828f60ae62e7219c7510` |
| `Molecula-Windows-Setup.exe` | `50ab40b9a0c11b4fcc5f06c53d73962a6d7f021a7ea44c5e9974d911e18e9f58` |

**macOS:** In Terminal, run `shasum -a 256 ~/Downloads/Molecula-macOS-universal.dmg` and compare the result with the Mac hash above. Open the disk image and drag Molecula to Applications. If Gatekeeper says the developer cannot be verified or Apple cannot check the app for malicious software, try opening it once, then go to **System Settings → Privacy & Security → Open Anyway** and confirm. You can also Control-click Molecula in Applications, choose **Open**, then confirm **Open**. This adds an exception for Molecula; it does not turn Gatekeeper off globally. These steps are for the unidentified-developer/notarization warning only. If macOS says the app will damage your Mac, do not override that alert; delete the download and check the release page and hash again. See [Apple’s instructions for safely opening apps](https://support.apple.com/en-us/102445).

**Windows:** In PowerShell, run `(Get-FileHash "$env:USERPROFILE\Downloads\Molecula-Windows-Setup.exe" -Algorithm SHA256).Hash` and compare the result with the Windows hash above. If the verified installer shows the blue **Windows protected your PC** SmartScreen warning for an unrecognized app, select **More info → Run anyway**. Do not turn SmartScreen off. If Microsoft Defender reports a specific threat, cancel the installation instead. See [Microsoft’s explanation of SmartScreen reputation](https://learn.microsoft.com/windows/apps/package-and-deploy/smartscreen-reputation).

These warnings can be avoided more reliably by using Apple Developer ID signing and notarization for macOS, or Microsoft Store distribution for Windows. Those routes involve paid developer services or extra distribution requirements. We are keeping this personal-use release unsigned by commercial publishers, so the one-time warning may remain. Apple explains its notarization and Gatekeeper checks in its [safe app opening guidance](https://support.apple.com/en-us/102445); Microsoft explains that new files can warn until reputation builds, even when signed, in its [SmartScreen guide](https://learn.microsoft.com/windows/apps/package-and-deploy/smartscreen-reputation).

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
