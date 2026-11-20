#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MODE="${1:-arm64}"
APP_NAME="Molecula"
STAGE="$(mktemp -d "${TMPDIR:-/tmp}/molecula-package.XXXXXX")"
trap 'rm -rf "$STAGE"' EXIT

APP="$STAGE/$APP_NAME.app"
mkdir -p "$ROOT/dist"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"
BUILD_HOME_DIR="${HOME:-}"
RUST_PATH_REMAP="--remap-path-prefix=$BUILD_HOME_DIR=~"

case "$MODE" in
  arm64)
    RUSTFLAGS="${RUSTFLAGS:+$RUSTFLAGS }$RUST_PATH_REMAP" cargo build --manifest-path "$ROOT/Cargo.toml" --release --locked --offline
    cp "$ROOT/target/release/molecula" "$APP/Contents/MacOS/MolecularSketch"
    OUTPUT="$ROOT/dist/Molecula-macOS-arm64.dmg"
    ZIP_OUTPUT="$ROOT/dist/Molecula-macOS-arm64.zip"
    ;;
  universal)
    ARM="$ROOT/target/aarch64-apple-darwin/release/molecula"
    INTEL="$ROOT/target/x86_64-apple-darwin/release/molecula"
    test -f "$ARM" && test -f "$INTEL" || {
      echo "Build both aarch64-apple-darwin and x86_64-apple-darwin release targets first." >&2
      exit 1
    }
    lipo -create "$ARM" "$INTEL" -output "$APP/Contents/MacOS/MolecularSketch"
    OUTPUT="$ROOT/dist/Molecula-macOS-universal.dmg"
    ZIP_OUTPUT="$ROOT/dist/Molecula-macOS-universal.zip"
    ;;
  *)
    echo "Usage: $0 [arm64|universal]" >&2
    exit 2
    ;;
esac

chmod 755 "$APP/Contents/MacOS/MolecularSketch"
iconutil -c icns "$ROOT/assets/Molecula.iconset" -o "$APP/Contents/Resources/Molecula.icns"
cat > "$APP/Contents/Info.plist" <<'PLIST'
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleDevelopmentRegion</key><string>en</string>
  <key>CFBundleExecutable</key><string>MolecularSketch</string>
  <key>CFBundleIconFile</key><string>Molecula.icns</string>
  <key>CFBundleIdentifier</key><string>dk.molecula.editor</string>
  <key>CFBundleInfoDictionaryVersion</key><string>6.0</string>
  <key>CFBundleName</key><string>Molecula</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleShortVersionString</key><string>1.0.13</string>
  <key>CFBundleVersion</key><string>1.0.13</string>
  <key>LSMinimumSystemVersion</key><string>11.0</string>
  <key>NSHighResolutionCapable</key><true/>
</dict>
</plist>
PLIST
printf 'APPL????' > "$APP/Contents/PkgInfo"
if [[ -n "$BUILD_HOME_DIR" ]] && strings "$APP/Contents/MacOS/MolecularSketch" | grep -Fq -- "$BUILD_HOME_DIR"; then
  echo "The packaged executable contains the build user's home path." >&2
  exit 1
fi
bash "$ROOT/scripts/verify_macos_bundle.sh" "$APP" "$MODE"
if [[ -n "${MOLECULA_SIGNING_IDENTITY:-}" ]]; then
  # Gatekeeper distribution requires a Developer ID Application certificate.
  codesign --force --deep --options runtime --timestamp \
    --sign "$MOLECULA_SIGNING_IDENTITY" "$APP"
  codesign --verify --deep --strict --verbose=2 "$APP"
else
  # Local builds remain runnable after Control-click > Open, but aren't
  # Developer ID signed and therefore cannot be notarized by Apple.
  codesign --force --deep --sign - --timestamp=none "$APP"
fi
ln -s /Applications "$STAGE/Applications"
if [[ -n "${MOLECULA_SIGNING_IDENTITY:-}" && -n "${APPLE_ID:-}" \
      && -n "${APPLE_APP_SPECIFIC_PASSWORD:-}" && -n "${APPLE_TEAM_ID:-}" ]]; then
  OPEN_NOTE="The app is signed with Developer ID and notarized by Apple."
else
  OPEN_NOTE="This build is not signed and notarized with Apple Developer ID."
fi
cat > "$STAGE/Install-Molecula.txt" <<EOF
Install Molecula

1. Drag Molecula.app onto the Applications shortcut in this window.
2. Open Molecula from Applications.
3. If macOS says the app cannot be verified, Control-click Molecula.app in Applications,
   choose Open, then confirm. $OPEN_NOTE
EOF

rm -f "$OUTPUT"
if hdiutil create -volname "$APP_NAME" -srcfolder "$STAGE" -ov -format UDZO "$OUTPUT"; then
  echo "Created $OUTPUT"
else
  rm -f "$OUTPUT"
  echo "Disk image creation is unavailable here; a drag-to-Applications archive will be created instead."
fi

# Always refresh both formats; otherwise the ZIP can silently remain from an
# older build when DMG creation succeeds.
rm -f "$ZIP_OUTPUT"
(cd "$STAGE" && /usr/bin/zip -qry -y "$ZIP_OUTPUT" .)
echo "Created $ZIP_OUTPUT"
if [[ "$MODE" == "arm64" ]]; then
  : > "$ROOT/dist/SHA256SUMS.txt"
  for ARTIFACT in "$OUTPUT" "$ZIP_OUTPUT"; do
    if [[ -f "$ARTIFACT" ]]; then
      shasum -a 256 "$ARTIFACT" | sed "s#  $ROOT/#  #" >> "$ROOT/dist/SHA256SUMS.txt"
    fi
  done
fi

if [[ -n "${MOLECULA_SIGNING_IDENTITY:-}" && -n "${APPLE_ID:-}" \
      && -n "${APPLE_APP_SPECIFIC_PASSWORD:-}" && -n "${APPLE_TEAM_ID:-}" ]]; then
  xcrun notarytool submit "$OUTPUT" --wait \
    --apple-id "$APPLE_ID" \
    --password "$APPLE_APP_SPECIFIC_PASSWORD" \
    --team-id "$APPLE_TEAM_ID"
  xcrun stapler staple "$OUTPUT"
  xcrun stapler validate "$OUTPUT"
fi
