#!/usr/bin/env bash
# Build the WebAssembly engine into web/pkg (requires wasm32-unknown-unknown and wasm-bindgen-cli 0.2.100).
set -euo pipefail
cd "$(dirname "$0")/.."
cargo build --release --target wasm32-unknown-unknown -p nextssr-wasm
wasm-bindgen --target web --no-typescript --out-dir web/pkg target/wasm32-unknown-unknown/release/nextssr_wasm.wasm
if command -v wasm-opt >/dev/null; then wasm-opt -Os -o web/pkg/nextssr_wasm_bg.wasm web/pkg/nextssr_wasm_bg.wasm; fi
ls -la web/pkg
