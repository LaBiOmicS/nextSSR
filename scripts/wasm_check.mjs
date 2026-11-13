// Usage: node scripts/wasm_check.mjs <fasta> <options-json> <out-prefix>
// Runs the WebAssembly build in Node and writes its outputs, for comparison with the CLI.
import { readFileSync, writeFileSync, createReadStream } from "node:fs";
import { fileURLToPath } from "node:url";
import path from "node:path";
import { initSync, Analyzer } from "../web/pkg/nextssr_wasm.js";

const here = path.dirname(fileURLToPath(import.meta.url));
initSync({ module: readFileSync(path.join(here, "../web/pkg/nextssr_wasm_bg.wasm")) });
const [fasta, opts, prefix] = process.argv.slice(2);
const t0 = performance.now();
const a = new Analyzer(opts || "{}", fasta);
for await (const chunk of createReadStream(fasta, { highWaterMark: 1 << 20 })) a.push(chunk);
const summary = JSON.parse(a.finish());
for (const f of ["tsv", "gff3", "misa", "statistics", "bed"]) writeFileSync(`${prefix}.${f}`, a.output(f));
console.log(`wasm: ${summary.stats.sequences} seqs, ${summary.stats.ssrs} SSRs, ${(performance.now() - t0).toFixed(0)} ms`);
