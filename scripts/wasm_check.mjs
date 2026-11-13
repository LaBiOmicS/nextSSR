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
const formats = ["tsv", "gff3", "misa", "statistics", "bed"];
const o = JSON.parse(opts || "{}");
if (o.primers) {
  formats.push("p3in");
  if (o.epcr) {
    // second pass: specificity of the designed pairs
    const e = a.epcr_scanner();
    for await (const chunk of createReadStream(fasta, { highWaterMark: 1 << 20 })) e.push(chunk);
    e.finish();
    a.set_epcr(e);
    writeFileSync(`${prefix}.epcr.tsv`, a.output("epcr"));
  }
  writeFileSync(`${prefix}.primers.tsv`, a.output("primers"));
}
for (const f of formats) writeFileSync(`${prefix}.${f}`, a.output(f));
console.log(`wasm: ${summary.stats.sequences} seqs, ${summary.stats.ssrs} SSRs, ${(performance.now() - t0).toFixed(0)} ms`);
