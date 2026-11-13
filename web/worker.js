// Web Worker: runs the nextSSR WebAssembly engine off the main thread.
import init, { Analyzer, EpcrScanner, version } from "./pkg/nextssr_wasm.js";

const ready = init();
let analyzer = null;
let scanner = null; // standalone e-PCR

self.onmessage = async (e) => {
  const msg = e.data;
  try {
    await ready;
    if (msg.type === "version") {
      self.postMessage({ type: "version", version: version() });
    } else if (msg.type === "run") {
      await run(msg);
    } else if (msg.type === "epcr") {
      await runEpcr(msg);
    } else if (msg.type === "output") {
      const src = msg.tool === "epcr" ? scanner : analyzer;
      const text = src ? src.output(msg.format) : "";
      self.postMessage({ type: "output", format: msg.format, tool: msg.tool, text });
    }
  } catch (err) {
    self.postMessage({ type: "error", message: String(err && err.message ? err.message : err) });
  }
};

/** Stream a Blob through `sink.push`, posting progress messages. */
async function stream(blob, sink, phase) {
  const total = blob.size;
  const reader = blob.stream().getReader();
  let read = 0;
  let lastPost = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    // feed in slices so progress stays responsive on large chunks
    for (let off = 0; off < value.length; off += 1 << 20) {
      sink.push(value.subarray(off, Math.min(off + (1 << 20), value.length)));
    }
    read += value.length;
    const now = performance.now();
    if (now - lastPost > 120) {
      lastPost = now;
      self.postMessage({ type: "progress", phase, read, total, stats: JSON.parse(sink.progress()) });
    }
  }
}

function pairStatuses(summary) {
  const m = {};
  for (const p of summary.pairs) m[p.id] = [p.status, p.amplicons];
  return m;
}

async function run({ file, text, name, options }) {
  if (analyzer) analyzer.free();
  analyzer = new Analyzer(JSON.stringify(options), name);
  const t0 = performance.now();
  const blob = file || new Blob([text]);
  await stream(blob, analyzer, "ssr");
  analyzer.finish();
  let pairs = null;
  if (options.primers && options.epcr) {
    // second pass over the same input: specificity of every designed pair
    const e = analyzer.epcr_scanner();
    try {
      await stream(blob, e, "epcr");
      pairs = pairStatuses(JSON.parse(e.finish()));
      analyzer.set_epcr(e);
    } finally {
      e.free();
    }
  }
  self.postMessage({
    type: "done",
    summary: JSON.parse(analyzer.summary()),
    preview: JSON.parse(analyzer.preview()),
    pairs,
    seconds: (performance.now() - t0) / 1000,
  });
}

async function runEpcr({ file, text, primers, options }) {
  if (scanner) scanner.free();
  scanner = null;
  const t0 = performance.now();
  const s = EpcrScanner.from_list(primers, JSON.stringify(options));
  scanner = s;
  await stream(file || new Blob([text]), s, "epcr");
  const summary = JSON.parse(s.finish());
  self.postMessage({ type: "epcr-done", summary, seconds: (performance.now() - t0) / 1000 });
}
