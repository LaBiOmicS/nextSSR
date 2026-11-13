// Web Worker: runs the nextSSR WebAssembly engine off the main thread.
import init, { Analyzer, version } from "./pkg/nextssr_wasm.js";

const ready = init();
let analyzer = null;

self.onmessage = async (e) => {
  const msg = e.data;
  try {
    await ready;
    if (msg.type === "version") {
      self.postMessage({ type: "version", version: version() });
    } else if (msg.type === "run") {
      await run(msg);
    } else if (msg.type === "output") {
      const text = analyzer ? analyzer.output(msg.format) : "";
      self.postMessage({ type: "output", format: msg.format, text });
    }
  } catch (err) {
    self.postMessage({ type: "error", message: String(err && err.message ? err.message : err) });
  }
};

async function run({ file, text, name, options }) {
  if (analyzer) analyzer.free();
  analyzer = new Analyzer(JSON.stringify(options), name);
  const t0 = performance.now();
  const blob = file || new Blob([text]);
  const total = blob.size;
  const reader = blob.stream().getReader();
  let read = 0;
  let lastPost = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    // feed in slices so progress stays responsive on large chunks
    for (let off = 0; off < value.length; off += 1 << 20) {
      analyzer.push(value.subarray(off, Math.min(off + (1 << 20), value.length)));
    }
    read += value.length;
    const now = performance.now();
    if (now - lastPost > 120) {
      lastPost = now;
      self.postMessage({ type: "progress", read, total, stats: JSON.parse(analyzer.progress()) });
    }
  }
  const summary = JSON.parse(analyzer.finish());
  const preview = JSON.parse(analyzer.preview());
  self.postMessage({
    type: "done",
    summary,
    preview,
    seconds: (performance.now() - t0) / 1000,
  });
}
