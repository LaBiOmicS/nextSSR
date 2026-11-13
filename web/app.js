// nextSSR web interface (main thread).
import { detectLang, setLang, t, lang } from "./i18n.js";

const $ = (id) => document.getElementById(id);
const worker = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });

const DEFAULT_REPEATS = { 1: 10, 2: 6, 3: 5, 4: 5, 5: 5, 6: 5 };
let input = null; // { file, name } | { text, name, labelKey }
let last = null; // { summary, preview, pairs, seconds }
let lastEpcr = null; // { summary, seconds } of the standalone e-PCR
let lastProgress = null;
let tool = "ssr"; // "ssr" | "epcr"

// ---------- language ----------
function applyLang(l) {
  setLang(l);
  document.querySelectorAll(".lang button").forEach((b) => b.classList.toggle("on", b.dataset.lang === lang()));
  renderUnitLabels();
  setTool(tool);
  if (input) $("file-name").textContent = inputLabel();
  if (tool === "ssr" && last) render();
  if (tool === "epcr" && lastEpcr) renderEpcr();
  if (lastProgress && !last) showProgress(lastProgress);
}
document.querySelectorAll(".lang button").forEach((b) => b.addEventListener("click", () => applyLang(b.dataset.lang)));

// ---------- parameters ----------
const units = $("units");
for (let u = 1; u <= 6; u++) {
  const l = document.createElement("label");
  l.innerHTML = `<input type="number" min="0" max="999" id="u${u}" value="${DEFAULT_REPEATS[u]}"> <span id="ul${u}"></span>`;
  units.appendChild(l);
}
function renderUnitLabels() {
  for (let u = 1; u <= 6; u++) {
    $(`ul${u}`).textContent = t(`u.short.${u}`);
    $(`ul${u}`).parentElement.title = t(`u.${u}`);
  }
}
function syncEpcrOpts() {
  $("epcr-opts").hidden = !(tool === "epcr" || ($("primers").checked && $("epcr").checked));
}
$("primers").addEventListener("change", () => {
  $("primer-opts").hidden = !$("primers").checked;
  syncEpcrOpts();
});
$("epcr").addEventListener("change", syncEpcrOpts);

// ---------- tool ----------
function setTool(v) {
  tool = v;
  document.querySelectorAll(".seg button").forEach((b) => b.classList.toggle("on", b.dataset.tool === tool));
  $("ssr-params").hidden = tool !== "ssr";
  $("epcr-params").hidden = tool !== "epcr";
  syncEpcrOpts();
  $("run-btn").textContent = t(tool === "ssr" ? "app.run" : "app.run.epcr");
  $("empty").textContent = t(tool === "ssr" ? "app.empty" : "app.empty.epcr");
  const has = tool === "ssr" ? last : lastEpcr;
  $("empty").hidden = !!has;
  $("results").hidden = tool !== "ssr" || !last;
  $("epcr-results").hidden = tool !== "epcr" || !lastEpcr;
}
document.querySelectorAll(".seg button").forEach((b) => b.addEventListener("click", () => setTool(b.dataset.tool)));
$("mode").addEventListener("change", () => {
  const misa = $("mode").value === "misa";
  $("imperfect").disabled = misa;
  if (misa) $("imperfect").checked = false;
});

function num(id) {
  const v = parseFloat($(id).value);
  return Number.isFinite(v) ? v : undefined;
}

function options() {
  const def = [];
  for (let u = 1; u <= 6; u++) {
    const r = parseInt($(`u${u}`).value, 10);
    if (r > 0) def.push(`${u}-${r}`);
  }
  return {
    definition: def.join(" "),
    mode: $("mode").value,
    compound: $("compound").checked,
    max_interruption: num("maxint"),
    imperfect: $("imperfect").checked,
    max_gap: num("maxgap"),
    min_purity: num("purity"),
    primers: $("primers").checked,
    product_min: num("pmin"), product_max: num("pmax"),
    tm_min: num("tmin"), tm_opt: num("topt"), tm_max: num("tmax"),
    size_min: num("smin"), size_opt: num("sopt"), size_max: num("smax"),
    gc_min: num("gmin"), gc_max: num("gmax"),
    num_pairs: num("npairs"),
    epcr: $("primers").checked && $("epcr").checked,
    ...epcrOptions(),
    flank: num("flank") || 0,
    preview_limit: 5000,
  };
}

function epcrOptions() {
  return {
    mismatches: num("mm"),
    three_prime_mismatches: num("mm3"),
    three_prime_len: num("len3"),
    epcr_max_product: num("emax"),
  };
}

// ---------- input ----------
function inputLabel() {
  if (!input) return "";
  if (input.file) return `${input.file.name} (${fmtBytes(input.file.size)})`;
  return t(input.labelKey);
}
function setInput(v) {
  input = v;
  $("file-name").textContent = inputLabel();
  $("run-btn").disabled = !input;
}
$("file").addEventListener("change", (e) => {
  const f = e.target.files[0];
  if (f) setInput({ file: f, name: f.name });
});
const drop = $("drop");
["dragenter", "dragover"].forEach((ev) => drop.addEventListener(ev, (e) => { e.preventDefault(); drop.classList.add("over"); }));
["dragleave", "drop"].forEach((ev) => drop.addEventListener(ev, (e) => { e.preventDefault(); drop.classList.remove("over"); }));
drop.addEventListener("drop", (e) => {
  const f = e.dataTransfer.files[0];
  if (f) setInput({ file: f, name: f.name });
});
$("paste").addEventListener("input", () => {
  const txt = $("paste").value;
  if (txt.trim()) setInput({ text: txt.trim().startsWith(">") ? txt : `>pasted\n${txt}`, name: "pasted.fasta", labelKey: "dyn.pasted" });
  else setInput(null);
});
$("example").addEventListener("click", async () => {
  const r = await fetch("examples/NC_000932.1_Arabidopsis_cp.fasta");
  setInput({ text: await r.text(), name: "NC_000932.1_Arabidopsis_cp.fasta", labelKey: "dyn.example" });
});

// ---------- run ----------
function showError(msg) {
  $("error").textContent = t("dyn.error", { m: msg });
  $("error").hidden = false;
}

$("run-btn").addEventListener("click", () => {
  if (!input) return;
  $("error").hidden = true;
  const src = { file: input.file || null, text: input.text || null, name: input.name };
  if (tool === "epcr") {
    const primers = $("primer-list").value;
    if (!primers.trim()) return showError(t("dyn.epcr.nopairs"));
    lastEpcr = null;
    worker.postMessage({ type: "epcr", ...src, primers, options: epcrOptions() });
  } else {
    last = null;
    worker.postMessage({ type: "run", ...src, options: options() });
  }
  $("progress").hidden = false;
  $("bar").style.width = "0%";
  $("progress-text").textContent = t("dyn.starting");
  $("run-btn").disabled = true;
  lastProgress = null;
});

function showProgress(m) {
  const pct = m.total ? (100 * m.read) / m.total : 0;
  $("bar").style.width = `${pct.toFixed(1)}%`;
  const vars = { read: fmtBytes(m.read), total: fmtBytes(m.total), seqs: fmt(m.stats.sequences), ssrs: fmt(m.stats.ssrs) };
  $("progress-text").textContent = t(m.phase === "epcr" ? "dyn.epcr.progress" : "dyn.progress", vars);
}

worker.onmessage = (e) => {
  const m = e.data;
  if (m.type === "version") $("version").textContent = `v${m.version}`;
  else if (m.type === "progress") { lastProgress = m; showProgress(m); }
  else if (m.type === "done") {
    $("bar").style.width = "100%";
    $("run-btn").disabled = false;
    last = { summary: m.summary, preview: m.preview, pairs: m.pairs, seconds: m.seconds };
    if (tool === "ssr") render();
  } else if (m.type === "epcr-done") {
    $("bar").style.width = "100%";
    $("run-btn").disabled = false;
    lastEpcr = { summary: m.summary, seconds: m.seconds };
    if (tool === "epcr") renderEpcr();
  } else if (m.type === "output") {
    download(m.text, m.format, m.tool);
  } else if (m.type === "error") {
    showError(m.message);
    $("run-btn").disabled = false;
  }
};

// ---------- rendering ----------
function fmtBytes(n) {
  if (n < 1024) return `${n} B`;
  const u = ["KB", "MB", "GB", "TB"];
  let i = -1;
  do { n /= 1024; i++; } while (n >= 1024 && i < u.length - 1);
  return `${n.toLocaleString(lang(), { maximumFractionDigits: n < 10 ? 1 : 0 })} ${u[i]}`;
}
function fmtBp(n) {
  if (n < 1e3) return `${n} bp`;
  const [d, s] = n >= 1e9 ? [1e9, "Gb"] : n >= 1e6 ? [1e6, "Mb"] : [1e3, "kb"];
  return `${(n / d).toLocaleString(lang(), { maximumFractionDigits: 2 })} ${s}`;
}
const fmt = (n) => Number(n).toLocaleString(lang());
const fmtD = (n, d) => Number(n).toLocaleString(lang(), { minimumFractionDigits: d, maximumFractionDigits: d });
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]);

function render() {
  const { summary: sum, seconds } = last;
  const s = sum.stats;
  $("progress-text").textContent = t("dyn.done", { s: fmtD(seconds, 2) });
  $("results").hidden = false;
  $("epcr-results").hidden = true;
  $("empty").hidden = true;
  const mb = s.total_bp / 1e6;
  const cards = [
    ["card.seqs", fmt(s.sequences)],
    ["card.size", fmtBp(s.total_bp)],
    ["card.ssrs", fmt(s.ssrs)],
    ["card.permb", mb > 0 ? fmtD(s.ssrs / mb, 1) : "–"],
    ["card.loci", fmt(s.loci)],
    ["card.compound", fmt(s.compounds)],
    ["card.imperfect", fmt(s.imperfect)],
  ];
  if (sum.config.primers) cards.push(["card.primers", `${fmt(s.primers_designed)} / ${fmt(s.loci)}`]);
  if (sum.epcr) cards.push(["card.locispec", `${fmt(sum.epcr.loci_with_specific_pair)} / ${fmt(s.primers_designed)}`]);
  cards.push(["card.time", `${fmtD(seconds, 2)} s`]);
  $("cards").innerHTML = cards.map(([k, v]) => `<div class="card"><b>${v}</b><span>${t(k)}</span></div>`).join("");

  const fmts = ["tsv", "gff3", "bed", "misa", "statistics", "json", "p3in"];
  if (sum.config.primers) fmts.push("primers");
  if (sum.epcr) fmts.push("epcr");
  downloadButtons($("downloads"), fmts, "ssr");

  const perUnit = Object.entries(s.per_unit).map(([u, c]) => [Number(u), c]);
  const maxU = Math.max(1, ...perUnit.map(([, c]) => c));
  $("unit-bars").innerHTML = perUnit
    .map(([u, c]) => `<div class="b"><span>${t(`u.${u}`)}</span><span class="tr"><i style="width:${(100 * c) / maxU}%"></i></span><span class="n">${fmt(c)}</span></div>`)
    .join("") || `<p>${t("dyn.none")}</p>`;

  const classes = Object.entries(sum.per_class)
    .map(([k, reps]) => [k, Object.values(reps).reduce((a, b) => a + b, 0)])
    .sort((a, b) => b[1] - a[1])
    .slice(0, 12);
  $("classes").innerHTML = `<tr><th>${t("th.class")}</th><th class="num">${t("th.ssrs")}</th><th class="num">%</th></tr>` +
    classes.map(([k, c]) => `<tr><td class="mono">${k}</td><td class="num">${fmt(c)}</td><td class="num">${fmtD((100 * c) / s.ssrs, 1)}</td></tr>`).join("");

  $("table-note").textContent = s.loci > last.preview.length ? t("dyn.first", { n: fmt(last.preview.length), total: fmt(s.loci) }) : `(${fmt(s.loci)})`;
  renderTable();
}

function renderTable() {
  if (!last) return;
  const q = $("filter").value.trim().toLowerCase();
  const withPrimers = last.summary.config.primers;
  const rows = last.preview.filter((r) => !q || `${r.seq_id} ${r.ssr_type} ${r.pattern}`.toLowerCase().includes(q)).slice(0, 2000);
  let h = `<tr><th>${t("th.seq")}</th><th class="num">#</th><th>${t("th.type")}</th><th>${t("th.ssr")}</th><th class="num">${t("th.size")}</th><th class="num">${t("th.start")}</th><th class="num">${t("th.end")}</th>`;
  const pairs = last.pairs;
  if (withPrimers) h += `<th>${t("th.fwd")}</th><th>${t("th.rev")}</th><th class="num">${t("th.product")}</th>`;
  if (withPrimers && pairs) h += `<th>${t("th.spec")}</th>`;
  h += `</tr>`;
  for (const r of rows) {
    h += `<tr><td>${esc(r.seq_id)}</td><td class="num">${r.nr}</td><td><span class="pill">${esc(r.ssr_type)}</span></td><td class="pat" title="${esc(r.pattern)}">${esc(r.pattern)}</td><td class="num">${fmt(r.end - r.start + 1)}</td><td class="num">${fmt(r.start)}</td><td class="num">${fmt(r.end)}</td>`;
    if (withPrimers) {
      const o = r.primers || {};
      const p = o.status === "designed" ? o.pairs[0] : null;
      h += p
        ? `<td class="mono" title="Tm ${p.left.tm} °C">${p.left.seq}</td><td class="mono" title="Tm ${p.right.tm} °C">${p.right.seq}</td><td class="num">${p.product_size}</td>`
        : `<td colspan="3" class="muted">${esc(o.reason || "–")}</td>`;
      if (pairs) {
        const st = p && pairs[`${r.seq_id}_${r.nr}_1`];
        h += st ? `<td>${statusPill(st[0], st[1])}</td>` : `<td class="muted">–</td>`;
      }
    }
    h += `</tr>`;
  }
  $("loci").innerHTML = h;
}
$("filter").addEventListener("input", renderTable);

function statusPill(status, amplicons) {
  return `<span class="pill st-${status}" title="${t("th.amplicons")}: ${fmt(amplicons)}">${t(`st.${status}`)}</span>`;
}

function downloadButtons(el, fmts, from) {
  el.innerHTML = fmts.map((f) => `<button type="button" data-f="${f}">⬇ ${t(`dl.${f}`)}</button>`).join("");
  el.querySelectorAll("button").forEach((b) => b.addEventListener("click", () => worker.postMessage({ type: "output", format: b.dataset.f, tool: from })));
}

// ---------- standalone e-PCR ----------
function renderEpcr() {
  const { summary: sum, seconds } = lastEpcr;
  const c = sum.pairs_by_status;
  $("progress-text").textContent = t("dyn.done", { s: fmtD(seconds, 2) });
  $("empty").hidden = true;
  $("results").hidden = true;
  $("epcr-results").hidden = false;
  const amplicons = sum.pairs.reduce((a, p) => a + p.amplicons, 0);
  const cards = [
    ["card.seqs", fmt(sum.sequences)],
    ["card.size", fmtBp(sum.bp)],
    ["card.primers", fmt(sum.pairs.length)],
    ["card.specific", fmt(c.specific || 0)],
    ["card.nonspecific", fmt(c.nonspecific || 0)],
    ["card.noproduct", fmt(c.no_product || 0)],
    ["card.amplicons", fmt(amplicons)],
    ["card.time", `${fmtD(seconds, 2)} s`],
  ];
  $("ecards").innerHTML = cards.map(([k, v]) => `<div class="card"><b>${v}</b><span>${t(k)}</span></div>`).join("");
  $("esat").hidden = !sum.pairs.some((p) => p.saturated);
  downloadButtons($("edownloads"), ["epcr", "pairs"], "epcr");
  renderPairs();
}

function renderPairs() {
  if (!lastEpcr) return;
  const q = $("efilter").value.trim().toLowerCase();
  const rows = lastEpcr.summary.pairs.filter((p) => !q || `${p.id} ${p.forward} ${p.reverse} ${p.status}`.toLowerCase().includes(q)).slice(0, 2000);
  let h = `<tr><th>${t("th.pair")}</th><th>${t("th.fwd")}</th><th>${t("th.rev")}</th><th>${t("th.status")}</th><th class="num">${t("th.amplicons")}</th></tr>`;
  for (const p of rows) {
    h += `<tr><td>${esc(p.id)}</td><td class="mono">${esc(p.forward)}</td><td class="mono">${esc(p.reverse)}</td><td>${statusPill(p.status, p.amplicons)}</td><td class="num">${fmt(p.amplicons)}</td></tr>`;
  }
  $("epairs").innerHTML = h;
}
$("efilter").addEventListener("input", renderPairs);

function download(text, format, from) {
  const ext = { json: "summary.json", primers: "primers.tsv", epcr: "epcr.tsv", pairs: "epcr_pairs.tsv" }[format] || format;
  const base = (input && input.name ? input.name : "nextssr").replace(/\.(gz|bgz)$/i, "");
  const url = URL.createObjectURL(new Blob([text], { type: "text/plain" }));
  const a = Object.assign(document.createElement("a"), { href: url, download: `${base}.${ext}` });
  document.body.appendChild(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 5000);
}

applyLang(detectLang());
worker.postMessage({ type: "version" });
