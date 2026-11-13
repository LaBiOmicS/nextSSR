// nextSSR web interface (main thread).
import { detectLang, setLang, t, lang } from "./i18n.js";

const $ = (id) => document.getElementById(id);
const worker = new Worker(new URL("./worker.js", import.meta.url), { type: "module" });

const DEFAULT_REPEATS = { 1: 10, 2: 6, 3: 5, 4: 5, 5: 5, 6: 5 };
let input = null; // { file, name } | { text, name, labelKey }
let last = null; // { summary, preview, seconds }
let lastProgress = null;

// ---------- language ----------
function applyLang(l) {
  setLang(l);
  document.querySelectorAll(".lang button").forEach((b) => b.classList.toggle("on", b.dataset.lang === lang()));
  renderUnitLabels();
  if (input) $("file-name").textContent = inputLabel();
  if (last) render();
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
$("primers").addEventListener("change", () => ($("primer-opts").hidden = !$("primers").checked));
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
    flank: num("flank") || 0,
    preview_limit: 5000,
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
$("run-btn").addEventListener("click", () => {
  if (!input) return;
  $("error").hidden = true;
  $("progress").hidden = false;
  $("bar").style.width = "0%";
  $("progress-text").textContent = t("dyn.starting");
  $("run-btn").disabled = true;
  last = null;
  lastProgress = null;
  worker.postMessage({ type: "run", file: input.file || null, text: input.text || null, name: input.name, options: options() });
});

function showProgress(m) {
  const pct = m.total ? (100 * m.read) / m.total : 0;
  $("bar").style.width = `${pct.toFixed(1)}%`;
  $("progress-text").textContent = t("dyn.progress", {
    read: fmtBytes(m.read), total: fmtBytes(m.total), seqs: fmt(m.stats.sequences), ssrs: fmt(m.stats.ssrs),
  });
}

worker.onmessage = (e) => {
  const m = e.data;
  if (m.type === "version") $("version").textContent = `v${m.version}`;
  else if (m.type === "progress") { lastProgress = m; showProgress(m); }
  else if (m.type === "done") {
    $("bar").style.width = "100%";
    $("run-btn").disabled = false;
    last = { summary: m.summary, preview: m.preview, seconds: m.seconds };
    render();
  } else if (m.type === "output") {
    download(m.text, m.format);
  } else if (m.type === "error") {
    $("error").textContent = t("dyn.error", { m: m.message });
    $("error").hidden = false;
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
  cards.push(["card.time", `${fmtD(seconds, 2)} s`]);
  $("cards").innerHTML = cards.map(([k, v]) => `<div class="card"><b>${v}</b><span>${t(k)}</span></div>`).join("");

  const fmts = ["tsv", "gff3", "bed", "misa", "statistics", "json", "p3in"];
  $("downloads").innerHTML = fmts.map((f) => `<button type="button" data-f="${f}">⬇ ${t(`dl.${f}`)}</button>`).join("");
  $("downloads").querySelectorAll("button").forEach((b) => b.addEventListener("click", () => worker.postMessage({ type: "output", format: b.dataset.f })));

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
  if (withPrimers) h += `<th>${t("th.fwd")}</th><th>${t("th.rev")}</th><th class="num">${t("th.product")}</th>`;
  h += `</tr>`;
  for (const r of rows) {
    h += `<tr><td>${esc(r.seq_id)}</td><td class="num">${r.nr}</td><td><span class="pill">${esc(r.ssr_type)}</span></td><td class="pat" title="${esc(r.pattern)}">${esc(r.pattern)}</td><td class="num">${fmt(r.end - r.start + 1)}</td><td class="num">${fmt(r.start)}</td><td class="num">${fmt(r.end)}</td>`;
    if (withPrimers) {
      const p = r.primers || {};
      h += p.status === "designed"
        ? `<td class="mono" title="Tm ${p.left.tm} °C">${p.left.seq}</td><td class="mono" title="Tm ${p.right.tm} °C">${p.right.seq}</td><td class="num">${p.product_size}</td>`
        : `<td colspan="3" class="muted">${esc(p.reason || "–")}</td>`;
    }
    h += `</tr>`;
  }
  $("loci").innerHTML = h;
}
$("filter").addEventListener("input", renderTable);

function download(text, format) {
  const ext = { json: "summary.json" }[format] || format;
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
