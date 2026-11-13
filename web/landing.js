// nextSSR presentation page: hero scanner, live MISA × nextSSR demo, pipeline, Tm calculator.
import { detectLang, setLang, t, lang } from "./i18n.js";
import init, { Analyzer, melting_temperature, version } from "./pkg/nextssr_wasm.js";

const $ = (id) => document.getElementById(id);
const reduced = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" })[c]);

// ---------- header ----------
const top = $("top");
const onScroll = () => top.classList.toggle("solid", window.scrollY > 40);
window.addEventListener("scroll", onScroll, { passive: true });
onScroll();

// ---------- hero scanner ----------
function rng(seed) {
  let x = seed >>> 0;
  return () => {
    x ^= x << 13; x >>>= 0; x ^= x >>> 17; x ^= x << 5; x >>>= 0;
    return x / 4294967296;
  };
}
const PLANT = [
  ["AG", 9], ["A", 14], ["CTT", 6], ["GATA", 5], ["AC", 8], ["T", 12], ["TTTA", 5],
  ["AAG", 7], ["TC", 10], ["CAG", 6], ["G", 11], ["ATT", 6], ["GT", 9], ["AATG", 5],
];
function buildScanner() {
  const rows = $("rows");
  const r = rng(20261113);
  const base = () => "ACGT"[Math.floor(r() * 4)];
  const nRows = 5;
  let k = 0;
  let html = "";
  for (let i = 0; i < nRows; i++) {
    let row = "";
    let len = 0;
    const width = 190;
    let next = 6 + Math.floor(r() * 30);
    while (len < width) {
      if (len >= next) {
        const [m, n] = PLANT[k++ % PLANT.length];
        const s = m.repeat(n);
        row += `<span class="ssr c${k % 3}" data-l="(${m})${n}">${s}<span class="tag">(${m})${n}</span></span>`;
        len += s.length;
        next = len + 14 + Math.floor(r() * 46);
      } else {
        row += base();
        len++;
      }
    }
    html += `<div class="row" style="transform:translateX(${-Math.floor(r() * 40)}ch)">${row}</div>`;
  }
  rows.innerHTML = html;
}
function runScanner() {
  const box = $("scanner");
  const beam = $("beam");
  const spans = [...box.querySelectorAll(".ssr")];
  const hud = $("hud-n");
  if (reduced) {
    spans.forEach((s) => s.classList.add("lit"));
    hud.textContent = spans.length;
    beam.style.display = "none";
    return;
  }
  const period = 7000;
  let t0 = performance.now();
  let visible = true;
  new IntersectionObserver(([e]) => (visible = e.isIntersecting)).observe(box);
  function frame(now) {
    requestAnimationFrame(frame);
    if (!visible) return;
    const w = box.clientWidth;
    let p = (now - t0) / period;
    if (p >= 1) {
      t0 = now;
      p = 0;
      spans.forEach((s) => s.classList.remove("lit"));
    }
    const x = p * (w + 80) - 40;
    beam.style.transform = `translateX(${x}px)`;
    const bx = box.getBoundingClientRect().left + x;
    let n = 0;
    for (const s of spans) {
      const r = s.getBoundingClientRect();
      if (r.left < bx && r.right > 0) {
        s.classList.add("lit");
      }
      if (s.classList.contains("lit") && r.left > -1 && r.right < box.getBoundingClientRect().right) n++;
    }
    hud.textContent = n;
  }
  requestAnimationFrame(frame);
}

// ---------- counters & reveal ----------
function countUp(el) {
  const target = Number(el.dataset.count);
  if (reduced || target === 0) { el.textContent = target; return; }
  const t0 = performance.now();
  const d = 1400;
  const step = (now) => {
    const p = Math.min(1, (now - t0) / d);
    el.textContent = Math.round(target * (1 - Math.pow(1 - p, 3)));
    if (p < 1) requestAnimationFrame(step);
  };
  requestAnimationFrame(step);
}
const io = new IntersectionObserver((entries) => {
  for (const e of entries) {
    if (!e.isIntersecting) continue;
    e.target.classList.add("in");
    if (e.target.dataset.count !== undefined) countUp(e.target);
    if (e.target.id === "bars") e.target.querySelectorAll(".bar-fill").forEach((b) => (b.style.width = `${Math.max(0.8, b.dataset.w)}%`));
    io.unobserve(e.target);
  }
}, { threshold: 0.2 });
document.querySelectorAll(".reveal, [data-count], #bars").forEach((el) => io.observe(el));

// ---------- pipeline ----------
const steps = [...document.querySelectorAll(".step")];
const stages = [...document.querySelectorAll(".stage")];
function showStep(i) {
  steps.forEach((s, j) => s.classList.toggle("on", i === j));
  stages.forEach((s, j) => (s.hidden = i !== j));
}
steps.forEach((s, i) => s.addEventListener("click", () => showStep(i)));

// ---------- live demo ----------
const CASES = {
  masking: "GATTCCGCCCCGAAAAAAAAAAAAAAAGAAGAAGAAGAAGCCTCGTACCATG",
  overlap: "GCATCGTTCCAAAAAAAAAAAAAGAGAGAGAGAGAGAGCGTTACGGATCCTAGC",
  imperfect: "TGCATCCAGAGAGAGAGAGAGTAGAGAGAGCCATTGCAGTC",
  // Arabidopsis thaliana chloroplast (NC_000932.1), 119058–119176
  distant: "TTTCTCTTTTTTTTTTTCTAAAGTGTATTTTCAACCTATCAATTTCAAGAAAGTAATTACGAATATATTATATTAAAAGAATGTGCCTTAATACTTAATATTATTATTATTATATAAAT",
};
let currentCase = "masking";
let wasmReady = false;

function analyze(seq, mode) {
  const a = new Analyzer(JSON.stringify({ mode, imperfect: mode === "native", preview_members: true }), "demo");
  a.push(new TextEncoder().encode(`>demo\n${seq}\n`));
  a.finish();
  const loci = JSON.parse(a.preview());
  a.free();
  return loci;
}

// Positions are percentages of the lane, so every track lines up with the
// bases whatever the font size of the labels inside it.
const pct = (x, n) => `${(100 * x) / n}%`;

function lane(loci, cls, n) {
  if (!loci.length) return `<div class="lane feat-lane"><span class="feat-none">${t("l.demo.none")}</span></div>`;
  let h = `<div class="lane feat-lane">`;
  for (const l of loci) {
    const len = l.end - l.start + 1;
    const mem = (l.members || []).length > 1
      ? l.members.map((m) => `<span class="mem" style="left:${pct(m.start - l.start, len)};width:${pct(m.end - m.start + 1, len)}"></span>`).join("")
      : "";
    h += `<div class="feat ${cls}" style="left:${pct(l.start - 1, n)};width:${pct(len, n)}" title="${esc(l.ssr_type)} ${esc(l.pattern)} · ${l.start}–${l.end}">${mem}<span class="lbl">${esc(l.pattern)}</span></div>`;
  }
  return h + "</div>";
}

function callList(loci) {
  if (!loci.length) return `<li>${t("l.demo.none")}</li>`;
  return loci.map((l) => `<li><span class="ty">${esc(l.ssr_type)}</span><span>${esc(l.pattern)}</span><span class="co">${l.start}–${l.end}</span></li>`).join("");
}

function renderDemo() {
  const own = currentCase === "own";
  $("own").classList.toggle("show", own);
  $("case-h").textContent = t(`l.case.${currentCase}.t`);
  $("case-p").innerHTML = t(`l.case.${currentCase}.p`);
  const v = $("verdict");
  if (!wasmReady) return;
  let seq = own ? $("own-seq").value : CASES[currentCase];
  seq = seq.toUpperCase().replace(/[^ACGTN]/g, "").slice(0, 400);
  const m = analyze(seq, "misa");
  const x = analyze(seq, "native");
  const n = seq.length;
  let ruler = "";
  for (let i = 10; i <= n; i += 10) ruler += `<i style="left:${pct(i - 0.5, n)}"></i><span style="left:${pct(i - 0.5, n)}">${i}</span>`;
  let bases = "";
  for (let i = 0; i < n; i++) bases += `<b class="${seq[i]}">${seq[i]}</b>`;
  const viewer = $("viewer");
  viewer.innerHTML = `<div class="track-grid">
    <div class="track-name"></div><div class="lane ruler">${ruler}</div>
    <div class="track-name">${t("l.demo.seq")}</div><div class="lane seq">${bases}</div>
    <div class="track-name misa">MISA</div>${lane(m, "misa", n)}
    <div class="track-name nx">nextSSR</div>${lane(x, "nx", n)}
  </div>`;
  fitViewer(n);
  $("call-misa").innerHTML = callList(m);
  $("call-nx").innerHTML = callList(x);
  const same = JSON.stringify(m.map((l) => [l.start, l.end, l.pattern])) === JSON.stringify(x.map((l) => [l.start, l.end, l.pattern]));
  v.className = `verdict ${same ? "ok" : "bad"}`;
  v.textContent = t(same ? "l.demo.same" : "l.demo.diff");
}
document.querySelectorAll("#cases button").forEach((b) =>
  b.addEventListener("click", () => {
    currentCase = b.dataset.case;
    document.querySelectorAll("#cases button").forEach((x) => x.classList.toggle("on", x === b));
    renderDemo();
  }),
);
let ownTimer = 0;
$("own-seq").addEventListener("input", () => {
  clearTimeout(ownTimer);
  ownTimer = setTimeout(renderDemo, 150);
});

// Scale the sequence so it fills the demo width (scrolls only when even the
// smallest readable size does not fit).
function fitViewer(n) {
  const viewer = $("viewer");
  if (!n) return;
  const avail = viewer.clientWidth - viewer.querySelector(".track-name").offsetWidth - 4;
  let fs = Math.max(11, Math.min(30, avail / (n * 0.6)));
  viewer.style.setProperty("--fs", `${fs.toFixed(2)}px`);
  // correct the estimate with the real width of the rendered bases
  const seq = viewer.querySelector(".seq");
  const w = seq ? seq.scrollWidth : 0;
  if (w > avail && fs > 11) {
    fs = Math.max(11, (fs * avail) / w);
    viewer.style.setProperty("--fs", `${fs.toFixed(2)}px`);
  }
}
window.addEventListener("resize", () => {
  const n = $("viewer").querySelectorAll(".seq b").length;
  fitViewer(n);
});

// ---------- Tm calculator ----------
function renderTm() {
  if (!wasmReady) return;
  const s = $("tm-in").value.toUpperCase().replace(/[^ACGT]/g, "");
  const tm = s.length >= 8 ? melting_temperature(s) : NaN;
  const gc = s.length ? (100 * (s.match(/[GC]/g) || []).length) / s.length : NaN;
  const f = (x, d) => (Number.isFinite(x) ? x.toLocaleString(lang(), { maximumFractionDigits: d, minimumFractionDigits: d }) : "–");
  $("tm-v").textContent = f(tm, 1);
  $("tm-gc").textContent = Number.isFinite(gc) ? `${f(gc, 0)}%` : "–";
  $("tm-len").textContent = s.length || "–";
  const pos = Number.isFinite(tm) ? Math.min(100, Math.max(0, ((tm - 45) / 30) * 100)) : 0;
  $("tm-g").style.left = `calc(${pos}% - 2px)`;
}
$("tm-in").addEventListener("input", renderTm);

// ---------- language ----------
function applyLang(l) {
  setLang(l);
  document.querySelectorAll(".lang button").forEach((b) => b.classList.toggle("on", b.dataset.lang === lang()));
  renderDemo();
  renderTm();
}
document.querySelectorAll(".lang button").forEach((b) => b.addEventListener("click", () => applyLang(b.dataset.lang)));

buildScanner();
runScanner();
applyLang(detectLang());
init().then(() => {
  wasmReady = true;
  $("version").textContent = `v${version()}`;
  $("cite-ver").textContent = version();
  renderDemo();
  renderTm();
});
