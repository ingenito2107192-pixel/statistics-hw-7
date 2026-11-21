/* ==========================================================
   script.js — Homework 7 · Statistics
   Shell (theme + navigation), Exercise 1 and Exercise 2.
   Vanilla JavaScript, no external libraries.
   ========================================================== */

/* ==========================================================
   SHELL — theme toggle and section highlighting in the nav
   ========================================================== */
(() => {
  'use strict';
  const root = document.documentElement;
  const themeBtn = document.getElementById('theme-toggle');
  function setTheme(t, notify) {
    root.dataset.theme = t;
    themeBtn.textContent = t === 'dark' ? '☀ Light' : '☾ Dark';
    try { localStorage.setItem('theme', t); } catch (e) { /* storage may be unavailable */ }
    if (notify) document.dispatchEvent(new Event('themechange'));
  }
  let saved = null;
  try { saved = localStorage.getItem('theme'); } catch (e) { /* ignore */ }
  setTheme(saved || (window.matchMedia && matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light'), false);
  themeBtn.addEventListener('click', () => setTheme(root.dataset.theme === 'dark' ? 'light' : 'dark', true));

  const links = Array.from(document.querySelectorAll('#nav a'));
  const map = new Map(links.map(a => [a.getAttribute('href').slice(1), a]));
  if ('IntersectionObserver' in window) {
    const io = new IntersectionObserver(entries => {
      entries.forEach(e => {
        if (e.isIntersecting) {
          links.forEach(l => l.classList.remove('active'));
          const a = map.get(e.target.id);
          if (a) a.classList.add('active');
        }
      });
    }, { rootMargin: '-20% 0px -70% 0px' });
    map.forEach((_, id) => { const s = document.getElementById(id); if (s) io.observe(s); });
  }
})();

/* ==========================================================
   EXERCISE 1 — Bayes' theorem & modes of convergence
   ========================================================== */
(() => {
'use strict';
/* ---------- pure maths helpers (testable without a browser) ---------- */
function lgamma(x) { // Lanczos approximation of log Gamma
  const c = [76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let y = x, t = x + 5.5;
  t -= (x + 0.5) * Math.log(t);
  let s = 1.000000000190015;
  for (let j = 0; j < 6; j++) { y += 1; s += c[j] / y; }
  return -t + Math.log(2.5066282746310005 * s / x);
}
function betaPdf(x, a, b) {
  if (x <= 0 || x >= 1) return 0;
  return Math.exp((a - 1) * Math.log(x) + (b - 1) * Math.log(1 - x) -
    (lgamma(a) + lgamma(b) - lgamma(a + b)));
}
function erf(x) { // Abramowitz-Stegun 7.1.26
  const s = Math.sign(x); x = Math.abs(x);
  const t = 1 / (1 + 0.3275911 * x);
  const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return s * y;
}
const normPdf = x => Math.exp(-x * x / 2) / Math.sqrt(2 * Math.PI);
const normCdf = x => 0.5 * (1 + erf(x / Math.SQRT2));
const cauchyPdf = x => 1 / (Math.PI * (1 + x * x));
const cauchyCdf = x => 0.5 + Math.atan(x) / Math.PI;

function bayesIDS(p, s, f) { // prior, sensitivity, false-alarm rate
  const ev = s * p + f * (1 - p);
  const post = s * p / ev;
  const postNoAlert = (1 - s) * p / (1 - ev);
  const priorOdds = p / (1 - p), LR = s / f;
  return { ev, post, postNoAlert, priorOdds, LR, postOdds: priorOdds * LR };
}

/* ----- convergence lab: sequences Y_n -> 0 (n = 1..N), stored in x[1..N] ----- */
const SEQ = {
  mean: {
    gen(N, x) { let s = 0; for (let n = 1; n <= N; n++) { s += Math.random() - 0.5; x[n] = s / n; } },
    msq: n => 1 / (12 * n),
    yr: [-0.5, 0.5],
    desc: 'Y_n = (1/n) Σ (U_i − 1/2). By the LLN this is the sample mean minus μ = 1/2, with σ² = 1/12.',
    ok: [true, true, true, true],
    why: [
      'Chebyshev: P(|Y_n| > ε) ≤ σ²/(nε²) → 0 (weak LLN).',
      'Strong LLN (Kolmogorov): almost every trajectory converges to 0.',
      'E[Y_n²] = σ²/n = 1/(12n) → 0.',
      'Implied by convergence in probability (limit is a constant).'
    ]
  },
  type: {
    gen(N, x) {
      const u = Math.random();
      for (let n = 1; n <= N; n++) {
        const k = Math.floor(Math.log2(n)), w = Math.pow(2, -k), lo = (n - Math.pow(2, k)) * w;
        x[n] = (u >= lo && u < lo + w) ? 1 : 0;
      }
    },
    msq: n => Math.pow(2, -Math.floor(Math.log2(n))),
    yr: [-0.1, 1.1],
    desc: 'Write n = 2^k + j with 0 ≤ j < 2^k, and let Y_n = 1 if U ∈ [j/2^k, (j+1)/2^k), else 0. The interval gets shorter, but sweeps [0,1) again and again.',
    ok: [true, false, true, true],
    why: [
      'P(Y_n = 1) = 2^(−k) → 0.',
      'For every ω, Y_n(ω) = 1 for infinitely many n, so Y_n(ω) does not converge.',
      'E[Y_n²] = P(Y_n = 1) = 2^(−k) → 0.',
      'Implied by convergence in probability.'
    ]
  },
  spike: {
    gen(N, x) { const u = Math.random(); for (let n = 1; n <= N; n++) x[n] = u < 1 / n ? n : 0; },
    msq: n => n,
    yr: [-1, 20],
    desc: 'Y_n = n if U < 1/n, else 0. Tall but rarer and rarer spikes (values above the chart are clipped).',
    ok: [true, true, false, true],
    why: [
      'P(Y_n > ε) = 1/n → 0.',
      'For each ω with U > 0: Y_n(ω) = 0 for all n ≥ 1/U, so the trajectory converges.',
      'E[Y_n²] = n² · (1/n) = n → ∞: the rare spikes are too big.',
      'Implied by convergence in probability.'
    ]
  }
};

function simulate(seq, M, N, eps, keep) {
  const x = new Float64Array(N + 1);
  const cnt = new Float64Array(N + 1), sq = new Float64Array(N + 1), sup = new Float64Array(N + 1);
  const paths = [];
  for (let m = 0; m < M; m++) {
    seq.gen(N, x);
    if (m < keep) paths.push(Array.from(x));
    let mx = 0;
    for (let n = N; n >= 1; n--) {
      const a = Math.abs(x[n]);
      if (a > eps) cnt[n]++;
      sq[n] += a * a;
      if (a > mx) mx = a;
      if (mx > eps) sup[n]++;
    }
  }
  for (let n = 1; n <= N; n++) { cnt[n] /= M; sq[n] /= M; sup[n] /= M; }
  return { paths, prob: cnt, msq: sq, sup, N, eps };
}

/* ----- CLT lab ----- */
const DIST = {
  unif:   { s: () => Math.random(), mu: 0.5, sg: Math.sqrt(1 / 12), pdf: normPdf, cdf: normCdf, std: true },
  exp:    { s: () => -Math.log(1 - Math.random()), mu: 1, sg: 1, pdf: normPdf, cdf: normCdf, std: true },
  bern:   { s: () => (Math.random() < 0.1 ? 1 : 0), mu: 0.1, sg: 0.3, pdf: normPdf, cdf: normCdf, std: true },
  cauchy: { s: () => Math.tan(Math.PI * (Math.random() - 0.5)), mu: 0, sg: 1, pdf: cauchyPdf, cdf: cauchyCdf, std: false }
};
const CLT_N = [1, 2, 3, 5, 10, 30, 100, 500];

function cltSample(d, n, K) {
  const D = DIST[d], z = new Float64Array(K);
  for (let k = 0; k < K; k++) {
    let s = 0;
    for (let i = 0; i < n; i++) s += D.s();
    z[k] = D.std ? (s - n * D.mu) / (D.sg * Math.sqrt(n)) : s / n; // Cauchy: plain mean
  }
  return z;
}
function ksDistance(z, cdf) {
  const a = Float64Array.from(z).sort(), K = a.length;
  let d = 0;
  for (let i = 0; i < K; i++) {
    const F = cdf(a[i]);
    d = Math.max(d, Math.abs((i + 1) / K - F), Math.abs(F - i / K));
  }
  return d;
}

if (typeof module !== 'undefined') {
  module.exports = { lgamma, betaPdf, erf, normCdf, bayesIDS, SEQ, simulate, cltSample, ksDistance, DIST };
}

/* ====================== browser part ====================== */
if (typeof document !== 'undefined') {
  const $ = id => document.getElementById('e1-' + id); // all Exercise-1 ids carry the e1- prefix
  const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  const pct = (x, d = 2) => (100 * x).toFixed(d) + ' %';

  /* ---- tiny plotting toolkit on <canvas> ---- */
  function ticks(a, b, n) {
    const raw = (b - a) / n, p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p;
    const step = (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * p, out = [];
    for (let v = Math.ceil(a / step - 1e-9) * step; v <= b + 1e-9; v += step) out.push(+v.toFixed(10));
    return out;
  }
  const nf = v => (v !== 0 && (Math.abs(v) >= 10000 || Math.abs(v) < 0.001)) ? v.toExponential(0) : String(+v.toFixed(3));

  function frame(canvas, o) {
    const dpr = window.devicePixelRatio || 1, w = canvas.clientWidth, h = canvas.clientHeight;
    canvas.width = Math.round(w * dpr); canvas.height = Math.round(h * dpr);
    const g = canvas.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0);
    const m = { l: 50, r: 14, t: 12, b: 40 }, pw = w - m.l - m.r, ph = h - m.t - m.b;
    const X = x => m.l + (x - o.xmin) / (o.xmax - o.xmin) * pw;
    const Y = y => m.t + ph - (y - o.ymin) / (o.ymax - o.ymin) * ph;
    g.font = '11px system-ui, sans-serif'; g.lineWidth = 1;
    g.strokeStyle = css('--grid'); g.fillStyle = css('--muted');
    g.textAlign = 'center'; g.textBaseline = 'top';
    ticks(o.xmin, o.xmax, 6).forEach(v => {
      g.beginPath(); g.moveTo(X(v), m.t); g.lineTo(X(v), m.t + ph); g.stroke();
      g.fillText(nf(v), X(v), m.t + ph + 4);
    });
    g.textAlign = 'right'; g.textBaseline = 'middle';
    ticks(o.ymin, o.ymax, 5).forEach(v => {
      g.beginPath(); g.moveTo(m.l, Y(v)); g.lineTo(m.l + pw, Y(v)); g.stroke();
      g.fillText(nf(v), m.l - 6, Y(v));
    });
    g.textAlign = 'center'; g.textBaseline = 'bottom';
    if (o.xlabel) g.fillText(o.xlabel, m.l + pw / 2, h - 2);
    if (o.ylabel) { g.save(); g.translate(12, m.t + ph / 2); g.rotate(-Math.PI / 2); g.textBaseline = 'top'; g.fillText(o.ylabel, 0, -4); g.restore(); }
    const clip = fn => { g.save(); g.beginPath(); g.rect(m.l, m.t, pw, ph); g.clip(); fn(); g.restore(); };
    return {
      g, X, Y, m, pw, ph,
      line(xs, ys, color, width = 2, dash = []) {
        clip(() => {
          g.strokeStyle = color; g.lineWidth = width; g.setLineDash(dash); g.beginPath();
          for (let i = 0; i < xs.length; i++) i ? g.lineTo(X(xs[i]), Y(ys[i])) : g.moveTo(X(xs[i]), Y(ys[i]));
          g.stroke();
        });
      },
      vline(x, color, dash = [5, 4]) { this.line([x, x], [o.ymin, o.ymax], color, 2, dash); },
      band(y1, y2, color) { clip(() => { g.fillStyle = color; g.fillRect(m.l, Y(y2), pw, Y(y1) - Y(y2)); }); },
      bars(edges, vals, color) {
        clip(() => {
          g.fillStyle = color;
          for (let i = 0; i < vals.length; i++) g.fillRect(X(edges[i]) + .5, Y(vals[i]), X(edges[i + 1]) - X(edges[i]) - 1, Y(0) - Y(vals[i]));
        });
      }
    };
  }

  /* ================= Example 1: IDS ================= */
  function renderIDS() {
    const p = +$('ids-prior').value / 100, s = +$('ids-sens').value / 100, f = +$('ids-fpr').value / 100;
    $('ids-prior-v').textContent = pct(p, 1);
    $('ids-sens-v').textContent = pct(s, 0);
    $('ids-fpr-v').textContent = pct(f, 1);
    const r = bayesIDS(p, s, f);
    $('ids-evidence').textContent = pct(r.ev);
    $('ids-post').textContent = pct(r.post, 1);
    $('ids-npv').textContent = pct(r.postNoAlert, 3);
    $('ids-odds').textContent = r.priorOdds.toFixed(3) + ' × ' + r.LR.toFixed(1) + ' = ' + r.postOdds.toFixed(3);
    $('bar-prior').style.width = (100 * p) + '%';
    $('bar-post').style.width = (100 * r.post) + '%';

    const N = 1000, A = Math.round(N * p), TP = Math.round(A * s), FN = A - TP;
    const FP = Math.round((N - A) * f), TN = N - A - FP;
    $('n-tp').textContent = TP; $('n-fn').textContent = FN; $('n-fp').textContent = FP; $('n-tn').textContent = TN;
    $('ids-msg').innerHTML = TP + FP > 0
      ? 'Out of <b>' + (TP + FP) + '</b> alerts, only <b>' + TP + '</b> are real attacks: ' +
        TP + ' / ' + (TP + FP) + ' ≈ <b>' + pct(TP / (TP + FP), 1) + '</b> (the exact posterior is ' + pct(r.post, 1) + ').'
      : 'No alerts in this sample of 1000 events.';

    const c = $('ids-dots'), dpr = window.devicePixelRatio || 1, w = c.clientWidth, h = c.clientHeight;
    c.width = Math.round(w * dpr); c.height = Math.round(h * dpr);
    const g = c.getContext('2d'); g.setTransform(dpr, 0, 0, dpr, 0, 0);
    const cols = 50, rows = 20, cw = w / cols, ch = h / rows, rad = Math.max(1.5, Math.min(cw, ch) / 2 - 1.2);
    const colors = [[TP, '--c-tp'], [FN, '--c-fn'], [FP, '--c-fp'], [TN, '--c-tn']];
    let i = 0;
    colors.forEach(([cnt, v]) => {
      g.fillStyle = css(v);
      for (let k = 0; k < cnt; k++, i++) {
        const cx = (i % cols + 0.5) * cw, cy = (Math.floor(i / cols) + 0.5) * ch;
        g.beginPath(); g.arc(cx, cy, rad, 0, 2 * Math.PI); g.fill();
      }
    });
  }
  ['ids-prior', 'ids-sens', 'ids-fpr'].forEach(id => $(id).addEventListener('input', renderIDS));

  /* ================= Example 2: Beta-Binomial ================= */
  const bb = { h: 0, t: 0 };
  const GRID = 400, xsB = Array.from({ length: GRID }, (_, i) => (i + 0.5) / GRID);

  function renderBB() {
    const [a, b] = $('bb-prior').value.split(',').map(Number);
    $('bb-true-v').textContent = (+$('bb-true').value).toFixed(2);
    const A = a + bb.h, B = b + bb.t, n = bb.h + bb.t;
    const pr = xsB.map(x => betaPdf(x, a, b));
    const lk = xsB.map(x => betaPdf(x, bb.h + 1, bb.t + 1));
    const po = xsB.map(x => betaPdf(x, A, B));
    const ymax = Math.max(...pr, ...lk, ...po) * 1.1;
    const f = frame($('bb-canvas'), { xmin: 0, xmax: 1, ymin: 0, ymax, xlabel: 'θ = P(heads)', ylabel: 'density' });
    f.line(xsB, pr, css('--muted'), 2, [6, 4]);
    f.line(xsB, lk, css('--accent2'), 2, [2, 3]);
    f.line(xsB, po, css('--accent'), 3);
    f.vline(+$('bb-true').value, css('--ok'));

    let cum = 0, lo = 0, hi = 1, gotLo = false;
    for (let i = 0; i < GRID; i++) {
      cum += po[i] / GRID;
      if (!gotLo && cum >= 0.025) { lo = xsB[i]; gotLo = true; }
      if (cum >= 0.975) { hi = xsB[i]; break; }
    }
    $('bb-data').textContent = bb.h + ' H / ' + bb.t + ' T';
    $('bb-post').textContent = 'Beta(' + A + ', ' + B + ')';
    $('bb-mean').textContent = (A / (A + B)).toFixed(3) + (n ? '  (MLE ' + (bb.h / n).toFixed(3) + ')' : '');
    $('bb-ci').textContent = '[' + lo.toFixed(3) + ', ' + hi.toFixed(3) + ']';
  }
  function flips(k) {
    const p = +$('bb-true').value;
    for (let i = 0; i < k; i++) Math.random() < p ? bb.h++ : bb.t++;
    renderBB();
  }
  $('bb-1').onclick = () => flips(1);
  $('bb-10').onclick = () => flips(10);
  $('bb-100').onclick = () => flips(100);
  $('bb-1000').onclick = () => flips(1000);
  $('bb-reset').onclick = () => { bb.h = bb.t = 0; renderBB(); };
  $('bb-prior').addEventListener('change', renderBB);
  $('bb-true').addEventListener('input', renderBB);

  /* ================= Convergence lab ================= */
  let labSim = null;
  const PATH_COLORS = ['#2f5bea', '#e8590c', '#1a8f4c', '#c2255c', '#7048e8', '#0c8599'];

  function labRun() {
    const key = $('lab-seq').value, N = +$('lab-N').value, eps = +$('lab-eps').value;
    $('lab-eps-v').textContent = eps.toFixed(2);
    $('lab-N-v').textContent = N;
    const M = 2000;
    labSim = simulate(SEQ[key], M, N, eps, 6);
    labSim.key = key;
    $('lab-info').textContent = M + ' simulated trajectories, N = ' + N;
    $('lab-desc').textContent = SEQ[key].desc;
    labRender();
    renderVerdicts(key);
  }
  function labRender() {
    if (!labSim) return;
    const { paths, prob, msq, sup, N, eps, key } = labSim, S = SEQ[key];
    const xs = Array.from({ length: N }, (_, i) => i + 1);

    let f = frame($('lab-paths'), { xmin: 1, xmax: N, ymin: S.yr[0], ymax: S.yr[1], xlabel: 'n', ylabel: 'Y_n(ω)' });
    f.band(-eps, eps, css('--grid'));
    paths.forEach((p, i) => f.line(xs, xs.map(n => p[n]), PATH_COLORS[i % 6], 1.3));

    f = frame($('lab-prob'), { xmin: 1, xmax: N, ymin: 0, ymax: 1, xlabel: 'n', ylabel: 'probability' });
    f.line(xs, xs.map(n => sup[n]), css('--accent2'), 2.5);
    f.line(xs, xs.map(n => prob[n]), css('--accent'), 2.5);

    const est = xs.map(n => msq[n]), ex = xs.map(n => S.msq(n));
    const ymax = Math.max(...est, ...ex) * 1.1;
    f = frame($('lab-l2'), { xmin: 1, xmax: N, ymin: 0, ymax, xlabel: 'n', ylabel: 'E[Y_n²]' });
    f.line(xs, ex, css('--muted'), 2, [6, 4]);
    f.line(xs, est, css('--accent'), 2);
  }
  function renderVerdicts(key) {
    const S = SEQ[key], names = ['In probability', 'Almost sure', 'Mean square (L²)', 'In distribution'];
    $('lab-verdicts').innerHTML = names.map((nm, i) =>
      '<div class="verdict ' + (S.ok[i] ? 'ok' : 'no') + '"><h4>' + nm + ' <span class="mark">' +
      (S.ok[i] ? '✔ converges' : '✘ does not converge') + '</span></h4>' + S.why[i] + '</div>').join('');
  }
  $('lab-run').onclick = labRun;
  $('lab-seq').addEventListener('change', labRun);
  $('lab-eps').addEventListener('change', labRun);
  $('lab-N').addEventListener('change', labRun);
  $('lab-eps').addEventListener('input', () => $('lab-eps-v').textContent = (+$('lab-eps').value).toFixed(2));
  $('lab-N').addEventListener('input', () => $('lab-N-v').textContent = $('lab-N').value);

  /* ================= CLT lab ================= */
  let cltZ = null;
  function cltRun() {
    const d = $('clt-dist').value, n = CLT_N[+$('clt-n').value];
    $('clt-n-v').textContent = 'n = ' + n;
    cltZ = { z: cltSample(d, n, 4000), d, n };
    cltRender();
  }
  function cltRender() {
    if (!cltZ) return;
    const { z, d, n } = cltZ, D = DIST[d], K = z.length;
    const lo = -4, hi = 4, nb = 40, bw = (hi - lo) / nb, cnt = new Array(nb).fill(0);
    let sum = 0, sum2 = 0;
    z.forEach(v => { sum += v; sum2 += v * v; const b = Math.floor((v - lo) / bw); if (b >= 0 && b < nb) cnt[b]++; });
    const dens = cnt.map(c => c / (K * bw)), edges = cnt.map((_, i) => lo + i * bw).concat([hi]);
    const grid = Array.from({ length: 201 }, (_, i) => lo + (hi - lo) * i / 200);
    const pdf = grid.map(D.pdf);
    const ymax = Math.max(...dens, ...pdf) * 1.15;
    const f = frame($('clt-canvas'), { xmin: lo, xmax: hi, ymin: 0, ymax, xlabel: D.std ? 'z  (standardized mean Z_n)' : 'x  (sample mean, NOT standardized)', ylabel: 'density' });
    f.bars(edges, dens, css('--accent') + '99');
    f.line(grid, pdf, css('--accent2'), 3);

    const mean = sum / K, v = sum2 / K - mean * mean;
    $('clt-ks').textContent = ksDistance(z, D.cdf).toFixed(3);
    $('clt-mean').textContent = D.std ? mean.toFixed(3) : 'n/a (no mean)';
    $('clt-var').textContent = D.std ? v.toFixed(3) : 'n/a (infinite)';
    $('clt-msg').textContent = D.std
      ? 'n = ' + n + ': the histogram approaches the bell curve as n grows (the Kolmogorov distance shrinks like 1/√n). Skewed or discrete laws need larger n.'
      : 'Cauchy: X̄_n has exactly the same Cauchy law for every n. There is no LLN (the mean does not exist) and no CLT: the distance does not decrease.';
  }
  $('clt-run').onclick = cltRun;
  $('clt-dist').addEventListener('change', cltRun);
  $('clt-n').addEventListener('input', cltRun);

  /* ================= start ================= */
  document.addEventListener('themechange', () => { renderIDS(); renderBB(); labRender(); cltRender(); });
  let tmr = null;
  window.addEventListener('resize', () => {
    clearTimeout(tmr);
    tmr = setTimeout(() => { renderIDS(); renderBB(); labRender(); cltRender(); }, 120);
  });
  renderIDS(); renderBB(); labRun(); cltRun();
}
})();

/* ==========================================================
   EXERCISE 2 — convergence, LLN & CLT, sqrt(n), Bayes
   ========================================================== */
/* ==========================================================
   script.js — interactive demos for the Statistics homework
   Modes of convergence, LLN & CLT, sqrt(n), Bayes' theorem.
   No external libraries.
   ========================================================== */
(() => {
  'use strict';

  /* ---------- tiny helpers ---------- */
  const $ = (s, r = document) => r.querySelector(s);
  const $$ = (s, r = document) => Array.from(r.querySelectorAll(s));
  const cssv = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  const clamp = (x, a, b) => Math.min(b, Math.max(a, x));
  const range = (a, b, s) => { const r = []; for (let v = a; v <= b + 1e-9; v += s) r.push(+v.toFixed(10)); return r; };
  const redraws = [];
  const redrawAll = () => redraws.forEach(f => f());

  function niceStep(x) {
    const p = Math.pow(10, Math.floor(Math.log10(x)));
    const f = x / p;
    return (f < 1.5 ? 1 : f < 3 ? 2 : f < 7 ? 5 : 10) * p;
  }

  document.addEventListener('themechange', redrawAll);

  let resizeTimer;
  window.addEventListener('resize', () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(redrawAll, 120); });

  /* ---------- canvas helpers ---------- */
  function prep(c) {
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const w = c.clientWidth, h = c.clientHeight;
    const W = Math.round(w * dpr), H = Math.round(h * dpr);
    if (c.width !== W || c.height !== H) { c.width = W; c.height = H; }
    const ctx = c.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h);
    ctx.font = '12px ' + cssv('--sans');
    return { ctx, w, h };
  }

  function axes(ctx, w, h, o) {
    const m = { l: 56, r: 14, t: 12, b: 42 };
    const iw = w - m.l - m.r, ih = h - m.t - m.b;
    const fx = o.xlog ? Math.log10 : v => v;
    const fy = o.ylog ? Math.log10 : v => v;
    const x0 = fx(o.xmin), x1 = fx(o.xmax), y0 = fy(o.ymin), y1 = fy(o.ymax);
    const X = v => m.l + (fx(v) - x0) / (x1 - x0) * iw;
    const Y = v => m.t + ih - (fy(v) - y0) / (y1 - y0) * ih;
    const line = cssv('--line'), txt = cssv('--muted');
    ctx.lineWidth = 1;
    ctx.strokeStyle = line; ctx.fillStyle = txt;
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    (o.xticks || []).forEach(t => {
      const x = X(t);
      ctx.beginPath(); ctx.moveTo(x, m.t); ctx.lineTo(x, m.t + ih); ctx.stroke();
      ctx.fillText(o.xfmt ? o.xfmt(t) : String(t), x, m.t + ih + 6);
    });
    ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    (o.yticks || []).forEach(t => {
      const y = Y(t);
      ctx.beginPath(); ctx.moveTo(m.l, y); ctx.lineTo(m.l + iw, y); ctx.stroke();
      ctx.fillText(o.yfmt ? o.yfmt(t) : String(t), m.l - 6, y);
    });
    ctx.strokeStyle = txt; ctx.strokeRect(m.l, m.t, iw, ih);
    ctx.fillStyle = txt; ctx.textAlign = 'center'; ctx.textBaseline = 'bottom';
    if (o.xlabel) ctx.fillText(o.xlabel, m.l + iw / 2, h - 2);
    if (o.ylabel) {
      ctx.save(); ctx.translate(13, m.t + ih / 2); ctx.rotate(-Math.PI / 2);
      ctx.textBaseline = 'top'; ctx.fillText(o.ylabel, 0, -4); ctx.restore();
    }
    return { X, Y, m, iw, ih };
  }
  const clipTo = (ctx, a) => { ctx.save(); ctx.beginPath(); ctx.rect(a.m.l, a.m.t, a.iw, a.ih); ctx.clip(); };

  /* ---------- random numbers & distributions ---------- */
  function randn() {
    let u = 0;
    while (u === 0) u = Math.random();
    return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * Math.random());
  }
  const DISTS = {
    bernoulli: { mu: 0.3, sd: Math.sqrt(0.21), draw: () => (Math.random() < 0.3 ? 1 : 0) },
    dice: { mu: 3.5, sd: Math.sqrt(35 / 12), draw: () => 1 + Math.floor(Math.random() * 6) },
    uniform: { mu: 0.5, sd: Math.sqrt(1 / 12), draw: () => Math.random() },
    exp: { mu: 1, sd: 1, draw: () => -Math.log(1 - Math.random()) },
    cauchy: { mu: null, sd: null, draw: () => Math.tan(Math.PI * (Math.random() - 0.5)) }
  };
  function erf(x) {
    const s = Math.sign(x); x = Math.abs(x);
    const t = 1 / (1 + 0.3275911 * x);
    const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
    return s * y;
  }
  const Phi = z => 0.5 * (1 + erf(z / Math.SQRT2));
  const phi = z => Math.exp(-z * z / 2) / Math.sqrt(2 * Math.PI);

  function lgamma(z) {
    const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
      -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
    if (z < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * z)) - lgamma(1 - z);
    z -= 1;
    let x = c[0];
    for (let i = 1; i < 9; i++) x += c[i] / (z + i);
    const t = z + 7.5;
    return 0.5 * Math.log(2 * Math.PI) + (z + 0.5) * Math.log(t) - t + Math.log(x);
  }
  const betaPdf = (x, a, b) =>
    Math.exp((a - 1) * Math.log(x) + (b - 1) * Math.log(1 - x) + lgamma(a + b) - lgamma(a) - lgamma(b));

  const fmt = (v, d = 3) => (Number.isFinite(v) ? v.toFixed(d) : '–');
  const pct = (v, d = 1) => (v * 100).toFixed(d) + '%';

  /* ==========================================================
     DEMO 1 — Law of large numbers
     ========================================================== */
  const lln = { paths: null, N: 0, eps: 0.1, center: 0, key: 'bernoulli', y: [0, 1] };

  function llnCompute() {
    const key = $('#lln-dist').value;
    const D = DISTS[key];
    const N = +$('#lln-n').value;
    const eps = +$('#lln-eps').value;
    const center = D.mu === null ? 0 : D.mu;
    const K = 8, M = 400;

    const paths = [];
    for (let k = 0; k < K; k++) {
      const arr = new Float32Array(N);
      let s = 0;
      for (let i = 0; i < N; i++) { s += D.draw(); arr[i] = s / (i + 1); }
      paths.push(arr);
    }

    // y-range
    let lo = Infinity, hi = -Infinity;
    paths.forEach(p => { for (let i = Math.min(4, N - 1); i < N; i++) { if (p[i] < lo) lo = p[i]; if (p[i] > hi) hi = p[i]; } });
    lo = Math.min(lo, center - eps * 1.6); hi = Math.max(hi, center + eps * 1.6);
    const cap = key === 'cauchy' ? 6 : 4 * D.sd;
    lo = Math.max(lo, center - cap); hi = Math.min(hi, center + cap);
    const pad = (hi - lo) * 0.06; lo -= pad; hi += pad;

    // P(|mean_n - mu| > eps) from M independent paths at checkpoints
    const cps = [...new Set([10, 100, 1000, N].filter(c => c <= N))].sort((a, b) => a - b);
    const counts = new Array(cps.length).fill(0);
    for (let m = 0; m < M; m++) {
      let s = 0, ci = 0;
      for (let i = 1; i <= N && ci < cps.length; i++) {
        s += D.draw();
        if (i === cps[ci]) { if (Math.abs(s / i - center) > eps) counts[ci]++; ci++; }
      }
    }
    $('#lln-table').innerHTML = cps.map((c, i) =>
      `<tr><td>${c.toLocaleString('en-US')}</td><td>${(counts[i] / M * 100).toFixed(1)}%</td></tr>`).join('');

    // last exit times
    const ul = $('#lln-exits');
    ul.innerHTML = '';
    const pal = ['--c1', '--c2', '--c3', '--c4', '--c5', '--c6', '--c7', '--c8'].map(cssv);
    paths.forEach((p, k) => {
      let last = -1;
      for (let i = N - 1; i >= 0; i--) if (Math.abs(p[i] - center) > eps) { last = i; break; }
      const li = document.createElement('li');
      const txt = last === -1 ? 'never outside the band'
        : last === N - 1 ? `still outside at <i>n</i> = ${N.toLocaleString('en-US')}`
          : `<i>n</i> = ${(last + 1).toLocaleString('en-US')}`;
      li.innerHTML = `<span style="color:${pal[k]}">●</span> path ${k + 1}: ${txt}`;
      ul.appendChild(li);
    });

    $('#lln-note').textContent = key === 'cauchy'
      ? 'Cauchy has no mean, so there is nothing to converge to: the average of n Cauchy variables is again Cauchy, the paths keep jumping, and the probability of leaving the band (around 0) does not shrink.'
      : 'Read the two panels together: the left one is an average over many paths at a fixed n (convergence in probability); the right one follows each single path forever (almost sure convergence). The SLLN says the right-hand list is finite for every path.';

    Object.assign(lln, { paths, N, eps, center, key, y: [lo, hi] });
    llnDraw();
  }

  function llnDraw() {
    if (!lln.paths) return;
    const { ctx, w, h } = prep($('#lln-canvas'));
    const { N, eps, center, y } = lln;
    const step = niceStep((y[1] - y[0]) / 5);
    const yticks = range(Math.ceil(y[0] / step) * step, y[1], step);
    const xticks = [1, 10, 100, 1000, 10000].filter(t => t <= N);
    const A = axes(ctx, w, h, {
      xmin: 1, xmax: N, ymin: y[0], ymax: y[1], xlog: true,
      xticks, yticks, xfmt: v => (v >= 1000 ? v / 1000 + 'k' : v),
      yfmt: v => (+v.toFixed(3)).toString(), xlabel: 'n (log scale)', ylabel: 'running mean'
    });
    clipTo(ctx, A);
    ctx.fillStyle = cssv('--acc-soft');
    ctx.fillRect(A.m.l, A.Y(center + eps), A.iw, A.Y(center - eps) - A.Y(center + eps));
    ctx.setLineDash([6, 5]); ctx.strokeStyle = cssv('--muted'); ctx.lineWidth = 1.2;
    ctx.beginPath(); ctx.moveTo(A.m.l, A.Y(center)); ctx.lineTo(A.m.l + A.iw, A.Y(center)); ctx.stroke();
    ctx.setLineDash([]);
    const pal = ['--c1', '--c2', '--c3', '--c4', '--c5', '--c6', '--c7', '--c8'].map(cssv);
    lln.paths.forEach((p, k) => {
      ctx.beginPath(); ctx.strokeStyle = pal[k]; ctx.lineWidth = 1.3; ctx.globalAlpha = 0.9;
      for (let i = 0; i < N; i++) { const x = A.X(i + 1), yy = A.Y(p[i]); i ? ctx.lineTo(x, yy) : ctx.moveTo(x, yy); }
      ctx.stroke();
    });
    ctx.globalAlpha = 1;
    ctx.restore();
  }

  $('#lln-eps').addEventListener('input', e => { $('#lln-eps-val').textContent = (+e.target.value).toFixed(2); llnCompute(); });
  $('#lln-dist').addEventListener('change', llnCompute);
  $('#lln-n').addEventListener('change', llnCompute);
  $('#lln-run').addEventListener('click', llnCompute);

  /* ==========================================================
     DEMO 2 — Central limit theorem
     ========================================================== */
  const clt = { hist: null, bins: 40, lo: -4, hi: 4 };

  function cltCompute() {
    const D = DISTS[$('#clt-dist').value];
    const n = +$('#clt-n').value;
    const R = 5000, B = clt.bins, lo = clt.lo, hi = clt.hi;
    const z = new Float64Array(R);
    const sn = Math.sqrt(n);
    for (let r = 0; r < R; r++) {
      let s = 0;
      for (let i = 0; i < n; i++) s += D.draw();
      z[r] = sn * (s / n - D.mu) / D.sd;
    }
    const counts = new Array(B).fill(0);
    const bw = (hi - lo) / B;
    let m1 = 0;
    for (let r = 0; r < R; r++) {
      m1 += z[r];
      const b = Math.floor((z[r] - lo) / bw);
      if (b >= 0 && b < B) counts[b]++;
    }
    m1 /= R;
    let m2 = 0, m3 = 0;
    for (let r = 0; r < R; r++) { const d = z[r] - m1; m2 += d * d; m3 += d * d * d; }
    m2 /= R; m3 /= R;
    const sorted = Float64Array.from(z).sort();
    let ks = 0;
    for (let i = 0; i < R; i++) {
      const F = Phi(sorted[i]);
      ks = Math.max(ks, Math.abs((i + 1) / R - F), Math.abs(i / R - F));
    }
    $('#clt-mean').textContent = fmt(m1);
    $('#clt-var').textContent = fmt(m2);
    $('#clt-skew').textContent = fmt(m3 / Math.pow(m2, 1.5));
    $('#clt-ks').textContent = fmt(ks);
    clt.hist = counts.map(c => c / (R * bw));
    cltDraw();
  }

  function cltDraw() {
    if (!clt.hist) return;
    const { ctx, w, h } = prep($('#clt-canvas'));
    const B = clt.bins, lo = clt.lo, hi = clt.hi, bw = (hi - lo) / B;
    const ymax = Math.max(0.45, Math.max(...clt.hist) * 1.1);
    const step = niceStep(ymax / 5);
    const A = axes(ctx, w, h, {
      xmin: lo, xmax: hi, ymin: 0, ymax,
      xticks: range(-4, 4, 1), yticks: range(0, ymax, step),
      yfmt: v => (+v.toFixed(2)).toString(), xlabel: 'z', ylabel: 'density'
    });
    clipTo(ctx, A);
    ctx.fillStyle = cssv('--c1'); ctx.globalAlpha = 0.45;
    clt.hist.forEach((d, i) => {
      const x0 = A.X(lo + i * bw), x1 = A.X(lo + (i + 1) * bw);
      ctx.fillRect(x0 + 0.5, A.Y(d), x1 - x0 - 1, A.Y(0) - A.Y(d));
    });
    ctx.globalAlpha = 1;
    ctx.beginPath(); ctx.strokeStyle = cssv('--c2'); ctx.lineWidth = 2.5;
    for (let i = 0; i <= 200; i++) {
      const x = lo + (hi - lo) * i / 200;
      i ? ctx.lineTo(A.X(x), A.Y(phi(x))) : ctx.moveTo(A.X(x), A.Y(phi(x)));
    }
    ctx.stroke();
    ctx.restore();
  }

  $('#clt-n').addEventListener('input', e => { $('#clt-n-val').textContent = e.target.value; cltCompute(); });
  $('#clt-dist').addEventListener('change', cltCompute);
  $('#clt-run').addEventListener('click', cltCompute);

  /* ==========================================================
     DEMO 3 — sqrt(n) and dependence
     ========================================================== */
  const sq = { n: [], sd: [], scn: 'indep', rho: 0.7 };

  function sqCompute() {
    const scn = $('#sq-scn').value;
    const R = 1000, K = 11, nmax = 1024, rho = 0.7, srho = Math.sqrt(1 - rho * rho);
    const cps = Array.from({ length: K }, (_, k) => Math.pow(2, k));
    const ss = new Float64Array(K);
    for (let r = 0; r < R; r++) {
      let x = randn(); const z0 = x; let s = 0, ci = 0;
      for (let t = 1; t <= nmax; t++) {
        if (t > 1) x = scn === 'indep' ? randn() : scn === 'ar1' ? rho * x + srho * randn() : z0;
        s += x;
        if (t === cps[ci]) { const m = s / t; ss[ci] += m * m; ci++; }
      }
    }
    const sd = Array.from(ss, v => Math.sqrt(v / R));
    $('#sq-table').innerHTML = cps.map((n, i) =>
      `<tr><td>${n}</td><td>${fmt(sd[i], 4)}</td><td>${fmt(1 / Math.sqrt(n), 4)}</td><td>${fmt(Math.sqrt(n) * sd[i], 3)}</td></tr>`).join('');
    $('#sq-note').textContent = {
      indep: 'Independent data: sd(X̄ₙ) follows σ/√n exactly (slope −½ on log–log axes), and √n · sd(X̄ₙ) stays ≈ 1 = σ. This is why √n is the right scale for the CLT.',
      ar1: 'Positively correlated data: the slope is still −½ for large n, but the curve sits higher — √n · sd(X̄ₙ) tends to σ·√((1+ρ)/(1−ρ)) ≈ 2.38. The naive standard error σ/√n would underestimate the real error by that factor.',
      dep: 'Perfectly dependent data: every X is the same number, so averaging gives nothing: sd(X̄ₙ) stays at σ = 1 for every n and √n · sd(X̄ₙ) grows like √n. Independence is what makes the errors cancel.'
    }[scn];
    Object.assign(sq, { n: cps, sd, scn, rho });
    sqDraw();
  }

  function sqDraw() {
    if (!sq.n.length) return;
    const { ctx, w, h } = prep($('#sq-canvas'));
    const A = axes(ctx, w, h, {
      xmin: 1, xmax: 1024, ymin: 0.02, ymax: 2, xlog: true, ylog: true,
      xticks: [1, 4, 16, 64, 256, 1024], yticks: [0.02, 0.05, 0.1, 0.2, 0.5, 1, 2],
      xlabel: 'n (log scale)', ylabel: 'sd of the sample mean (log)'
    });
    clipTo(ctx, A);
    const curve = (fn, color, dash, width) => {
      ctx.beginPath(); ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash);
      for (let i = 0; i <= 100; i++) {
        const n = Math.pow(1024, i / 100);
        i ? ctx.lineTo(A.X(n), A.Y(fn(n))) : ctx.moveTo(A.X(n), A.Y(fn(n)));
      }
      ctx.stroke(); ctx.setLineDash([]);
    };
    curve(n => 1 / Math.sqrt(n), cssv('--muted'), [7, 5], 1.8);
    if (sq.scn === 'ar1') {
      const c = Math.sqrt((1 + sq.rho) / (1 - sq.rho));
      curve(n => c / Math.sqrt(n), cssv('--c4'), [2, 4], 2);
    }
    ctx.beginPath(); ctx.strokeStyle = cssv('--c1'); ctx.lineWidth = 2;
    sq.n.forEach((n, i) => { const x = A.X(n), y = A.Y(sq.sd[i]); i ? ctx.lineTo(x, y) : ctx.moveTo(x, y); });
    ctx.stroke();
    ctx.fillStyle = cssv('--c1');
    sq.n.forEach((n, i) => { ctx.beginPath(); ctx.arc(A.X(n), A.Y(sq.sd[i]), 4, 0, 2 * Math.PI); ctx.fill(); });
    ctx.restore();
  }

  $('#sq-scn').addEventListener('change', sqCompute);
  $('#sq-run').addEventListener('click', sqCompute);

  /* ==========================================================
     EXAMPLE A — base rate / IDS alert
     ========================================================== */
  function alertUpdate() {
    const p = Math.pow(10, +$('#al-prior').value);
    const sens = +$('#al-sens').value / 100;
    const fpr = +$('#al-fpr').value / 100;
    $('#al-prior-val').textContent = (p * 100 < 1 ? (p * 100).toFixed(3) : (p * 100).toFixed(2)) + '%';
    $('#al-sens-val').textContent = (sens * 100).toFixed(1) + '%';
    $('#al-fpr-val').textContent = (fpr * 100).toFixed(1) + '%';
    const post = sens * p / (sens * p + fpr * (1 - p));
    $('#al-post').textContent = pct(post, post < 0.1 ? 2 : 1);
    const N = 100000;
    const attacks = N * p, tp = attacks * sens, fp = (N - attacks) * fpr;
    $('#al-tp').style.width = (tp / (tp + fp) * 100) + '%';
    $('#al-fp').style.width = (fp / (tp + fp) * 100) + '%';
    const r = v => Math.round(v).toLocaleString('en-US');
    $('#al-text').innerHTML =
      `Out of ${N.toLocaleString('en-US')} events there are about <b>${r(attacks)}</b> real attacks; <b>${r(tp)}</b> of them raise an alert. ` +
      `Among the ${r(N - attacks)} normal events, <b>${r(fp)}</b> raise a false alert. ` +
      `So only ${r(tp)} of the ${r(tp + fp)} alerts are real: ` +
      `P(attack | alert) = ${sens.toFixed(3)}·${p.toPrecision(2)} / (${sens.toFixed(3)}·${p.toPrecision(2)} + ${fpr.toFixed(3)}·${(1 - p).toPrecision(4)}) = <b>${pct(post, 2)}</b>.`;
  }
  ['#al-prior', '#al-sens', '#al-fpr'].forEach(s => $(s).addEventListener('input', alertUpdate));
  alertUpdate();

  /* ==========================================================
     EXAMPLE B — Beta–Binomial
     ========================================================== */
  function betaState() {
    const a = +$('#be-a').value, b = +$('#be-b').value, n = +$('#be-n').value;
    const kEl = $('#be-k');
    kEl.max = n;
    if (+kEl.value > n) kEl.value = n;
    const k = +kEl.value;
    $('#be-a-val').textContent = a; $('#be-b-val').textContent = b;
    $('#be-n-val').textContent = n; $('#be-k-val').textContent = k;
    return { a, b, n, k };
  }

  function betaDraw() {
    const { a, b, n, k } = betaState();
    const ap = a + k, bp = b + n - k;

    // numbers
    const priorMean = a / (a + b);
    const postMean = ap / (ap + bp);
    const mle = n > 0 ? k / n : NaN;
    const w = n / (a + b + n);
    $('#be-pm').textContent = fmt(priorMean);
    $('#be-mle').textContent = n > 0 ? fmt(mle) : '–';
    $('#be-post').textContent = fmt(postMean);
    $('#be-map').textContent = (ap > 1 && bp > 1) ? fmt((ap - 1) / (ap + bp - 2)) : 'at a boundary';
    $('#be-w').textContent = pct(w, 0);

    // grid, credible interval
    const G = 1000;
    const xs = [], post = [], prior = [], like = [];
    for (let i = 0; i < G; i++) {
      const x = (i + 0.5) / G;
      xs.push(x); prior.push(betaPdf(x, a, b)); post.push(betaPdf(x, ap, bp)); like.push(betaPdf(x, k + 1, n - k + 1));
    }
    let cum = 0, lo = 0, hi = 1, gotLo = false;
    const tot = post.reduce((s, v) => s + v, 0);
    for (let i = 0; i < G; i++) {
      cum += post[i] / tot;
      if (!gotLo && cum >= 0.025) { lo = xs[i]; gotLo = true; }
      if (cum >= 0.975) { hi = xs[i]; break; }
    }
    $('#be-cri').textContent = `[${fmt(lo)}, ${fmt(hi)}]`;
    if (n > 0) {
      const se = Math.sqrt(mle * (1 - mle) / n);
      $('#be-wald').textContent = se > 0
        ? `[${fmt(Math.max(0, mle - 1.96 * se))}, ${fmt(Math.min(1, mle + 1.96 * se))}]`
        : `degenerate at ${fmt(mle)}`;
    } else $('#be-wald').textContent = 'undefined (no data)';

    // plot
    const { ctx, w: cw, h: ch } = prep($('#be-canvas'));
    const ymax = Math.min(30, Math.max(...prior.filter(isFinite), ...post, ...like.filter(isFinite)) * 1.1) || 3;
    const step = niceStep(ymax / 5);
    const A = axes(ctx, cw, ch, {
      xmin: 0, xmax: 1, ymin: 0, ymax,
      xticks: range(0, 1, 0.1), yticks: range(0, ymax, step),
      xfmt: v => v.toFixed(1), yfmt: v => (+v.toFixed(1)).toString(), xlabel: 'p', ylabel: 'density'
    });
    clipTo(ctx, A);
    const draw = (arr, color, width, dash, fill) => {
      ctx.beginPath(); ctx.strokeStyle = color; ctx.lineWidth = width; ctx.setLineDash(dash);
      xs.forEach((x, i) => { const y = A.Y(Math.min(arr[i], ymax * 3)); i ? ctx.lineTo(A.X(x), y) : ctx.moveTo(A.X(x), y); });
      ctx.stroke(); ctx.setLineDash([]);
      if (fill) {
        ctx.lineTo(A.X(xs[G - 1]), A.Y(0)); ctx.lineTo(A.X(xs[0]), A.Y(0)); ctx.closePath();
        ctx.globalAlpha = 0.18; ctx.fillStyle = color; ctx.fill(); ctx.globalAlpha = 1;
      }
    };
    draw(prior, cssv('--c8'), 2.2, [6, 5], false);
    draw(post, cssv('--c1'), 3, [], true);
    draw(like, cssv('--c2'), 2.4, [2, 4], false);
    ctx.restore();
  }

  ['#be-a', '#be-b', '#be-n', '#be-k'].forEach(s => $(s).addEventListener('input', betaDraw));
  $$('.be-preset').forEach(btn => btn.addEventListener('click', () => {
    $('#be-a').value = btn.dataset.a; $('#be-b').value = btn.dataset.b; betaDraw();
  }));

  /* ==========================================================
     EXAMPLE C — naive Bayes spam / phishing filter
     ========================================================== */
  const WORDS = [
    { w: 'free', s: 0.30, h: 0.02, on: true },
    { w: 'winner', s: 0.10, h: 0.005, on: true },
    { w: 'urgent', s: 0.20, h: 0.03, on: false },
    { w: 'password', s: 0.15, h: 0.03, on: false },
    { w: 'invoice', s: 0.08, h: 0.06, on: false },
    { w: 'meeting', s: 0.02, h: 0.25, on: false },
    { w: 'thesis', s: 0.002, h: 0.04, on: false }
  ];
  const spamBody = $('#sp-table tbody');
  WORDS.forEach((o, i) => {
    const tr = document.createElement('tr');
    const lr = o.s / o.h;
    tr.innerHTML =
      `<td><label><input type="checkbox" data-i="${i}" ${o.on ? 'checked' : ''}> “${o.w}”</label></td>` +
      `<td>${o.s}</td><td>${o.h}</td><td>${lr >= 1 ? '×' + (+lr.toFixed(2)) : '×' + (+lr.toFixed(3))}</td>`;
    spamBody.appendChild(tr);
  });

  function spamUpdate() {
    const prior = +$('#sp-prior').value / 100;
    $('#sp-prior-val').textContent = Math.round(prior * 100) + '%';
    let odds = prior / (1 - prior);
    const steps = [`prior <b>${pct(prior, 1)}</b>`];
    $$('input[data-i]', spamBody).forEach(cb => {
      if (!cb.checked) return;
      const o = WORDS[+cb.dataset.i];
      odds *= o.s / o.h;
      steps.push(`after “${o.w}” <b>${pct(odds / (1 + odds), 1)}</b>`);
    });
    const post = odds / (1 + odds);
    $('#sp-post').textContent = pct(post, 1);
    $('#sp-chain').innerHTML = 'Sequential updating: ' + steps.join(' → ') +
      '. <br>Only words that appear in the e-mail are used here (a full Bernoulli naive-Bayes model would also use the absent words).';
  }
  $('#sp-prior').addEventListener('input', spamUpdate);
  spamBody.addEventListener('change', spamUpdate);
  spamUpdate();

  /* ---------- initial run & redraw registration ---------- */
  redraws.push(llnDraw, cltDraw, sqDraw, betaDraw);
  llnCompute();
  cltCompute();
  sqCompute();
  betaDraw();
})();
