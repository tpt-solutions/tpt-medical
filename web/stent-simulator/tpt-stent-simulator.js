/**
 * tpt-stent-simulator — white-label, zero-cloud stent deployment web component.
 *
 * Pure-Rust Nitinol deployment model compiled to WebAssembly
 * (`tpt-med-wasm`); everything runs in the visitor's browser. Drop the
 * module and the `<tpt-stent-simulator>` element into any page, set the
 * branding attributes, and ship.
 *
 * Engine build (once per checkout):
 *   ../../scripts/build-web.ps1   (Windows)
 *   ../../scripts/build-web.sh    (macOS/Linux)
 *
 * Events: `tpt-deploy` (detail: { params, result, sweep }),
 *         `tpt-error`  (detail: { message }).
 * License: MIT OR Apache-2.0 (same as the workspace).
 */

const DEFAULTS = Object.freeze({
  "expanded-diameter": 6.0,
  "crimped-diameter": 1.8,
  "n-crowns": 12,
  "crown-stiffness": 0.5,
  "lumen-diameter": 4.6,
  "vessel-compliance": 6.0,
  "pressure-kpa": 100,
});

const SWEEP_KPA = [60, 80, 100, 120, 140, 160];

const STYLES = `
:host {
  display: block;
  --accent: #2c7ef8;
  --bg: #ffffff;
  --fg: #1b2733;
  --muted: #5d6c7c;
  --line: #dde4ec;
  --field: #f4f7fa;
  color: var(--fg);
  font: 14px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif;
}
:host([theme="dark"]) {
  --bg: #141b24;
  --fg: #dce5ee;
  --muted: #8fa3b8;
  --line: #27313e;
  --field: #0e141c;
}
.card {
  background: var(--bg); border: 1px solid var(--line); border-radius: 10px;
  overflow: hidden; box-shadow: 0 1px 3px rgba(15, 23, 42, .08);
}
header {
  display: flex; align-items: center; gap: 10px;
  padding: 12px 16px; border-bottom: 2px solid var(--accent);
}
header img { height: 26px; max-width: 130px; object-fit: contain; }
header .brand { font-weight: 700; letter-spacing: .01em; }
header .brand a { color: inherit; text-decoration: none; }
header .brand a:hover { text-decoration: underline; }
header .product { color: var(--muted); font-size: 13px; margin-left: auto; }
.status {
  font-size: 12px; color: var(--muted); padding: 6px 16px 0;
}
.status.ready { color: #1a9b6c; }
.status.error { color: #c4433b; }
.body { display: flex; flex-wrap: wrap; gap: 16px; padding: 12px 16px 16px; }
.controls { flex: 1 1 210px; min-width: 200px; }
figure { flex: 1 1 300px; min-width: 260px; margin: 0; }
figure svg { width: 100%; height: auto; display: block; }
figcaption { color: var(--muted); font-size: 12px; text-align: center; margin-top: 4px; }
.grid { display: grid; grid-template-columns: 1fr 1fr; gap: 2px 10px; }
.grid label { font-size: 12px; color: var(--muted); margin-top: 6px; }
.grid input {
  width: 100%; box-sizing: border-box; background: var(--field); color: var(--fg);
  border: 1px solid var(--line); border-radius: 5px; padding: 5px 7px; font: inherit;
}
button.deploy {
  margin-top: 12px; width: 100%; padding: 9px 12px; border: 0; border-radius: 6px;
  background: var(--accent); color: #fff; font: inherit; font-weight: 600; cursor: pointer;
}
button.deploy:hover { filter: brightness(1.08); }
button.deploy:disabled { opacity: .55; cursor: default; }
.metrics {
  display: grid; grid-template-columns: repeat(auto-fit, minmax(120px, 1fr));
  gap: 8px; padding: 0 16px 12px;
}
.metric {
  background: var(--field); border: 1px solid var(--line); border-radius: 8px;
  padding: 8px 10px;
}
.metric .k { font-size: 11px; color: var(--muted); text-transform: uppercase; letter-spacing: .06em; }
.metric .v { font-size: 17px; font-weight: 650; font-variant-numeric: tabular-nums; margin-top: 2px; }
.metric .v em { font-style: normal; font-size: 12px; font-weight: 500; color: var(--muted); }
.chart { padding: 0 16px; }
.chart svg { width: 100%; height: auto; display: block; }
.legend { display: flex; gap: 14px; font-size: 12px; color: var(--muted); margin-top: 2px; }
.legend i { display: inline-block; width: 14px; height: 3px; vertical-align: middle; margin-right: 5px; border-radius: 2px; }
.note {
  padding: 10px 16px 14px; font-size: 12px; color: var(--muted);
  border-top: 1px solid var(--line);
}
.note a { color: var(--accent); }
`;

const TEMPLATE = `
<style>${STYLES}</style>
<div class="card">
  <header>
    <img class="logo" alt="" hidden>
    <span class="brand"><a class="brand-link" href="#" tabindex="-1"><span class="brand-name"></span></a></span>
    <span class="product">Stent Deployment Simulator</span>
  </header>
  <div class="status">loading engine…</div>
  <div class="body">
    <div class="controls">
      <div class="grid">
        <div><label>Expanded ⌀ (mm)</label><input class="p" data-key="expanded-diameter" type="number" step="0.1" min="0"></div>
        <div><label>Crimped ⌀ (mm)</label><input class="p" data-key="crimped-diameter" type="number" step="0.1" min="0"></div>
        <div><label>Crowns</label><input class="p" data-key="n-crowns" type="number" step="1" min="3"></div>
        <div><label>Stiffness (N/mm)</label><input class="p" data-key="crown-stiffness" type="number" step="0.05" min="0"></div>
        <div><label>Lumen ⌀ (mm)</label><input class="p" data-key="lumen-diameter" type="number" step="0.1" min="0"></div>
        <div><label>Compliance (mm/MPa)</label><input class="p" data-key="vessel-compliance" type="number" step="0.5" min="0"></div>
        <div><label>Pressure (kPa)</label><input class="p" data-key="pressure-kpa" type="number" step="10" min="0"></div>
      </div>
      <button class="deploy" type="button" disabled>Deploy stent</button>
    </div>
    <figure>
      <svg class="fig" viewBox="0 0 320 300" role="img" aria-label="Stent deployment cross-section"></svg>
      <figcaption class="figcap">cross-section — dashed: crimped; thin: nominal; solid: deployed</figcaption>
    </figure>
  </div>
  <div class="metrics">
    <div class="metric"><div class="k">Equilibrium ⌀</div><div class="v" data-out="diameter">–</div></div>
    <div class="metric"><div class="k">Radial force</div><div class="v" data-out="force">–</div></div>
    <div class="metric"><div class="k">Contact pressure</div><div class="v" data-out="contact">–</div></div>
    <div class="metric"><div class="k">Acute recoil</div><div class="v" data-out="recoil">–</div></div>
    <div class="metric"><div class="k">Dogboning</div><div class="v" data-out="dogboning">–</div></div>
  </div>
  <div class="chart">
    <svg class="sweep" viewBox="0 0 320 130" role="img" aria-label="Radial force and diameter vs pressure"></svg>
    <div class="legend">
      <span><i class="l-force"></i>radial force [N]</span>
      <span><i class="l-diam"></i>equilibrium ⌀ [mm]</span>
    </div>
  </div>
  <div class="note"><span class="report-note"></span></div>
</div>
`;

function attrNumber(el, name, fallback) {
  const v = parseFloat(el.getAttribute(name));
  return Number.isFinite(v) ? v : fallback;
}

function ringPoints(cx, cy, rOut, rIn, crowns) {
  const pts = [];
  const n = Math.max(3, crowns | 0) * 2;
  for (let i = 0; i < n; i++) {
    const a = (i * Math.PI * 2) / n - Math.PI / 2;
    const r = i % 2 === 0 ? rOut : rIn;
    pts.push(`${(cx + r * Math.cos(a)).toFixed(2)},${(cy + r * Math.sin(a)).toFixed(2)}`);
  }
  return pts.join(" ");
}

function polyline(values, x0, y0, w, h) {
  if (!values.length) return "";
  const min = Math.min(...values), max = Math.max(...values);
  const span = max - min || 1;
  return values.map((v, i) => {
    const x = x0 + (values.length === 1 ? w / 2 : (i * w) / (values.length - 1));
    const y = y0 + h - ((v - min) / span) * h;
    return `${x.toFixed(1)},${y.toFixed(1)}`;
  }).join(" ");
}

export class TptStentSimulator extends HTMLElement {
  static observedAttributes = [
    "brand-name", "brand-url", "logo-url", "accent-color", "theme",
    "wasm-base", "report-note",
    ...Object.keys(DEFAULTS),
  ];

  #mod = null;
  #enginePromise = null;
  #anim = null;
  #lastResult = null;

  constructor() {
    super();
    this.attachShadow({ mode: "open" }).innerHTML = TEMPLATE;
    this.#fillParams();
  }

  connectedCallback() {
    this.#applyBranding();
    this.#wire();
    this.#ensureEngine().then(
      () => {
        this.#setStatus("engine ready — runs fully in-browser", "ready");
        this.#setEnabled(true);
        this.deploy();
      },
      (err) => this.#fail("engine not found: " + (err && err.message || err)),
    );
  }

  attributeChangedCallback() {
    this.#applyBranding();
    this.#fillParams();
  }

  /** Last deployment result: { diameter, radialForce, contactPressure, recoil, dogboning }. */
  get result() { return this.#lastResult; }

  /** Runs a deployment with the current parameters and updates the UI. */
  async deploy() {
    const btn = this.$("button.deploy");
    if (btn) btn.disabled = true;
    this.#setStatus("deploying…");
    try {
      const mod = await this.#ensureEngine();
      const p = this.#params();
      const r = mod.wasm_deploy_stent(
        p.expandedDiameter, p.crimpedDiameter, p.nCrowns, p.crownStiffness,
        p.lumenDiameter, p.pressureKpa / 1000, p.vesselCompliance,
      );
      const result = {
        diameter: r.diameter,
        radialForce: r.radial_force,
        contactPressure: r.contact_pressure,
        recoil: r.recoil,
        dogboning: r.dogboning,
      };
      if (typeof r.free === "function") r.free();

      const sweep = [];
      for (const kpa of SWEEP_KPA) {
        const s = mod.wasm_deploy_stent(
          p.expandedDiameter, p.crimpedDiameter, p.nCrowns, p.crownStiffness,
          p.lumenDiameter, kpa / 1000, p.vesselCompliance,
        );
        sweep.push({
          pressureKpa: kpa,
          diameter: s.diameter,
          radialForce: s.radial_force,
          contactPressure: s.contact_pressure,
          recoil: s.recoil,
        });
        if (typeof s.free === "function") s.free();
      }

      this.#lastResult = result;
      this.#renderMetrics(result);
      this.#renderSweep(sweep);
      this.#animateDeploy(result.diameter, p);
      this.#setStatus("deployed at " + p.pressureKpa.toFixed(0) + " kPa", "ready");
      this.dispatchEvent(new CustomEvent("tpt-deploy", {
        bubbles: true, composed: true, detail: { params: p, result, sweep },
      }));
      return result;
    } catch (err) {
      this.#fail(String((err && err.message) || err));
      this.dispatchEvent(new CustomEvent("tpt-error", {
        bubbles: true, composed: true, detail: { message: String(err) },
      }));
      return null;
    } finally {
      if (btn) btn.disabled = !this.#mod;
    }
  }

  $(sel) { return this.shadowRoot.querySelector(sel); }
  $$(sel) { return [...this.shadowRoot.querySelectorAll(sel)]; }

  #params() {
    const g = (key, fb) => attrNumber(this, key, fb);
    const field = (key, fb) => {
      const input = this.$(`input[data-key="${key}"]`);
      const v = input ? parseFloat(input.value) : NaN;
      return Number.isFinite(v) ? v : (Number.isFinite(g(key, NaN)) ? g(key, fb) : fb);
    };
    return {
      expandedDiameter: field("expanded-diameter", DEFAULTS["expanded-diameter"]),
      crimpedDiameter: field("crimped-diameter", DEFAULTS["crimped-diameter"]),
      nCrowns: Math.max(3, Math.round(field("n-crowns", DEFAULTS["n-crowns"]))),
      crownStiffness: field("crown-stiffness", DEFAULTS["crown-stiffness"]),
      lumenDiameter: field("lumen-diameter", DEFAULTS["lumen-diameter"]),
      vesselCompliance: field("vessel-compliance", DEFAULTS["vessel-compliance"]),
      pressureKpa: field("pressure-kpa", DEFAULTS["pressure-kpa"]),
    };
  }

  #fillParams() {
    for (const input of this.$$("input.p")) {
      if (document.activeElement === input) continue;
      const key = input.dataset.key;
      input.value = attrNumber(this, key, DEFAULTS[key]);
    }
  }

  #applyBranding() {
    const name = this.getAttribute("brand-name") || "TPT Solutions";
    const url = this.getAttribute("brand-url");
    const logo = this.getAttribute("logo-url");
    const accent = this.getAttribute("accent-color");
    const note = this.getAttribute("report-note");

    this.$(".brand-name").textContent = name;
    const link = this.$(".brand-link");
    if (url) { link.href = url; link.removeAttribute("tabindex"); }
    else { link.removeAttribute("href"); link.tabIndex = -1; }

    const img = this.$(".logo");
    if (logo) { img.src = logo; img.alt = name; img.hidden = false; }
    else { img.hidden = true; }

    if (accent) this.style.setProperty("--accent", accent);

    this.$(".report-note").textContent = note ||
      "Research use only. Not FDA cleared for clinical decision making. " +
      "Simulations run entirely in the visitor's browser — no geometry or " +
      "patient data is uploaded. Powered by tpt-medical (pure Rust → WASM).";
  }

  #wire() {
    this.$("button.deploy").addEventListener("click", () => this.deploy());
    for (const input of this.$$("input.p")) {
      input.addEventListener("change", () => {
        const v = parseFloat(input.value);
        if (Number.isFinite(v)) this.setAttribute(input.dataset.key, String(v));
        this.deploy();
      });
    }
  }

  #setStatus(msg, cls) {
    const el = this.$(".status");
    el.textContent = msg;
    el.className = "status" + (cls ? " " + cls : "");
  }

  #fail(msg) {
    this.#setStatus(msg, "error");
    this.#setEnabled(false);
  }

  #setEnabled(on) {
    this.$("button.deploy").disabled = !on;
    for (const input of this.$$("input.p")) input.disabled = !on;
  }

  async #ensureEngine() {
    if (this.#mod) return this.#mod;
    if (!this.#enginePromise) {
      const base = this.getAttribute("wasm-base") || "../pkg/";
      const url = new URL("tpt_med_wasm.js", new URL(base, import.meta.url));
      this.#enginePromise = import(url.href).then(async (mod) => {
        await mod.default();
        this.#mod = mod;
        return mod;
      });
    }
    return this.#enginePromise;
  }

  #renderMetrics(r) {
    const set = (key, value, unit) => {
      this.$(`[data-out="${key}"]`).innerHTML =
        `${value} <em>${unit}</em>`;
    };
    set("diameter", r.diameter.toFixed(3), "mm");
    set("force", r.radialForce.toFixed(2), "N");
    set("contact", (r.contactPressure * 1000).toFixed(1), "kPa");
    set("recoil", (100 * r.recoil).toFixed(2), "%");
    set("dogboning", (100 * r.dogboning).toFixed(2), "%");
  }

  #renderSweep(sweep) {
    const svg = this.$(".sweep");
    const x0 = 34, y0 = 8, w = 270, h = 96;
    const forces = sweep.map(s => s.radialForce);
    const diams = sweep.map(s => s.diameter);
    const kp = sweep.map(s => s.pressureKpa);
    const fMin = Math.min(...forces), fMax = Math.max(...forces);
    const dMin = Math.min(...diams), dMax = Math.max(...diams);
    const fPath = polyline(forces, x0, y0, w, h);
    const dPath = polyline(diams, x0, y0, w, h);
    const accent = getComputedStyle(this).getPropertyValue("--accent").trim() || "#2c7ef8";
    svg.innerHTML = `
      <rect x="${x0}" y="${y0}" width="${w}" height="${h}" fill="none"
            stroke="currentColor" opacity=".25" rx="4"/>
      <polyline points="${dPath}" fill="none" stroke="#39c9a0" stroke-width="2"/>
      <polyline points="${fPath}" fill="none" stroke="${accent}" stroke-width="2"/>
      <text x="${x0}" y="${y0 + h + 16}" font-size="11" fill="currentColor" opacity=".7">${kp[0]} kPa</text>
      <text x="${x0 + w}" y="${y0 + h + 16}" font-size="11" text-anchor="end"
            fill="currentColor" opacity=".7">${kp[kp.length - 1]} kPa</text>
      <text x="${x0 - 6}" y="${y0 + 8}" font-size="11" text-anchor="end"
            fill="currentColor" opacity=".7">${fMax.toFixed(1)} N</text>
      <text x="${x0 - 6}" y="${y0 + h}" font-size="11" text-anchor="end"
            fill="currentColor" opacity=".7">${fMin.toFixed(1)} N</text>
      <text x="${x0 + w + 4}" y="${y0 + 8}" font-size="11"
            fill="#39c9a0">${dMax.toFixed(2)}</text>
      <text x="${x0 + w + 4}" y="${y0 + h}" font-size="11"
            fill="#39c9a0">${dMin.toFixed(2)}</text>`;
    this.$(".legend .l-force").style.background = accent;
    this.$(".legend .l-diam").style.background = "#39c9a0";
  }

  #animateDeploy(targetDiameter, p) {
    if (this.#anim) cancelAnimationFrame(this.#anim);
    const from = p.crimpedDiameter;
    const to = targetDiameter;
    const t0 = performance.now();
    const dur = 650;
    const step = (now) => {
      const t = Math.min(1, (now - t0) / dur);
      const e = 1 - Math.pow(1 - t, 3);
      this.#drawFigure(from + (to - from) * e, p);
      if (t < 1) this.#anim = requestAnimationFrame(step);
      else this.#anim = null;
    };
    this.#drawFigure(from, p);
    this.#anim = requestAnimationFrame(step);
  }

  #drawFigure(diameterNow, p) {
    const svg = this.$(".fig");
    const cx = 160, cy = 150;
    const worldMax = Math.max(p.expandedDiameter, p.lumenDiameter +
      p.vesselCompliance * (p.pressureKpa / 1000)) * 1.25 || 1;
    const px = 120 / worldMax; // px per mm
    const lumenD = p.lumenDiameter + p.vesselCompliance * (p.pressureKpa / 1000);
    const strut = Math.max(0.12, 5 / px);
    const rDeploy = (diameterNow / 2) * px;
    const rIn = Math.max(2, rDeploy - strut * px);
    const accent = getComputedStyle(this).getPropertyValue("--accent").trim() || "#2c7ef8";
    const rCrimp = (p.crimpedDiameter / 2) * px;
    const rNominal = (p.expandedDiameter / 2) * px;
    const rLumen = (lumenD / 2) * px;
    const rWall = rLumen + Math.max(4, 0.8 * px);

    svg.innerHTML = `
      <circle cx="${cx}" cy="${cy}" r="${rWall.toFixed(1)}" fill="#e8dcd0" stroke="#cbb9a6"/>
      <circle cx="${cx}" cy="${cy}" r="${rLumen.toFixed(1)}" fill="#f7f9fc"/>
      <circle cx="${cx}" cy="${cy}" r="${rNominal.toFixed(1)}" fill="none"
              stroke="currentColor" stroke-width="1" opacity=".45"/>
      <circle cx="${cx}" cy="${cy}" r="${rCrimp.toFixed(1)}" fill="none"
              stroke="currentColor" stroke-width="1" stroke-dasharray="4 4" opacity=".6"/>
      <polygon points="${ringPoints(cx, cy, rDeploy, rIn, p.nCrowns)}"
               fill="none" stroke="${accent}" stroke-width="2.4"
               stroke-linejoin="round"/>
      <circle cx="${cx}" cy="${cy}" r="2" fill="${accent}"/>
      <text x="${cx}" y="${cy - 8}" font-size="11" text-anchor="middle"
            fill="currentColor" opacity=".75">⌀ ${diameterNow.toFixed(2)} mm</text>`;
  }
}

customElements.define("tpt-stent-simulator", TptStentSimulator);
