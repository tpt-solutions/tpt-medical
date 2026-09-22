# @tpt-solutions/stent-simulator

White-label, zero-cloud **stent deployment simulator** as a native web
component. The Nitinol superelastic deployment model (ASTM F2394-style
metrics) is pure Rust compiled to WebAssembly — every calculation runs in
the visitor's browser, so no geometry, pressure setting, or patient data
ever leaves the page.

```html
<script type="module" src="./tpt-stent-simulator.js"></script>

<tpt-stent-simulator
  brand-name="Acme Vascular"
  brand-url="https://acme.example"
  accent-color="#0e7c86"
  theme="light"
  wasm-base="../pkg/">
</tpt-stent-simulator>
```

## Engine prerequisite

The component imports the shared wasm-bindgen glue from `web/pkg/` (one
build per checkout, ~190 KB `.wasm`):

```console
# from the repository root
./scripts/build-web.sh        # macOS / Linux
# or
.\scripts\build-web.ps1       # Windows
```

Serve over HTTP from the repository root (ES modules and `fetch` do not
work from `file://`):

```console
python -m http.server 8080
# open http://localhost:8080/web/stent-simulator/
```

For npm distribution, publish this directory together with the `web/pkg/`
artifacts and point `wasm-base` at wherever you host them (any static
origin works, CORS-permissive).

## Attributes (branding)

| Attribute | Default | Purpose |
|---|---|---|
| `brand-name` | `TPT Solutions` | Header brand text |
| `brand-url` | — | Makes the brand text a link |
| `logo-url` | — | Optional logo image in the header |
| `accent-color` | `#2c7ef8` | Button, header rule, stent stroke, chart |
| `theme` | `light` | `light` or `dark` card palette |
| `report-note` | built-in disclaimer | Footer / report note line |
| `wasm-base` | `../pkg/` | Base URL of the wasm-bindgen glue (resolved against this module) |

## Attributes (simulation parameters)

| Attribute | Default | Unit |
|---|---|---|
| `expanded-diameter` | 6.0 | mm |
| `crimped-diameter` | 1.8 | mm |
| `n-crowns` | 12 | – |
| `crown-stiffness` | 0.5 | N/mm per crown |
| `lumen-diameter` | 4.6 | mm |
| `vessel-compliance` | 6.0 | mm/MPa (`D(p) = lumen + compliance·p`) |
| `pressure-kpa` | 100 | kPa |

Parameter attributes also work as two-way seeds: editing a field updates
its attribute; changing the attribute updates the field.

## Methods & properties

| Member | Description |
|---|---|
| `deploy()` | Runs the deployment + pressure sweep, updates the UI, returns the result (async) |
| `result` | Last result: `{ diameter, radialForce, contactPressure, recoil, dogboning }` or `null` |

## Events

| Event | `detail` |
|---|---|
| `tpt-deploy` | `{ params, result, sweep }` — one row per sweep pressure (60–160 kPa) |
| `tpt-error` | `{ message }` |

```js
sim.addEventListener("tpt-deploy", ({ detail }) => {
  console.log(detail.result.radialForce, "N");
});
```

## What it computes

One call to the shared `wasm_deploy_stent` binding evaluates the ring
deployment equilibrium (chronic outward force, wall contact pressure,
acute recoil, dogboning) plus a 6-point pressure sweep plotted inline —
the same model exercised by `cargo run -p tpt-med-examples --bin
stent-deployment`. Fidelity is the Level-1 radial-ring model described in
`rfcs/0004` and `docs/book/src/stents.md`; 3D FEM contact is the
documented upgrade path.

## Support

- Issues: https://github.com/tpt-solutions/tpt-medical/issues
- White-label packaging, SLA, and integration support: TPT Solutions
- License: MIT OR Apache-2.0 (dual), same as the workspace

**Research use only. Not FDA cleared for clinical decision making.**
