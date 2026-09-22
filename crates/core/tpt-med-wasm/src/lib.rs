//! WebAssembly bindings for the tpt-medical zero-cloud simulation stack.
//!
//! Exposes the imaging→mesh→solve pipeline and the stent deployment model
//! to JavaScript. All computation happens inside the browser sandbox: a
//! patient's DICOM bytes enter the module and only derived results (mesh
//! buffers, stress scalars) leave it.
//!
//! Build for the browser (see `scripts/build-web.ps1` / `.sh`):
//!
//! ```console
//! cargo build -p tpt-med-wasm --target wasm32-unknown-unknown --release
//! wasm-bindgen --out-dir web/pkg --target web target/wasm32-unknown-unknown/release/tpt_med_wasm.wasm
//! ```
//!
//! The generated `web/pkg` module is shared by `web/viewer` (end-to-end
//! CT → mesh → FEM demo) and `web/stent-simulator` (white-label stent
//! deployment component).

#![forbid(unsafe_code)]
#![warn(missing_docs)]

use wasm_bindgen::prelude::*;

/// Serializes a pipeline failure into a JS-visible error.
fn to_js<E: core::fmt::Display>(e: E) -> JsValue {
    JsValue::from_str(&e.to_string())
}

/// Frame helper: a payload of concatenated DICOM files, each prefixed with
/// a 4-byte little-endian length (the JS side builds this with DataView).
fn split_framed(payload: &[u8]) -> Vec<&[u8]> {
    let mut files = Vec::new();
    let mut off = 0usize;
    while off + 4 <= payload.len() {
        let len = u32::from_le_bytes([
            payload[off],
            payload[off + 1],
            payload[off + 2],
            payload[off + 3],
        ]) as usize;
        off += 4;
        if off + len > payload.len() {
            break;
        }
        files.push(&payload[off..off + len]);
        off += len;
    }
    files
}

/// Parses a DICOM series from raw file bytes, segments bone at
/// `threshold_hu`, and returns mesh statistics.
#[wasm_bindgen]
pub struct WasmMeshPipeline {
    node_count: usize,
    element_count: usize,
    mean_modulus: f64,
    max_modulus: f64,
    mesh: tpt_med_meshing::VoxelHexMesh,
}

#[wasm_bindgen]
impl WasmMeshPipeline {
    /// Builds the pipeline from framed DICOM bytes and an HU threshold.
    #[wasm_bindgen(constructor)]
    pub fn new(dicom_payload: Vec<u8>, threshold_hu: f64) -> Result<WasmMeshPipeline, JsValue> {
        let mut slices = Vec::new();
        for bytes in split_framed(&dicom_payload) {
            slices.push(tpt_med_dicom::DicomParser::parse_bytes(bytes).map_err(to_js)?);
        }
        let series = tpt_med_dicom::DicomSeries::from_slices(slices).map_err(to_js)?;
        let mask = tpt_med_meshing::SegmentationMask::threshold_hu(&series, threshold_hu);
        let mesh = tpt_med_meshing::MedicalMesher::default()
            .voxels_to_hex_mesh(&mask)
            .map_err(to_js)?;
        Ok(Self {
            node_count: mesh.nodes.len(),
            element_count: mesh.elements.len(),
            mean_modulus: mesh.mean_modulus(),
            max_modulus: mesh.max_modulus(),
            mesh,
        })
    }

    /// Number of mesh nodes.
    #[wasm_bindgen(getter)]
    pub fn node_count(&self) -> usize {
        self.node_count
    }

    /// Number of hexahedral elements.
    #[wasm_bindgen(getter)]
    pub fn element_count(&self) -> usize {
        self.element_count
    }

    /// Mean element Young's modulus (MPa).
    #[wasm_bindgen(getter)]
    pub fn mean_modulus(&self) -> f64 {
        self.mean_modulus
    }

    /// Maximum element Young's modulus (MPa).
    #[wasm_bindgen(getter)]
    pub fn max_modulus(&self) -> f64 {
        self.max_modulus
    }

    /// Flat node position buffer `[x0, y0, z0, x1, ...]` (mm).
    pub fn node_positions(&self) -> Vec<f32> {
        self.mesh
            .nodes
            .iter()
            .flat_map(|p| [p.x as f32, p.y as f32, p.z as f32])
            .collect()
    }

    /// Flat element connectivity: 8 node indices per hex, element-major.
    pub fn element_nodes(&self) -> Vec<u32> {
        self.mesh.elements.iter().flatten().copied().collect()
    }

    /// Per-element Young's modulus (MPa) for colour mapping.
    pub fn element_moduli(&self) -> Vec<f32> {
        self.mesh
            .materials
            .iter()
            .map(|m| m.youngs_modulus as f32)
            .collect()
    }
}

/// In-browser static FEM solve: pinned bottom face, traction on the top
/// band, returning displacement and von Mises statistics.
#[wasm_bindgen]
pub struct WasmStressResult {
    max_displacement: f64,
    max_von_mises: f64,
    mean_von_mises: f64,
    iterations: usize,
}

#[wasm_bindgen]
impl WasmStressResult {
    /// Maximum displacement magnitude (mm).
    #[wasm_bindgen(getter)]
    pub fn max_displacement(&self) -> f64 {
        self.max_displacement
    }

    /// Peak von Mises stress (MPa).
    #[wasm_bindgen(getter)]
    pub fn max_von_mises(&self) -> f64 {
        self.max_von_mises
    }

    /// Mean von Mises stress (MPa).
    #[wasm_bindgen(getter)]
    pub fn mean_von_mises(&self) -> f64 {
        self.mean_von_mises
    }

    /// CG iterations used.
    #[wasm_bindgen(getter)]
    pub fn iterations(&self) -> usize {
        self.iterations
    }
}

/// Runs the linear FEM solve on a pipeline's mesh.
///
/// Nodes on the bottom band (`z <= min_z`) are fully fixed; a total
/// compressive force of `load_newtons` is shared uniformly over the nodes
/// of the top band (`z >= min_z + 0.85·Δz`).
#[wasm_bindgen]
pub fn wasm_solve_stance_load(
    pipeline: &WasmMeshPipeline,
    load_newtons: f64,
) -> Result<WasmStressResult, JsValue> {
    let model = tpt_med_biomechanics::BiomechanicsModel::from_voxel_mesh(&pipeline.mesh);
    let (mut min_z, mut max_z) = (f64::INFINITY, f64::NEG_INFINITY);
    for n in &model.nodes {
        min_z = min_z.min(n.z);
        max_z = max_z.max(n.z);
    }
    let mut bc = tpt_med_biomechanics::BoundaryConditions::default();
    let mut fixed = Vec::new();
    let mut loaded = Vec::new();
    for (i, n) in model.nodes.iter().enumerate() {
        if n.z <= min_z + 1e-6 {
            fixed.push(i as u32);
        }
        if n.z >= min_z + 0.85 * (max_z - min_z) {
            loaded.push(i as u32);
        }
    }
    for i in fixed {
        bc.fix_nodes([i]);
    }
    let per = load_newtons / loaded.len().max(1) as f64;
    for i in loaded {
        bc.add_force(i, tpt_med_geometry::Vec3::new(0.0, 0.0, -per));
    }
    let r = model.solve(&bc, 1e-8, 20_000).map_err(to_js)?;
    Ok(WasmStressResult {
        max_displacement: r.max_displacement(),
        max_von_mises: r.max_von_mises(),
        mean_von_mises: r.mean_von_mises(),
        iterations: r.stats.iterations,
    })
}

/// In-browser stent deployment screening (radial-force model).
#[wasm_bindgen]
pub struct WasmStentResult {
    diameter: f64,
    radial_force: f64,
    contact_pressure: f64,
    recoil: f64,
    dogboning: f64,
}

#[wasm_bindgen]
impl WasmStentResult {
    /// Equilibrium diameter (mm).
    #[wasm_bindgen(getter)]
    pub fn diameter(&self) -> f64 {
        self.diameter
    }

    /// Radial force on the vessel (N).
    #[wasm_bindgen(getter)]
    pub fn radial_force(&self) -> f64 {
        self.radial_force
    }

    /// Mean wall contact pressure (MPa).
    #[wasm_bindgen(getter)]
    pub fn contact_pressure(&self) -> f64 {
        self.contact_pressure
    }

    /// Acute recoil fraction.
    #[wasm_bindgen(getter)]
    pub fn recoil(&self) -> f64 {
        self.recoil
    }

    /// Dogboning fraction (0 for the uniform ring model).
    #[wasm_bindgen(getter)]
    pub fn dogboning(&self) -> f64 {
        self.dogboning
    }
}

/// Deploys a stent ring into a compliant vessel.
///
/// * `expanded_diameter`, `crimped_diameter`, `lumen_diameter` — mm
/// * `crown_stiffness` — N/mm per crown
/// * `pressure_mpa` — intraluminal pressure at deployment (MPa)
/// * `vessel_compliance_mm_per_mpa` — pressure–diameter law slope,
///   `D(p) = lumen_diameter + compliance · p` (0 ⇒ rigid lumen)
#[wasm_bindgen]
pub fn wasm_deploy_stent(
    expanded_diameter: f64,
    crimped_diameter: f64,
    n_crowns: u32,
    crown_stiffness: f64,
    lumen_diameter: f64,
    pressure_mpa: f64,
    vessel_compliance_mm_per_mpa: f64,
) -> WasmStentResult {
    let stent = tpt_med_stents::StentModel {
        expanded_diameter,
        crimped_diameter,
        n_crowns,
        crown_stiffness,
    };
    let r = tpt_med_stents::simulate_deployment(
        &stent,
        &tpt_med_stents::NitinolParams::default(),
        |p_mpa| lumen_diameter + vessel_compliance_mm_per_mpa * p_mpa,
        tpt_med_units::Pressure::from_mpa(pressure_mpa),
    );
    WasmStentResult {
        diameter: r.diameter,
        radial_force: r.radial_force,
        contact_pressure: r.contact_pressure,
        recoil: r.recoil,
        dogboning: r.dogboning,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn framed_payload_splits_into_files() {
        let mut payload = Vec::new();
        for bytes in [b"abc".as_slice(), b"defghij"] {
            payload.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
            payload.extend_from_slice(bytes);
        }
        let files = split_framed(&payload);
        assert_eq!(files, vec![&b"abc"[..], &b"defghij"[..]]);
    }

    #[test]
    fn framed_payload_tolerates_trailing_garbage() {
        let mut payload = Vec::new();
        payload.extend_from_slice(&1u32.to_le_bytes());
        payload.push(b'x');
        payload.extend_from_slice(&99u32.to_le_bytes()); // length overruns EOF
        payload.extend_from_slice(&[1, 2]);
        let files = split_framed(&payload);
        assert_eq!(files, vec![&b"x"[..]]);
    }

    #[test]
    fn deploy_produces_contact_metrics_and_dogboning() {
        // Mirrors the Phase 5 example: 6 mm stent, compliant artery
        // D = 4.6 + 6.0·p, deployed at 0.1 MPa.
        let r = wasm_deploy_stent(6.0, 1.8, 12, 0.5, 4.6, 0.1, 6.0);
        assert!(r.radial_force > 0.0);
        assert!(r.contact_pressure > 0.0);
        assert!(r.recoil > 0.0 && r.recoil < 0.2);
        assert!(r.diameter > 4.6 && r.diameter <= 6.0);
        assert_eq!(r.dogboning, 0.0);
    }

    #[test]
    fn rigid_vessel_yields_no_contact_when_stent_undersized() {
        let r = wasm_deploy_stent(4.0, 1.5, 8, 0.5, 6.0, 0.01, 0.0);
        assert_eq!(r.radial_force, 0.0);
        assert_eq!(r.diameter, 4.0);
    }

    #[test]
    fn synthetic_ct_mesh_solves_stance_load() {
        let dir = concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/../../../test-data/dicom/synthetic_ct"
        );
        let mut entries: Vec<_> = std::fs::read_dir(dir)
            .expect("synthetic CT test data present")
            .filter_map(|e| e.ok())
            .map(|e| e.path())
            .filter(|p| p.extension().map(|x| x == "dcm").unwrap_or(false))
            .collect();
        entries.sort();
        assert!(!entries.is_empty(), "synthetic CT slices found");

        let mut payload = Vec::new();
        for path in entries {
            let bytes = std::fs::read(&path).expect("read slice");
            payload.extend_from_slice(&(bytes.len() as u32).to_le_bytes());
            payload.extend_from_slice(&bytes);
        }

        let pipeline = WasmMeshPipeline::new(payload, 200.0).expect("pipeline builds");
        assert!(pipeline.element_count() > 0);
        assert_eq!(pipeline.element_nodes().len(), 8 * pipeline.element_count());
        assert_eq!(pipeline.element_moduli().len(), pipeline.element_count());
        assert_eq!(pipeline.node_positions().len(), 3 * pipeline.node_count());

        let r = wasm_solve_stance_load(&pipeline, 100.0).expect("solve succeeds");
        assert!(r.max_displacement() > 0.0);
        assert!(r.max_von_mises() > 0.0);
        assert!(r.mean_von_mises() <= r.max_von_mises() + 1e-9);
        assert!(r.iterations() > 0 && r.iterations() <= 20_000);
    }
}
