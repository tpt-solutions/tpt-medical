//! Anatomical taxonomy and models.
//!
//! Region enums cover the structures the simulation stack models. Variants
//! are intentionally closed: adding an anatomical region is a semver-minor
//! API change with an RFC, not a stringly-typed free-for-all, because
//! downstream regulatory tooling must enumerate regions exhaustively.

use std::collections::BTreeMap;

use tpt_med_geometry::{Aabb, Vec3};

/// Skeletal structures modelled by the stack.
#[derive(Debug, Clone, PartialEq)]
pub enum BoneType {
    /// Femur (thigh bone).
    Femur,
    /// Tibia (shin bone).
    Tibia,
    /// Fibula.
    Fibula,
    /// Patella (kneecap).
    Patella,
    /// Pelvis.
    Pelvis,
    /// Humerus (upper arm).
    Humerus,
    /// Radius.
    Radius,
    /// Ulna.
    Ulna,
    /// Vertebra at a named level, e.g. `L4`.
    Vertebra {
        /// Spinal level label (e.g. `"C3"`, `"T12"`, `"L5"`).
        level: String,
    },
    /// Skull (calvaria + base).
    Skull,
    /// Mandible.
    Mandible,
    /// User-defined structure (phantoms, explants, custom segments).
    Custom(String),
}

impl core::fmt::Display for BoneType {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            BoneType::Vertebra { level } => write!(f, "vertebra:{}", level.to_lowercase()),
            BoneType::Custom(name) => write!(f, "custom:{name}"),
            other => write!(f, "{}", format!("{other:?}").to_lowercase()),
        }
    }
}

/// Vascular structures.
#[derive(Debug, Clone, PartialEq)]
pub enum VesselType {
    /// Systemic artery (carotid, femoral, …).
    Artery {
        /// Clinical name, e.g. `"left_carotid"`.
        name: String,
    },
    /// Systemic vein.
    Vein {
        /// Clinical name.
        name: String,
    },
    /// Aorta segment.
    Aorta {
        /// Anatomical segment.
        segment: AortaSegment,
    },
    /// Coronary artery branch (LAD, RCA, …).
    Coronary {
        /// Branch name.
        branch: String,
    },
}

/// Aortic segments (De Bakey/SVS nomenclature, simplified).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AortaSegment {
    /// Aortic root and sinuses of Valsalva.
    Root,
    /// Ascending aorta.
    Ascending,
    /// Aortic arch.
    Arch,
    /// Descending thoracic aorta.
    DescendingThoracic,
    /// Abdominal aorta (infrarenal by default).
    Abdominal,
}

/// Hollow or solid organs relevant to planning.
#[derive(Debug, Clone, PartialEq)]
pub enum OrganType {
    /// Heart (whole-organ).
    Heart,
    /// Lung (left or right).
    Lung {
        /// `true` for the left lung.
        left: bool,
    },
    /// Liver.
    Liver,
    /// Kidney.
    Kidney {
        /// `true` for the left kidney.
        left: bool,
    },
    /// Brain.
    Brain,
}

/// Soft-tissue structures.
#[derive(Debug, Clone, PartialEq)]
pub enum SoftTissueType {
    /// Arterial wall (for HGO inflation models).
    ArterialWall,
    /// Skin and subcutis.
    Skin,
    /// Skeletal muscle.
    Muscle {
        /// Muscle name (e.g. `"vastus_medialis"`).
        name: String,
    },
    /// Articular cartilage.
    ArticularCartilage,
    /// Intervertebral disc.
    Disc {
        /// Adjacent spinal level.
        level: String,
    },
    /// Ligament.
    Ligament {
        /// Ligament name (e.g. `"acl"`).
        name: String,
    },
    /// Tendon.
    Tendon {
        /// Tendon name (e.g. `"achilles"`).
        name: String,
    },
}

/// Implant families (see also `tpt-med-stents`, `tpt-med-orthopedics`).
#[derive(Debug, Clone, PartialEq)]
pub enum ImplantType {
    /// Total hip arthroplasty components.
    TotalHip {
        /// Femoral stem family.
        stem: String,
        /// Acetabular cup family.
        cup: String,
    },
    /// Total knee arthroplasty components.
    TotalKnee {
        /// Femoral component family.
        femoral: String,
        /// Tibial component family.
        tibial: String,
    },
    /// Vascular stent.
    Stent {
        /// Nominal expanded diameter in mm.
        diameter_mm: f64,
        /// Nominal length in mm.
        length_mm: f64,
    },
    /// Fracture fixation plate.
    Plate {
        /// Bone the plate is applied to.
        bone: BoneType,
        /// Number of cortical screws.
        screw_count: u32,
    },
}

/// An anatomical region: the union of everything the pipeline can segment.
#[derive(Debug, Clone, PartialEq)]
pub enum AnatomicalRegion {
    /// Bony structure.
    Bone {
        /// Which bone.
        bone_type: BoneType,
    },
    /// Organ.
    Organ {
        /// Which organ.
        organ_type: OrganType,
    },
    /// Vessel lumen or vessel wall.
    Vessel {
        /// Which vessel.
        vessel_type: VesselType,
    },
    /// Implant volume.
    Implant {
        /// Which implant.
        implant_type: ImplantType,
    },
    /// Soft tissue.
    SoftTissue {
        /// Which tissue.
        tissue_type: SoftTissueType,
    },
}

impl AnatomicalRegion {
    /// Short stable key used in masks, labels and audit records.
    pub fn key(&self) -> String {
        match self {
            AnatomicalRegion::Bone { bone_type } => format!("bone:{}", bone_type),
            AnatomicalRegion::Organ { organ_type } => {
                format!("organ:{}", format!("{organ_type:?}").to_lowercase())
            }
            AnatomicalRegion::Vessel { vessel_type } => match vessel_type {
                VesselType::Artery { name } => format!("vessel:artery:{name}"),
                VesselType::Vein { name } => format!("vessel:vein:{name}"),
                VesselType::Aorta { segment } => {
                    format!("vessel:aorta:{}", format!("{segment:?}").to_lowercase())
                }
                VesselType::Coronary { branch } => format!("vessel:coronary:{branch}"),
            },
            AnatomicalRegion::Implant { implant_type } => {
                format!("implant:{}", format!("{implant_type:?}").to_lowercase())
            }
            AnatomicalRegion::SoftTissue { tissue_type } => {
                format!("tissue:{}", format!("{tissue_type:?}").to_lowercase())
            }
        }
    }
}

/// A named anatomical landmark in patient coordinates (e.g. medial femoral
/// epicondyle). Landmarks drive implant sizing and alignment planning.
#[derive(Debug, Clone, PartialEq)]
pub struct Landmark {
    /// Stable landmark key (e.g. `"femur:medial_epicondyle"`).
    pub key: String,
    /// Position in patient coordinates (mm).
    pub position: Vec3,
    /// Optional confidence in [0, 1] (automated landmark detection).
    pub confidence: Option<f64>,
}

/// Assembled anatomical model: regions with bounding extents plus landmarks.
#[derive(Debug, Clone, Default)]
pub struct AnatomicalModel {
    /// Regions keyed by [`AnatomicalRegion::key`].
    pub regions: BTreeMap<String, AnatomicalRegion>,
    /// Region bounding boxes in patient coordinates (mm).
    pub extents: BTreeMap<String, Aabb>,
    /// Anatomical landmarks.
    pub landmarks: Vec<Landmark>,
}

impl AnatomicalModel {
    /// Registers a region with its extent.
    pub fn add_region(&mut self, region: AnatomicalRegion, extent: Aabb) {
        let key = region.key();
        self.extents.insert(key.clone(), extent);
        self.regions.insert(key, region);
    }

    /// Looks up a region by key.
    pub fn region(&self, key: &str) -> Option<&AnatomicalRegion> {
        self.regions.get(key)
    }

    /// Adds or replaces a landmark (last write wins, per landmark key).
    pub fn add_landmark(&mut self, landmark: Landmark) {
        if let Some(existing) = self.landmarks.iter_mut().find(|l| l.key == landmark.key) {
            *existing = landmark;
        } else {
            self.landmarks.push(landmark);
        }
    }

    /// Looks up a landmark by key.
    pub fn landmark(&self, key: &str) -> Option<&Landmark> {
        self.landmarks.iter().find(|l| l.key == key)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn region_keys_are_stable_and_distinct() {
        let femur = AnatomicalRegion::Bone {
            bone_type: BoneType::Femur,
        };
        let vertebra = AnatomicalRegion::Bone {
            bone_type: BoneType::Vertebra { level: "L4".into() },
        };
        assert_eq!(femur.key(), "bone:femur");
        assert_eq!(vertebra.key(), "bone:vertebra:l4");
        assert_ne!(femur.key(), vertebra.key());
    }

    #[test]
    fn vessel_keys_include_names() {
        let carotid = AnatomicalRegion::Vessel {
            vessel_type: VesselType::Artery {
                name: "left_carotid".into(),
            },
        };
        assert_eq!(carotid.key(), "vessel:artery:left_carotid");
    }

    #[test]
    fn anatomical_model_upserts_landmarks() {
        let mut model = AnatomicalModel::default();
        model.add_landmark(Landmark {
            key: "femur:te_line_mid".into(),
            position: Vec3::new(1.0, 2.0, 3.0),
            confidence: Some(0.9),
        });
        model.add_landmark(Landmark {
            key: "femur:te_line_mid".into(),
            position: Vec3::new(4.0, 5.0, 6.0),
            confidence: Some(0.99),
        });
        assert_eq!(model.landmarks.len(), 1);
        assert_eq!(model.landmark("femur:te_line_mid").unwrap().position.x, 4.0);
    }

    #[test]
    fn model_regions_roundtrip() {
        let mut model = AnatomicalModel::default();
        let region = AnatomicalRegion::Implant {
            implant_type: ImplantType::Stent {
                diameter_mm: 6.0,
                length_mm: 40.0,
            },
        };
        let key = region.key();
        model.add_region(
            region.clone(),
            Aabb {
                min: Vec3::ZERO,
                max: Vec3::new(6.0, 6.0, 40.0),
            },
        );
        assert_eq!(model.region(&key), Some(&region));
    }
}
