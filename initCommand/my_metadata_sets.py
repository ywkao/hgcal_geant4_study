#!/usr/bin/env python

parameters = {
    'set1' : {
        "directory"    : "batch_CloseByParticle_Photon_ERZRanges",
        "particle_gun" : "CloseByParticle_Photon_ERZRanges_cfi",
        "nevents"      : "10",
        "conditions"   : "auto:phase2_realistic_T15",
        "beamspot"     : "HGCALCloseBy",
        "geometry"     : "Extended2026D49",
        "era"          : "Phase2C9",
    },

    'set2' : {
        "directory"    : "batch_muon_gun_v1",
        "particle_gun" : "SingleMuFlatPt2To100_cfi",
        "nevents"      : "10",
        "conditions"   : "auto:phase2_realistic_T21",
        "beamspot"     : "HLLHC14TeV",
        "geometry"     : "Extended2026D86",
        "era"          : "Phase2C11I13M9",
    },

    # SingleElectronFlatPt2To100_cfi.py
    'set3' : {
        "directory"    : "batch_positron_gun_D83", # v15
        "particle_gun" : "SingleElectronFlatPt2To100_cfi",
        "nevents"      : "10",
        "conditions"   : "auto:phase2_realistic_T21",
        "beamspot"     : "HLLHC14TeV",
        "geometry"     : "Extended2026D83",
        "era"          : "Phase2C11I13M9",
    },

    'set4' : {
        "directory"    : "batch_positron_gun_D86", # v16
        "particle_gun" : "SingleElectronFlatPt2To100_cfi",
        "nevents"      : "10",
        "conditions"   : "auto:phase2_realistic_T21",
        "beamspot"     : "HLLHC14TeV",
        "geometry"     : "Extended2026D86",
        "era"          : "Phase2C11I13M9",
    },
}

