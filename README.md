# hgcal_geant4_study

## Procedure of HGCAL GEANT4 simulation in cmssw

```
#--------------------------------------------------
# look for workflow of interest
#--------------------------------------------------
$ runTheMatrix.py -n -w upgrade | grep -i pythia | grep -v PU | grep D86 | grep -i electron

# 38602.0 2026D86+SingleElectronPt35_pythia8_GenSimHLBeamSpot+DigiTrigger+RecoGlobal+HARVESTGlobal

#--------------------------------------------------
# obtain scripts for the selected workflow
#--------------------------------------------------
$ time runTheMatrix.py -w upgrade -l 38602.0

#--------------------------------------------------
# proper modification
#--------------------------------------------------
$ mv cmdLog exe.sh
$ chmod +x exe.sh
$ vim exe.sh # cmsRun
$ vim SingleElectronPt35_pythia8_cfi_GEN_SIM.py # Nevents = 1000, Eta, etc.
$ vim step2_DIGI_L1TrackTrigger_L1_DIGI2RAW_HLT.py # Nevents = -1
$ vim step3_RAW2DIGI_L1Reco_RECO_RECOSIM_PAT_VALIDATION_DQM.py # Nevents = -1
$ vim step4_HARVESTING.py # Nevents = -1

#--------------------------------------------------
# execution
#--------------------------------------------------
$ time ./exe.sh
```

## Study of coverting ADC counts to MIPs

Details in SimTrackAna

## Reference
[1] https://hgcal.web.cern.ch/

