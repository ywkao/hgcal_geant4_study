# Study of coverting ADC counts to MIPs

Initial steps for running code 
```
# cd to your $CMSSW_RELEASE_BASE/src
git clone git@github.com:ywkao/hgcal_geant4_study.git
cd hgcal_geant4_study/SimTrackAna/
cmsenv
time scram b -j 10

# modify proper input (step2.root) in python/digiHit_cfi.py before executing the following command
time cmsRun python/digiHit_cfg.py
```

The output should be error-free messages and a root file (geantoutput.root).

To implement the conversion from ADC counts to MIPs, relevant official code in cmssw is referred:

- [HGCalUncalibRecHitProducer.cc](https://github.com/cms-sw/cmssw/blob/master/RecoLocalCalo/HGCalRecProducers/plugins/HGCalUncalibRecHitProducer.cc#L21-L22)
    - C++ module for the creation of the un-calibrated rechits
    - Convert the input HGCAL Digis (ADC count) to rechits (MIPs)
    - Amplitude of rechits is expressed in terms of average number of MIPs
    - "algo" is specified in the configuration file: [HGCalUncalibRecHit\_cfi.py](https://github.com/cms-sw/cmssw/blob/master/RecoLocalCalo/HGCalRecProducers/python/HGCalUncalibRecHit_cfi.py#L71)
- [HGCalUncalibRecHitWorkerWeights.cc](https://github.com/cms-sw/cmssw/blob/master/RecoLocalCalo/HGCalRecProducers/plugins/HGCalUncalibRecHitWorkerWeights.cc)
    - Define methods: set / runHGCEE / runHGCHEsil / runHGCHEscint / runHGCHFNose
    - Push back uncalibrated rechit into corresponding containers
- [HGCalUncalibRecHitRecWeightsAlgo.h](https://github.com/cms-sw/cmssw/blob/master/RecoLocalCalo/HGCalRecAlgos/interface/HGCalUncalibRecHitRecWeightsAlgo.h#L51-L100)
    - Define a core method: makeRecHit
    - Evaluate amplitude and jitter

One key object is "worker\_", which is an object with auxiliary methods defined in the first reference code.
Related lines are implemented in python/digiHit\_cfg.py and plugins/digiHit.cc. 

However, if uncommenting this line, [plugins/digiHit.cc#L418](https://github.com/ywkao/hgcal_geant4_study/blob/main/SimTrackAna/plugins/digiHit.cc#L418),
the following fatal messages will appear:
```
----- Begin Fatal Exception 29-Dec-2021 06:36:13 CET-----------------------
An exception of category 'NoProxyException' occurred while
[0] Processing  Event run: 1 lumi: 1 event: 1 stream: 0
[1] Running path 'p'
[2] Calling method for module DigiSim/'prodEE_DigiSim'
Exception Message:
No data of type "HGCalGeometry" with label "HGCalEESensitive" in record "IdealGeometryRecord"
Please add an ESSource or ESProducer to your job which can deliver this data.
----- End Fatal Exception -------------------------------------------------
```

Investigation to the treatment is still on going.

[UPDATE] The problem is resolved by adding a line to python/digiHit\_cfg.py:
```
process.load('Configuration.Geometry.GeometryExtended2026D86Reco_cff')
```

## Reference
[1] https://hgcal.web.cern.ch/HitCalibration/hitCalibration/

