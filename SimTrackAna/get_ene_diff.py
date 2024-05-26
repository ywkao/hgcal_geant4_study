#!/usr/bin/env python2
import json

observable = "mean" # energy scale
observable = "resolution" # relative resolution
# available keys of observable: "error_mean", "error_sigma", "mean", "resolution", "sigma", "uncertainty_resolution"

def check(ene, fname):
    print(">>> %s"%fname)
    with open(fname, 'r') as f:
        data = json.load(f)
        data_set1 = data[ene]["set1"]
        data_set2 = data[ene]["set2"]

    for e in [20, 60, 100, 175, 225, 300]:
        tag = "E%d" % e
        m1 = data_set1[tag][observable]
        m2 = data_set2[tag][observable]
        diff = (m1 - m2)/m2
        print "%s %s: %.3f" % (ene, tag, diff)
    
if __name__ == "__main__":
    check("MIP", "toolbox/test_resolution_unclustered.json")
    check("MIP", "toolbox/test_resolution_clustered.json")
    check("ENE", "toolbox/test_resolution_unclustered.json")
    check("ENE", "toolbox/test_resolution_clustered.json")
