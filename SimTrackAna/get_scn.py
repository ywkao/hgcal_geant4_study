#!/usr/bin/env python2
import json

observable = "mean" # energy scale
observable = "resolution" # relative resolution
# available keys of observable: "error_mean", "error_sigma", "mean", "resolution", "sigma", "uncertainty_resolution"

def get_scn(data):
    c, ec = 1.*data["fit_intercept"] , 1.*data["fit_error_intercept"]
    s, es = 1.*data["fit_slope"]     , 1.*data["fit_error_slope"]
    n, en = 1.*data["fit_noise"]       , 1.*data["fit_error_noise"]
    print "$%.2f \pm %.2f$ & $%.3f \pm %.3f$ & $%.1f \pm %.1f$\\\\" % (s, es, c, ec, n, en)
    #print "$(%.1f \pm %.1f)\\%%$ & $(%.1f \pm %.1f)\\%%$ & $%.1f \pm %.1f$\\\\" % (s, es, c, ec, n, en)

def check(fname):
    print(">>> %s"%fname)
    with open(fname, 'r') as f:
        data = json.load(f)
        get_scn(data["resolution_unclustered_FIT_set1_set2_MIPs_regression"]["set1"])
        get_scn(data["resolution_unclustered_FIT_set1_set2_MIPs_regression"]["set2"])
        get_scn(data["resolution_unclustered_FIT_set1_set2_MeV"]["set1"])
        get_scn(data["resolution_unclustered_FIT_set1_set2_MeV"]["set2"])
        get_scn(data["resolution_clustered_FIT_set1_set2_MIPs_regression"]["set1"])
        get_scn(data["resolution_clustered_FIT_set1_set2_MIPs_regression"]["set2"])
        get_scn(data["resolution_clustered_FIT_set1_set2_MeV"]["set1"])
        get_scn(data["resolution_clustered_FIT_set1_set2_MeV"]["set2"])

if __name__ == "__main__":
    check("toolbox/resolution_fit_snc_E.json")

# toolbox/resolution_fit_linear.json
# toolbox/resolution_fit_sc.json
# toolbox/resolution_fit_snc.json

