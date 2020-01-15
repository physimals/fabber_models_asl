#!/bin/env python
"""
This script generates a set of test images for ASL models - running Fabber with matching
options on this test data should return the expected parameter values!
"""

import os, sys
import traceback

import numpy as np
import nibabel as nib

FSLDIR = os.environ["FSLDIR"]
sys.path.insert(0, FSLDIR + "/lib/python")
from fabber import self_test

pi = 3.141592654

TEST_DATA = {
    "asl_multiphase" : [
        {"repeats" : 1, "nph" : 8, "modfn" : "fermi", "alpha" : 55, "beta" : 12}, # Model options
        {"mag" : [100, 200, 400],
         "phase" : [-pi, -pi/2, -pi/4, 0, pi/4, pi/2, pi],
         "offset" : [500, 1000]}, # Parameter values - at most 3 can vary
        {"nt" : 8} # Other options
    ]
}

try:
    for model, test_data in TEST_DATA.items():
        rundata, params, kwargs = test_data
        ret = self_test(model, rundata, params, noise=50, patchsize=10, save_input=True, save_output=True, **kwargs)
except:
    traceback.print_exc()

