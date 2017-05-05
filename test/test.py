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

TEST_DATA = {
    "asl_multiphase" : [
        {"repeats" : 1, "nph" : 8, "modfn" : "fermi", "alpha" : 66, "beta" : 21}, # Model options
        {"mag" : [10, 20, 40, 80, 160, 320],
         "phase" : [0, 0.1, 0.2, 0.3, 0.4, 0.5],
         "offset" : [500, 1000, 1500, 2000]}, # Parameter values - at most 3 can vary
        {"nt" : 8} # Other options
    ]
}

try:
    for model, test_data in TEST_DATA.items():
        rundata, params, kwargs = test_data
        self_test(model, rundata, params, noise=0.01, **kwargs)
except:
    traceback.print_exc()

