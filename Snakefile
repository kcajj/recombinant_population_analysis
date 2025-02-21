import yaml
import pathlib
import numpy as np
import json
import sys


#  check that config file is provided
def check_configfile():
    try:
        msg = f"""
----------- loading config file -----------
--- run config:
{config["run_config"]}

--- alignments:
{config["alignments"]}

--- HMM parameters:
{config["HMM_parameters"]}

--- likelihood optimization of hyperparameter:
{config["optimization_recombination_parameter"]}

--- plots
{config["plots"]}
-------------------------------------------

"""
        print(msg)
    except:
        raise Exception(
            "config file not specified. Please specify with --configfile flag."
        )


check_configfile()

# make create log folder, required for cluster execution
pathlib.Path("log").mkdir(exist_ok=True)


# extract folder structure specification
run_config = config["run_config"]
in_fld = run_config["input"].removesuffix("/")
out_fld = run_config["output"].removesuffix("/")
HMM = run_config["HMM"]

# print run options
print("----------- run configuration ------------")
print("input folder:", in_fld)
print("\noutput folder:", out_fld)
print("\nHMM:", json.dumps(HMM, indent=2))
print("------------------------------------------")


include: "rules/HMM.smk"
include: "rules/plots.smk"


rule all:
    input:
        rules.HMM_all.input,
        rules.plot_all.input