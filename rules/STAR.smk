include: "Snakefile_utils.py"
import os, sys

""" We are going to separate samples into groups of specified sizes """
""" Then we are going to perform STAR in each group while keeping the """
""" genome index in memory. This saves a lot of loading time """
rule star:
    input: SEEDSFILE
    output: rootdir/STAR_wrapper.done

rule star_group:
    input: 

