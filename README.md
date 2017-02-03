# Petersen_lab_tools
## Some general-use scripts for the Petersen lab
### This repo is maintained by [Mathew H. Evans](mailto:mathew.evans@manchester.ac.uk "I have the internets")

Below is a (probably incomplete) list of scripts in this folder. Look at the individual files for usage information.

To use these scripts, download the folder using the link at screen right ==>

A better idea is to install git (see how [HERE](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git#Installing-on-Windows)), then type in a command window/terminal:

`git clone https://github.com/mathewzilla/Petersen_lab_tools.git`

You will need git installed on your machine to be able to track/upload/download changes in the code anyway.

After that, if you find any bugs in the code open an Issue here:
[https://github.com/mathewzilla/Petersen_lab_tools/issues](https://github.com/mathewzilla/Petersen_lab_tools/issues)

Or by clicking the Issues tab at the right (it looks like an exclamation mark in a circle) ==> 
I will then see all the bugs in one place and will (hopefully) fix them in due course. When posting issues please provide as much info as possible, at a minimum the error message/a description of the problem and the name/location of the file where this problem occurred.

## Current content of this repo:
### chimera.m
Matlab script for post-processing whisker tracking data from the Janelia whisker tracker

### avi_tester.m
A script for testing different avi compression levels for whisker tracking

### boot_tuning.m
Bootstrapped tuning curves

### contact_detector.m
Work out distances from whiskers to poles

### contact_detector_parallel.m
Parallel computation version of contact_detector.m

### contact_estimation.m
Use distances from poles to whiskers from contact_detector.m to estimate when contact occurs

### extract_trials.m
Extract a subset of trials from an experiment based on some criteria (VPm data)

### poletracker.m
Find a pole in a whisker tracking video

### poletracker2.m
Updated/experimental version of poletracker.m

### poletracker_avi.m
.avi video version of poletracker.m (buggy)

### PTTH.m
Peri-touch-time-histogram generation code (VPm data)


Any problems email [Mat](mailto:mathew.evans@manchester.ac.uk "but don't waste my time")
