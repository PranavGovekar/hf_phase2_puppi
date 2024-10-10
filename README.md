
# HF_phase2_puppi
This repo documents the developments towards adding HF to the Phase2 L1 Correlato

### Getting Started 
Setup CMSSW as described in [correlator-common](https://gitlab.cern.ch/cms-cactus/phase2/firmware/correlator-common) 
If you already have CMSSW setup elsewhere set ``PATH_CMSWS`` in scripts/build.tcl to its path 

Run ``./run.sh``
```
Usage: ./run.sh [OPTIONS]
Options:
  -init       Initialize Submodules
  -build      Set up Project
  -csim       Run C Simulation
  -csyn       Run C Synthesis
  -rsyn       Run RTL Synthesis
  -impl       Run Implementation
  -csimrp     View C Synthesis Report
  -implrp     View Implementation Report
  -axis       Build AXIS Wrapper
  -zcu        Build AXIMM Wrapper

```
Before running C Simulation make a directory named inputs ``mkdir inputs``  and place all you input files inside.


### Helper Scripts

**Pattern File Generation**
Generate HF pattern files with a customizable firing pattern. Please edit the files to customize the Et distribution of the clusters. For geenrating different configuration with a specific file , please change the seed of the generator by speciying it in the argument line through `-s`

```
[]$ python scripts/makeHFRawEvent.py  -s 799 -t hfV0
```
**Valitation Plots**
Dependancies : `numpy` , `matplotlib`
Generates the plots for PF-Clusters,HFTowers,Jets in HF geometry.
```
# First pipe the CSIM output to a file.
./run.sh  -csim | tee output.log

# run the plotter
python3 makeValidationPlots.py  -i output.log --tag v1

```
TODO :
 * validation plots for regionization steps
 * validation plots for EGs
 * validation plots for Energy Sums
 * validation plots for puppi objects
