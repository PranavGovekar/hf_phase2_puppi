
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
