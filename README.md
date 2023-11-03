
# HF_phase2_puppi
This repo documents the developments towards adding HF to the Phase2 L1 Correlato

### Getting Started 
Setup CMSSW as described in [correlator-common](https://gitlab.cern.ch/cms-cactus/phase2/firmware/correlator-common) 
If you already have CMSSW setup elsewhere set ``PATH_CMSWS`` in scripts/build.tcl to its path 

Run ``./run.sh``
```
Usage: ./run.sh [-csm] [-cyn] [-build] [-rp]
  -build  Set up Project
  -csim   Run C Simulation
  -csyn   Run C Synthesis
  -rp     View Synthesis Report 

```
Before running C Simulation make a directory named inputs ``mkdir inputs``  and place all you input files inside.
