# HF_phase2_puppi
This repo documents the developments towards adding HF to the Phase2 L1 Correlato

### Getting Started 

To build the project
```bash
git clone https://gitlab.cern.ch/pgovekar/hf_phase2_puppi.git
cd hf_phase2_puppi
git checkout clusterizer
./run.sh -build
```
---
## Usage
```bash
./run.sh [OPTIONS]
```

### Options

| **Option**  | **Description**                       |
|-------------|---------------------------------------|
| `-build`    | Set up the project environment.      |
| `-csim`     | Run C simulation.                    |
| `-csyn`     | Run C synthesis.                     |
| `-rsyn`     | Run RTL synthesis.                   |
| `-impl`     | Run FPGA implementation.             |
| `-csimrp`   | View C synthesis report.             |
| `-implrp`   | View implementation report.          |
| `-axis`     | Build an AXI Stream (AXIS) wrapper.  |
| `-zcu`      | Build an AXI Memory Mapped (AXIMM) wrapper. |

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

# run the plotter v0
python3 scripts/validateOutput.py  -i output.log -t v1
# old
python3 scripts/makeValidationPlots.py  -i output.log --tag v1

```

