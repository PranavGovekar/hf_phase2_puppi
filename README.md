# HF_phase2_puppi

This repository documents the developments towards adding HF (Hadronic Forward) calorimeter information to the Phase-2 Level-1 Correlator Trigger.

---

## Getting Started

Clone the repository:

```bash
git clone https://gitlab.cern.ch/pgovekar/hf_phase2_puppi.git
cd hf_phase2_puppi
```

To build and launch the project in Vitis HLS:

```bash
source /opt/Xilinx/Vivado/2022.1/settings64.sh
./run.sh -build
./run.sh -gui
```

Once Vitis HLS is launched:
- Click on **Open Project**
- Then, click on **Open** (no path changes required)

---

## Generate IP

To run HLS synthesis and implementation and generate the IP:

```bash
./run.sh -build -csyn -impl
```

---

## Validation Plots

### Step 1: Generate Input Test Vectors

Use the `makeHFRawEvent.py` script to generate customizable input test vectors.

You can choose any combination of cluster density using the `-m` option. The available modes are:
- `s` / `small`
- `m` / `medium`
- `l` / `large`
- `e` / `extreme`

The total number of vectors (given by the second value in the `-s` argument) is evenly divided among the selected modes.

**Example:**

```bash
python3 scripts/makeHFRawEvent.py -m s,m,l -s 20 6 -p iv
```

This generates 6 input files with seed offset 20 and 2 vectors each for small (`s`), medium (`m`), and large (`l`) density of clusters.

---

### Step 2: Run Simulation and Plot Results

Pipe the simulation output and generate validation plots:

```bash
./run.sh -csim -plot
```

The plots for PFClusters, HFTowers, and Jets in HF geometry are saved in the `./outputs` directory.

---

## Usage

### `./run.sh` Options

| **Option**      | **Description**                                      |
|------------------|------------------------------------------------------|
| `-build`         | Set up project structure and generate test vectors.  |
| `-csim`          | Run C simulation.                                    |
| `-plot`          | Generate validation plots from simulation results.   |
| `-csyn`          | Run C synthesis.                                     |
| `-rsyn`          | Run RTL synthesis.                                   |
| `-impl`          | Run FPGA implementation.                             |
| `-rst_sim`       | Reset simulation environment (keep project files).   |
| `-reset`         | Full project reset (clean all generated files).      |
| `-axis`          | Build AXI Stream (AXIS) wrapper.                     |
| `-zcu`           | Build AXI Memory Mapped (AXIMM) wrapper.            |
| `-csimrp`        | View C synthesis report.                             |
| `-implrp`        | View implementation report.                          |
| `-gui`           | Launch Vitis HLS project GUI.                        |
| `-gui_vivado`    | Launch Vivado project GUI.                           |
---

### `makeHFRawEvent.py` Options

| **Option**         | **Description**                                                    |
|--------------------|--------------------------------------------------------------------|
| `-p PREFIX`        | Prefix for the output files.                                       |
| `-m MODES`         | Comma-separated list of cluster density modes (`s`, `m`, `l`, `e`).        |
| `-s SEED N`        | Seed offset and number of vectors to generate.                  |

**Example Usage:**

```bash
python3 scripts/makeHFRawEvent.py -m s,m,l -s 20 6 -p iv
```

This creates 6 vectors with 3 cluster density types (`s`, `m`, `l`), generating 2 of each type.
