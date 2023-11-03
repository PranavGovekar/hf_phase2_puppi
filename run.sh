#!/bin/bash

usage() {
    echo "Usage: $0 [-csm] [-cyn] [-build] [-rp]"
    echo "  -build  Set up Project"
    echo "  -csim   Run C Simulation"
    echo "  -csyn   Run C Synthesis"
    echo "  -rp     View Synthesis Report "
    exit 1
}

flag_provided=false

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -csim)
            echo "Running C Simulation...."
            vitis_hls -f scripts/run_csim.tcl 
            flag_provided=true
            ;;
        -csyn)
            echo "Running C Synthesis...."
            vitis_hls -f scripts/run_csynth.tcl 
            flag_provided=true
            ;;
        -build)
            echo "Setting up Project...."
            vitis_hls -f scripts/build.tcl 
            flag_provided=true
            ;;
        -rp)
            echo "Synthesis Report"
            cat projects/HF_CSIM/csim_solution/syn/report/algo_topIP1_csynth.rpt
            cp projects/HF_CSIM/csim_solution/syn/report/algo_topIP1_csynth.rpt ./report.rpt
            flag_provided=true
            ;;
        *)
            # Unknown flag
            usage
            ;;
    esac
    shift
done

# If no flags provided, show usage
if ! $flag_provided; then
    usage
fi
