#!/bin/bash

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -init   Initialize Submodules"
    echo "  -build  Set up Project"
    echo "  -csim   Run C Simulation"
    echo "  -csyn   Run C Synthesis"
    echo "  -rsyn   Run RTL Synthesis"
    echo "  -impl   Run Implementation"
    echo "  -rp     View Synthesis Report"
    exit 1
}

flag_provided=false

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -init)
            echo "Initializing correlator-common...."
            date
            git submodule init
            git submodule update
            flag_provided=true
        ;;
        -csim)
            echo "Running C Simulation...."
            date
            start=`date +%s`
            vitis_hls -f scripts/run_csim.tcl 
            date
            end=`date +%s`
            echo Execution time was `expr $end - $start` seconds.
            flag_provided=true
            ;;
        -csyn)    
	    echo "Copying the codes for provenenace ! "
            mkdir -p  project/HF_algotopMultiCopy/code
	    cp -r src/* project/HF_algotopMultiCopy/code
            echo "Running C Synthesis...."
            date
            start=`date +%s`
            vitis_hls -f scripts/run_csynth.tcl 
            date
            end=`date +%s`
            echo Execution time was `expr $end - $start` seconds.
            flag_provided=true
            ;;
        -rsyn)
            echo "Running RTL Synthesis...."
            date
            start=`date +%s`
            vitis_hls -f scripts/run_rsynth.tcl 
            date
            end=`date +%s`
            echo Execution time was `expr $end - $start` seconds.
            flag_provided=true
            ;;
         -impl)
            echo "Running Implementation ...."
            date
            start=`date +%s`
            vitis_hls -f scripts/run_impl.tcl 
            date
            end=`date +%s`
            echo Execution time was `expr $end - $start` seconds.
            flag_provided=true
            ;;
        -build)
            echo "Setting up Project...."
            date
            start=`date +%s`
            sed -i -e '3i\#define REG_HF\' correlator-common/puppi/firmware/linpuppi.h
            vitis_hls -f scripts/build.tcl 
            date
            end=`date +%s`
            echo Execution time was `expr $end - $start` seconds.
            flag_provided=true
            ;;
        -rp)
            echo "Synthesis Report"
            date
            cat project/HF_algotopMultiCopy/csim_solution/syn/report/algo_topIP1_csynth.rpt
            cp project/HF_algotopMultiCopy/csim_solution/syn/report/algo_topIP1_csynth.rpt ./report.rpt
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
