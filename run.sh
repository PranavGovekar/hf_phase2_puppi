#!/bin/bash

usage() {
    echo "Usage: $0 [OPTIONS]"
    echo "Options:"
    echo "  -build        Set up Project structure and generate test vectors"
    echo "  -csim         Run C Simulation"
    echo "  -plot         Generate validation plots from simulation results"
    echo "  -csyn         Run C Synthesis"
    echo "  -rsyn         Run RTL Synthesis"
    echo "  -impl         Run Implementation"
    echo "  -rst_sim      Reset simulation environment (keep project files)"
    echo "  -reset        Full project reset (clean all generated files)"
    echo "  -axis         Build AXIS Wrapper"
    echo "  -zcu          Build AXIMM Wrapper"
    echo "  -csimrp       View C Synthesis Report"
    echo "  -implrp       View Implementation Report"
    echo "  -gui          Launch Vitis HLS GUI"
    echo "  -gui_vivado   Launch Vivado GUI"
    exit 1
}

split_log() {
    input_file="./vitis_hls.log"
    output_dir="./outputs/csim_outputs/"

    if [ ! -f "$input_file" ]; then
        echo "Warning: No vitis_hls.log file found to split."
        return
    fi

    mkdir -p "$output_dir"

    echo "Splitting vitis_hls.log into parts..."

    count=$(awk -v RS="@@@NEW@@@" -v out="$output_dir" '
    {
        if (length($0) > 0) {
            filename = out "csim_" count ".txt"
            print "Writing to: " filename > "/dev/stderr"
            print $0 > filename
            count++
        }
    }
    END { print count }' count=0 "$input_file")

    if [ "$count" -gt 0 ]; then
        last_file="${output_dir}csim_$((count-1)).txt"
        if [ -f "$last_file" ]; then
            rm "$last_file"
        fi
    fi

    echo "All parts saved in $output_dir"
}
make_plots() {
    local input_dir="./outputs/csim_outputs"
    local script_path="./scripts/makeValidationPlots.py"

    for input_file in "$input_dir"/*.txt; do
        local base_name=$(basename "$input_file" .txt)
        
        python3 "$script_path" \
            -i "$input_file" \
            -t "${base_name}_tag" \
            --processST \
            -p \
            -c \

        echo "Plots generated for : $base_name"
    done

    echo "All plots made. Results saved to ./outputs"
}

flag_provided=false

while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -rst_sim)
            rm -rf ./data
            rm -rf ./outputs
            rm ./vitis_hls.log

            mkdir -p ./data/inputPatternFiles
            mkdir -p ./outputs/clusterSummary
            mkdir -p ./outputs/caloClusterSummary
            mkdir -p ./outputs/stJets

            echo "Simulation environment reset."

            flag_provided=true
            ;;
        -reset)
            rm -rf ./data
            rm -rf ./outputs
            rm ./vitis_hls.log
            rm -rf ./project

             echo "Project environment reset."

            flag_provided=true
            ;;
        -build)
            echo "Setting up Project...."
            date
            start=`date +%s`
            mkdir -p ./data/inputPatternFiles
            mkdir -p ./outputs/clusterSummary
            mkdir -p ./outputs/caloClusterSummary
            mkdir -p ./outputs/stJets
            vitis_hls -f scripts/build.tcl 

            date
            end=`date +%s`
            echo Execution time was `expr $end - $start` seconds.
            flag_provided=true
            ;;
        -csim)
            input_dir="./data/inputPatternFiles"
            
            if [ -d "$input_dir" ] && [ -n "$(ls -A "$input_dir" 2>/dev/null)" ]; then
                echo "Running C Simulation...."
                date
                start=$(date +%s)
                vitis_hls -f scripts/run_csim.tcl 
                date
                end=$(date +%s)
                echo "Execution time was $((end - start)) seconds."
            else
                echo "No input files found. Generating sample vector..."

                mkdir -p "$input_dir"

                if python3 scripts/makeHFRawEvent.py; then
                    if [ -n "$(ls -A "$input_dir" 2>/dev/null)" ]; then
                        echo "Generation successful. Running simulation..."
                        date
                        start=$(date +%s)
                        vitis_hls -f scripts/run_csim.tcl 
                        date
                        end=$(date +%s)
                        echo "Execution time was $((end - start)) seconds."
                    else
                        echo "Error: Vector generation completed but no files created!"
                    fi
                else
                    echo "Error: Test vector generation failed! Skipping simulation."
                fi
            fi
            
            flag_provided=true
            ;;
        -plot)
            echo "Making Validation Plots...."
            date
            start=$(date +%s)

            split_log
            make_plots

            date
            end=$(date +%s)
            echo "Execution time was $((end - start)) seconds."
            flag_provided=true
            ;;
        -csyn)    
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
            #source /home/bitsmith/Xilinx/Vivado/2023.1/settings64.sh
            date
            start=`date +%s`
            rm -rf exported_IP
            mkdir exported_IP
            vitis_hls -f scripts/run_impl.tcl 
            date
            end=`date +%s`
            echo Execution time was `expr $end - $start` seconds.
            flag_provided=true
            ;;
        -axis)
            echo "Setting up Project...."
            date
            start=`date +%s`
            vitis_hls -f scripts/build_AXIS.tcl 
            date
            end=`date +%s`
            echo Execution time was `expr $end - $start` seconds.
            flag_provided=true
            ;;
        -zcu)
            echo "Setting up Project...."
            date
            start=`date +%s`
            vitis_hls -f scripts/build_ZCU.tcl 
            date
            end=`date +%s`
            echo Execution time was `expr $end - $start` seconds.
            flag_provided=true
            ;;
        -csimrp)
            echo "Synthesis Report"
            date
            less project/HF_PUPPI/csim_solution/syn/report/algo_topIP1_csynth.rpt
            flag_provided=true
            ;;
        -implrp)
            echo "Implementation Report"
            date
            cat project/HF_PUPPI/csim_solution/impl/report/vhdl/algo_topIP1_export.rpt
            flag_provided=true
            ;;
	    -gui)
            echo "Launching Vitis HLS"
            date
            cd project/HF_PUPPI/
	    vitis_hls
	    cd ../..
            flag_provided=true
            ;;
	    -gui_vivado)
            echo "Launching Vivado"
            date
            cat project/HF_PUPPI/csim_solution/impl/vhdl
	    vivado project.xpr
	    cd ../../../../../
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
