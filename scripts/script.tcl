############################################################
## This file is generated automatically by Vitis HLS.
## Please DO NOT edit it.
## Copyright 1986-2022 Xilinx, Inc. All Rights Reserved.
## Copyright 2022-2023 Advanced Micro Devices, Inc. All Rights Reserved.
############################################################
set TOPDIR [pwd]
file mkdir projects
cd projects 
open_project HF_CSIM
add_files ../include/algo_topIP1.h
add_files ../src/algo_top_IP.cpp -cflags "-I${TOPDIR}/include"
add_files -tb ../testbench/algo_top_tb_IP.cpp -cflags "-I${TOPDIR}/include"
open_solution "csim_solution" -flow_target vivado
set_part {xcvu13p-flga2577-2-e}
create_clock -period 10 -name default
#source "./PF_clusterizer/solution1/directives.tcl"
csim_design
#csynth_design
#cosim_design
#export_design -format ip_catalog
exit
