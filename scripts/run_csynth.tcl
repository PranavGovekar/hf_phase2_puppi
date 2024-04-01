cd project
open_project HF_algotopMultiCopy
open_solution "csim_solution" -flow_target vivado
#csim_design
csynth_design
#cosim_design
#export_design -format ip_catalog
exit

