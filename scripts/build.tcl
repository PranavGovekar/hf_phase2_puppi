set TOPDIR [pwd]
set PATH_CMSWS ""

if { $PATH_CMSWS == "" } {
	set PATH_CMSWS "${TOPDIR}/correlator-common/CMSSW_12_5_2_patch1/src"
}

file mkdir project
file mkdir outputs
cd projects 

open_project HF_CSIM
set_top algo_topIP1
add_files ${TOPDIR}/src/algo_topIP1.h -cflags "-I$PATH_CMSWS -I${TOPDIR}/correlator-common/puppi/firmware" 
add_files ${TOPDIR}/src/algo_top_IP.cpp -cflags "-I$PATH_CMSWS -I${TOPDIR}/correlator-common/puppi/firmware"
add_files ${TOPDIR}/correlator-common/puppi/firmware/linpuppi.cpp -cflags "-I${TOPDIR}/correlator-common/CMSSW_12_5_2_patch1/src"
add_files ${TOPDIR}/correlator-common/puppi/firmware/linpuppi.h -cflags "-I${TOPDIR}/correlator-common/CMSSW_12_5_2_patch1/src"
add_files ${TOPDIR}/correlator-common/puppi/firmware/linpuppi_bits.h -cflags "-I${TOPDIR}/correlator-common/CMSSW_12_5_2_patch1/src"

add_files -tb ${TOPDIR}/testbench/algo_top_tb_IP.cpp -cflags "-I${TOPDIR}/src"
open_solution "csim_solution" -flow_target vivado
set_part {xcvu13p-flga2577-3-e}
create_clock -period 10 -name default

exit

