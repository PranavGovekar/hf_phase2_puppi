set TOPDIR [pwd]
set PATH_CMSWS "/home/bitsmith/CMSSW/CMSSW_12_5_2_patch1/src"

if { $PATH_CMSWS == "" } {
	set PATH_CMSWS "${TOPDIR}/correlator-common/CMSSW_12_5_2_patch1/src"
}

file mkdir project
cd project

open_project HF_algotopMultiCopy
set_top algo_topIP1
set_top packLinks
set_top unpackLinks
set_top compute
set_top algo_topIP1
add_files ${TOPDIR}/src/algo_topIP1.h -cflags "-I$PATH_CMSWS -I${TOPDIR}/correlator-common/puppi/firmware" 
add_files ${TOPDIR}/src/algo_top_IP.cpp -cflags "-I$PATH_CMSWS -I${TOPDIR}/correlator-common/puppi/firmware"
add_files ${TOPDIR}/correlator-common/puppi/firmware/linpuppi.cpp -cflags "-I$PATH_CMSWS"
add_files ${TOPDIR}/correlator-common/puppi/firmware/linpuppi.h -cflags "-I$PATH_CMSWS"
add_files ${TOPDIR}/correlator-common/puppi/firmware/linpuppi_bits.h -cflags "-I$PATH_CMSWS"

add_files -tb ${TOPDIR}/testbench/algo_top_tb_IP_all_events.cpp -cflags "-I${TOPDIR}/src -I$PATH_CMSWS -I${TOPDIR}/correlator-common/puppi/firmware"
open_solution "csim_solution" -flow_target vivado
set_part {xcvu13p-flga2577-1-e}
create_clock -period 2.77 -name default

exit

