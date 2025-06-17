set TOPDIR [pwd]
set PATH_CMSWS "/home/pranav/CMSSW_12_5_2_patch1/src"

if { $PATH_CMSWS == "" } {
	set PATH_CMSWS "${TOPDIR}/correlator-common/CMSSW_12_5_2_patch1/src"
}

file mkdir project
cd project

open_project AXIS_Interface

set_top AXIStream_wrapper
add_files ${TOPDIR}/hls_wrappers/axisWrapper.h -cflags "-I$PATH_CMSWS -I${TOPDIR}/correlator-common/puppi/firmware -I${TOPDIR}/src"
add_files ${TOPDIR}/hls_wrappers/axisWrapper.cpp -cflags "-I$PATH_CMSWS -I${TOPDIR}/correlator-common/puppi/firmware -I${TOPDIR}/src"
add_files ${TOPDIR}/src/algo_topIP1.h -cflags "-I$PATH_CMSWS -I${TOPDIR}/correlator-common/puppi/firmware -I${TOPDIR}/src" 
add_files ${TOPDIR}/src/algo_top_IP.cpp -cflags "-I$PATH_CMSWS -I${TOPDIR}/correlator-common/puppi/firmware -I${TOPDIR}/src"
add_files ${TOPDIR}/correlator-common/puppi/firmware/linpuppi.h -cflags "-I$PATH_CMSWS -I${TOPDIR}/correlator-common/puppi/firmware -I${TOPDIR}/src"
add_files ${TOPDIR}/correlator-common/puppi/firmware/linpuppi.cpp -cflags "-I$PATH_CMSWS -I${TOPDIR}/correlator-common/puppi/firmware -I${TOPDIR}/src"
add_files ${TOPDIR}/correlator-common/puppi/firmware/linpuppi_bits.h -cflags "-I$PATH_CMSWS -I${TOPDIR}/correlator-common/puppi/firmware -I${TOPDIR}/src"


add_files -tb ${TOPDIR}/testbench/axiStream_tb.cpp -cflags "-I$PATH_CMSWS -I${TOPDIR}/correlator-common/puppi/firmware -I${TOPDIR}/src"
open_solution "csim_solution" -flow_target vivado
set_part {xcvu13p-flga2577-1-e}
create_clock -period 2.77 -name default

exit

