set TOPDIR [pwd]

file mkdir project
cd project

open_project HF_PUPPI

set_top algo_topIP1
add_files ${TOPDIR}/src/algo_topIP1.h -cflags "-I${TOPDIR}/common"
add_files ${TOPDIR}/src/algo_top_IP.cpp -cflags "-I${TOPDIR}/common"
add_files ${TOPDIR}/src/bitonicSort16.cpp -cflags "-I${TOPDIR}/src"
add_files ${TOPDIR}/src/bitonicSort16.h
add_files ${TOPDIR}/src/pfCluster.h
add_files ${TOPDIR}/common/firmware/linpuppi.h -cflags "-I${TOPDIR}/common"
add_files ${TOPDIR}/common/firmware/linpuppi.cpp -cflags "-I${TOPDIR}/common"
add_files ${TOPDIR}/common/firmware/linpuppi_bits.h -cflags "-I${TOPDIR}/common"

add_files -tb ${TOPDIR}/testbench/patternFileTestbench.cpp -cflags "-I${TOPDIR}/common -I${TOPDIR}/src"
open_solution "csim_solution" -flow_target vivado
set_part {xcvu13p-flga2577-1-e}
create_clock -period 2.77 -name default

exit
