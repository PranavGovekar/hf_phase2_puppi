#ifndef ALGO_TOPIP1_H
#define ALGO_TOPIP1_H

#include <iostream>
#include "ap_int.h"
#include "ap_utils.h"
#include <algorithm>
#include <utility>
#include <stdint.h>

#include "caloObjects.h"
#include "pfClusterering.h"

void algo_topIP1(ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], ap_uint<576> link_out[10]);


#endif
