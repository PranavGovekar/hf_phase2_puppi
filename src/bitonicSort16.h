#ifndef _BITONICSORT16_H_
#define _BITONICSORT16_H_

#include <iostream>
#include "ap_int.h"
#include "pfCluster.h"

using namespace std;

typedef PFcluster din_t;
typedef ap_uint<6> dloop_t;

class GreaterSmaller
{
public:
    din_t greater, smaller;
};

void bitonicSort16(din_t in[16], din_t out[16]);

#endif

