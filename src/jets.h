#ifndef JETS_H
#define JETS_H

#include "ap_int.h"

class jets
{
public:
    ap_uint<12> ET;
    ap_uint<3> Eta;
    ap_uint<7> Phi;
    ap_uint<14> seedET;

    jets()
    {
        ET = 0;
        Eta = 0;
        Phi = 0;
        seedET = 0;
    }
    ap_uint<64> data()
    {
        ap_uint<64> out = ET | ((ap_uint<64>)Eta << 12) | ((ap_uint<64>)Phi << 15) | ((ap_uint<64>)seedET << 22) ;
        return out ;
    }

};

#endif
