#ifndef PF_CLUSTER_H
#define PF_CLUSTER_H

#include "ap_int.h"

class PFcluster {
public:
    ap_uint<12> ET;
    ap_uint<5> Eta;
    ap_uint<8> Phi;
    ap_uint<39> Spare;
    ap_uint<64> all;

    PFcluster() {
        ET = 0;
        Eta = 0;
        Phi = 0;
        Spare = 0;
    }

    void getPFcluster(ap_uint<64> i) {
        this->ET = i.range(11,0);
        this->Eta = i.range(16,12);
        this->Phi = i.range(24,17);
        this->Spare = i.range(63,25);
    }

    void getdata() {
        all = ET | ((ap_uint<64>)Eta << 12) | ((ap_uint<64>)Phi << 17) | ((ap_uint<64>)Spare << 25) ;
    }

    ap_uint<64> data() {
        ap_uint<64> out = ET | ((ap_uint<64>)Eta << 12) | ((ap_uint<64>)Phi << 17) | ((ap_uint<64>)Spare << 25) ;
        return out ;
    }

    PFcluster(const PFcluster& rhs) {
        ET = rhs.ET ;
        Eta = rhs.Eta ;
        Phi = rhs.Phi ;
        Spare = rhs.Spare ;
        all = rhs.all ;
    }

    PFcluster& operator=(const PFcluster& rhs) {
        this->ET = rhs.ET ;
        this->Eta = rhs.Eta ;
        this->Phi = rhs.Phi ;
        this->Spare = rhs.Spare ;
        this->all = rhs.all ;
        return *this ;
    }

};

#endif
