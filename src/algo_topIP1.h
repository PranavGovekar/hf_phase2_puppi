#ifndef ALGO_TOPIP1_H
#define ALGO_TOPIP1_H

#include <iostream>
#include "ap_int.h"
#include <algorithm>
#include <utility>
#include <stdint.h>


#include "DataFormats/L1TParticleFlow/interface/layer1_objs.h"
#include "DataFormats/L1TParticleFlow/interface/pf.h"
#include "DataFormats/L1TParticleFlow/interface/puppi.h"
#include "../../dataformats/layer1_multiplicities.h"
#include "linpuppi_bits.h"

//-for-streaming-interface---------------------//
#include "hls_stream.h"                        //
#include "ap_axi_sdata.h"                      //
#define N_WORDS 9                              // 
#define BIT_WIDTH 64                           //
typedef ap_axis <BIT_WIDTH,0,0,1> axi_stream;  //
//---------------------------------------------//

#define N_INPUT_LINKS   6
#define N_OUTPUT_LINKS  6

#define N_PF_LINK 8
#define N_PUPPI_LINK 8
#define N_SECTORS 6
#define N_PF 48

using namespace std;
typedef ap_uint<10> loop;

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

void algo_topIP1(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]);


void AXIStream_wrapper(
    hls::stream<axi_stream> &inputStream0,
    hls::stream<axi_stream> &inputStream1,
    hls::stream<axi_stream> &inputStream2,
    hls::stream<axi_stream> &inputStream3,
    hls::stream<axi_stream> &inputStream4,
    hls::stream<axi_stream> &inputStream5,

    hls::stream<axi_stream> &outputStream0,
    hls::stream<axi_stream> &outputStream1,
    hls::stream<axi_stream> &outputStream2,
    hls::stream<axi_stream> &outputStream3,
    hls::stream<axi_stream> &outputStream4,
    hls::stream<axi_stream> &outputStream5
    );

#endif
