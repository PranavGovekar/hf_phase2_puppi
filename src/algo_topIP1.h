#ifndef ALGO_TOPIP1_H
#define ALGO_TOPIP1_H

#include "linpuppi.h"
#include <iostream>
#include "ap_int.h"
#include <algorithm>
#include <utility>
#include <stdint.h>
#include <hls_stream.h>
#include<hls_task.h>

#include "DataFormats/L1TParticleFlow/interface/layer1_objs.h"
#include "DataFormats/L1TParticleFlow/interface/pf.h"
#include "DataFormats/L1TParticleFlow/interface/puppi.h"
#include "../../dataformats/layer1_multiplicities.h"
#include "linpuppi_bits.h"

#define N_INPUT_LINKS   6
#define N_OUTPUT_LINKS  6

#define N_PF_LINK 8
#define N_PUPPI_LINK 8
#define N_SECTORS 6
#define N_PF 48
#define N_EXTRA (NCALO - NNEUTRALS)

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

void Regionizer(const ap_uint<576> link_center,
			const ap_uint<576> link_left,
			const ap_uint<576> link_right,
			const ap_uint<3> sector,
			l1ct::HadCaloObj puppiIn[NCALO]);

void pack(	l1ct::PuppiObj pfselne[NNEUTRALS],
			ap_uint<576> &link_out);

void fillCenterLink(const ap_uint<576> &link,
					const l1ct::PFRegion &region,
					l1ct::HadCaloObj puppiIn[NCALO]);

void fillExtra(	const ap_uint<576> &link_left,
				const ap_uint<576> &link_right,
				const l1ct::PFRegion &region,
				const ap_uint<3> N_REGION,
				l1ct::HadCaloObj puppiIn[NCALO]);

void mergeSort(l1ct::HadCaloObj leftStream[N_EXTRA],
			l1ct::HadCaloObj rightStream[N_EXTRA],
			l1ct::HadCaloObj puppiIn[NCALO]);

void getInside(const ap_uint<576> &link,
				const int &phi_offset,
				const l1ct::PFRegion &region,
				l1ct::HadCaloObj outstream[N_EXTRA]);

void puppi(	l1ct::PFRegion region[N_SECTORS],
			l1ct::HadCaloObj puppiIn[N_SECTORS][NCALO],
			l1ct::PuppiObj pfselne[N_SECTORS][NNEUTRALS]);

#endif
