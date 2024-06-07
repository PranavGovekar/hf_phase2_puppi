#ifndef ALGO_TOPIP1_H
#define ALGO_TOPIP1_H

#include "firmware/linpuppi.h"
#include <iostream>
#include "ap_int.h"
#include <algorithm>
#include <utility>
#include <stdint.h>
#include <hls_stream.h>
#include <hls_task.h>

#include "layer1_objs.h"
#include "pf.h"
#include "puppi.h"
#include "layer1_multiplicities.h"
#include "firmware/linpuppi_bits.h"

#define TOWERS_IN_ETA 12
#define TOWERS_IN_PHI 36
#define MIN_CLUSTER_SEED_ENERGY 5

#define LINK_WIDTH 220

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

class hftower{
    public:
    ap_uint<8> energy;
    ap_uint<2> fb;

    hftower(){
        energy = 0;
        fb = 0;
    }

    ap_uint<10> gettower(void){
    	ap_uint<10> data;
        data  =
	((ap_uint<10>)energy & 0xFF) |
	((ap_uint<10>)fb<<8 & 0x300) ;
	return data ;
    }

    void fillhftower(ap_uint<10> i){
    	this->energy = i.range(7, 0);
    	this->fb = i.range(9, 8);
    }

    hftower(ap_uint<10> i){
    	this->energy = i.range(7, 0);
    	this->fb = i.range(9, 8);
    }

    hftower(const hftower& rhs){
    energy=rhs.energy;
    fb=rhs.fb;
    }

    hftower& operator=(const hftower& rhs){
    this->energy=rhs.energy;
    this->fb=rhs.fb;
    return *this;
    }

    ap_uint<8> Energy(void) {return energy;}
    ap_uint<2> Fb(void) {return fb;}
};


void algo_topIP1(ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], ap_uint<LINK_WIDTH> link_out[N_OUTPUT_LINKS]);

#endif
