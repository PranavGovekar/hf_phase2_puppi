#include "linpuppi.h"
#include "algo_topIP1.h"
#include <cstdint>

void compute(const ap_uint<576> &link_center,
			const ap_uint<576> &link_left,
			const ap_uint<576> &link_right,
			const l1ct::PFRegion &region,
			ap_uint<576> &link_out
				){
	l1ct::HadCaloObj puppiIn[NCALO];
#pragma HLS ARRAY_PARTITION dim=1 type=complete variable=puppiIn
	l1ct::PuppiObj pfselne[NNEUTRALS];
#pragma HLS ARRAY_PARTITION dim=1 type=complete variable=pfselne


	const int N_REGION = (region.hwPhiCenter - 5)/12;

	fillCenterLink(link_center, region, puppiIn);
	fillExtra(link_left, link_right, region, N_REGION, puppiIn);

#ifndef __SYNTHESIS__
	for( loop j=0 ;j < NCALO ; j++){
		std::cout<<"HCALOobj "<<" , "<<j
		  <<" : phi "<<puppiIn[j].hwPhi
		  <<" : eta "<<puppiIn[j].hwEta
		  <<" : et  "<<puppiIn[j].hwPt
		  <<"\n";
	}
#endif

	fwdlinpuppi(region, puppiIn, pfselne);

	pack(pfselne, link_out);
}


void pack(	l1ct::PuppiObj pfselne[NNEUTRALS],
			ap_uint<576> &link_out){
	const ap_uint<8> BW = 64;
	ap_uint<12> start = 0;

    for(loop idx=0; idx<NNEUTRALS; idx++) {
        link_out.range(start+(BW-1),start) = pfselne[idx].pack();

        start = start + BW;
    }
}


void fillCenterLink(const ap_uint<576> &link,
					const l1ct::PFRegion &region,
					l1ct::HadCaloObj puppiIn[NCALO]
				){
	ap_uint<576> word = link;
	const int eta_offset = region.hwEtaCenter;
	const int phi_offset = region.hwPhiCenter;

	for(loop j=0; j<N_PF_LINK; j++) {
		puppiIn[j].hwPt = word.range(11,0);
		puppiIn[j].hwEta = (word.range(16,12) - eta_offset);
		puppiIn[j].hwPhi = (word.range(24,17) - phi_offset);
		word=word>>64;

	}
}

void fillExtra(	const ap_uint<576> &link_left,
				const ap_uint<576> &link_right,
				const l1ct::PFRegion &region,
				const int N_REGION,
				l1ct::HadCaloObj puppiIn[NCALO]
				){

	int region_l = (N_REGION==0) 				? (N_INPUT_LINKS-1):(N_REGION-1);
	int region_r = (N_REGION==N_INPUT_LINKS-1) 	? (0):(N_REGION+1);

	int phi_offset_l = 12*N_REGION + 5;
	int phi_offset_r = 12*N_REGION + 5;

	if (  N_REGION == 0              and    region_l == (N_INPUT_LINKS-1)) phi_offset_l+=72;
	if (  N_REGION == (N_SECTORS-1)  and    region_r == 0                ) phi_offset_r-=72;

	l1ct::HadCaloObj leftStream[4];
	l1ct::HadCaloObj rightStream[4];

	for(loop idx=0; idx<4; idx++) {
		leftStream[idx].clear();
		rightStream[idx].clear();
	}

	getStream(link_left, phi_offset_l, region, leftStream);
	getStream(link_right, phi_offset_r, region, rightStream);

	mergeSort(leftStream, rightStream, puppiIn);

}

void mergeSort(	l1ct::HadCaloObj leftStream[4],
				l1ct::HadCaloObj rightStream[4],
				l1ct::HadCaloObj puppiIn[NCALO]
	){

	l1ct::HadCaloObj left;
	l1ct::HadCaloObj right;

	int idx_left = 0;
	int idx_right = 0;

	left = leftStream[idx_left];
	right = rightStream[idx_right];

	for(loop idx=0; idx<4; idx++) {
		if(left.hwPt < right.hwPt){
			puppiIn[8+idx] = right;
			idx_right = idx_right + 1;
			right = rightStream[idx_right];
		}
		else {
			puppiIn[8+idx] = left;
			idx_left = idx_left + 1;
			left = leftStream[idx_left];
		}
	}
}

//void clearStream(hls::stream<l1ct::HadCaloObj> &stream) {
//	l1ct::HadCaloObj clear;
//	while(!stream.empty()){
//		clear = stream.read();
//	}
//}

void getStream(const ap_uint<576> &link,
				const int &phi_offset,
				const l1ct::PFRegion &region,
				l1ct::HadCaloObj outstream[4]
				){
	ap_uint<576> word = link;
	int count = 0;

	for(loop j=0; j<N_PF_LINK; j++) {
		ap_uint<8>  ETA = (word.range(16,12) - 12);
		ap_int<9>   PHI = (word.range(24,17) - phi_offset);

		bool isInside = region.isInside(ETA, PHI);
		if(isInside){
			outstream[count].hwPt = word.range(11,0);
			outstream[count].hwEta = ETA;
			outstream[count].hwPhi = PHI;

			count = count + 1;
		}
		if (count == 4) break;
		word=word>>64;
	}
}




void algo_topIP1(
    ap_uint<576> link_in[N_INPUT_LINKS],
    ap_uint<576> link_out[N_OUTPUT_LINKS]
    ) {
#pragma HLS PIPELINE II=8

#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS INTERFACE ap_ctrl_hs port=return

//	ap_uint<576> link_in2[N_INPUT_LINKS];
//#pragma HLS ARRAY_PARTITION variable=link_in2 complete dim=0
//	for(int i=0 ; i < N_INPUT_LINKS ; i++) {
//		link_in2[i] = link_in[i];
//	}

    // define the  sector boundaries and overlaps
    l1ct::PFRegion region[N_SECTORS];
regions_init:
    for(int i=0 ; i < N_SECTORS ; i++) {
        region[i].hwEtaCenter = l1ct::glbeta_t(12);
        region[i].hwEtaHalfWidth = l1ct::eta_t(12);
        region[i].hwEtaExtra = l1ct::eta_t(0);
        region[i].hwPhiExtra = l1ct::phi_t(2);
        region[i].hwPhiHalfWidth = l1ct::phi_t(6);
        region[i].hwPhiCenter = l1ct::glbphi_t(12*i+5);
    }

    for(int region_n=0 ; region_n < N_SECTORS ; region_n++) {
    	int region_l = (region_n==0) 				? (N_INPUT_LINKS-1):(region_n-1);
    	int region_r = (region_n==N_INPUT_LINKS-1) 	? (0):(region_n+1);

		compute(link_in[region_n], link_in[region_l], link_in[region_r], region[region_n], link_out[region_n]);
		#ifndef __SYNTHESIS__
				std::cout<< region_l << " | " << region_n << " | " << region_r << "\n";
		#endif

    }

//    compute(link_in[0], link_in[5], link_in[1], region[0], link_out[0]);
//    compute(link_in2[1], link_in2[0], link_in2[2], region[1], link_out[1]);
//    compute(link_in[2], link_in[1], link_in[3], region[2], link_out[2]);
//    compute(link_in2[3], link_in2[2], link_in2[4], region[3], link_out[3]);
//    compute(link_in[4], link_in[3], link_in[5], region[4], link_out[4]);
//    compute(link_in2[5], link_in2[4], link_in2[0], region[5], link_out[5]);

}


