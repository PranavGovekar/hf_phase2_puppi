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

	hls::stream<l1ct::HadCaloObj> mainStream;
	hls::stream<l1ct::HadCaloObj> extraStream;

	const int N_REGION = (region.hwPhiCenter - 5)/12;

	fillCenterLink(link_center, region, mainStream);
	fillExtra(link_left, link_right, region, N_REGION, extraStream);
	fill(mainStream, extraStream, puppiIn);

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

void fill(hls::stream<l1ct::HadCaloObj> &mainStream,
		hls::stream<l1ct::HadCaloObj> &extraStream,
		l1ct::HadCaloObj puppiIn[NCALO]){

	l1ct::HadCaloObj tempObj;

	tempObj.hwPt = 0;
	tempObj.hwEta = 0;
	tempObj.hwPhi = 0;

	for(loop idx=0; idx<N_PF_LINK; idx++) {
		puppiIn[idx] = mainStream.read();
	}

	for(loop idx=0; idx<(NCALO-N_PF_LINK); idx++) {
		if(extraStream.empty()) {
			puppiIn[N_PF_LINK+idx] = tempObj;
		}
		else {
			puppiIn[N_PF_LINK+idx] = extraStream.read();
		}
	}
}

void fillCenterLink(const ap_uint<576> &link,
					const l1ct::PFRegion &region,
					hls::stream<l1ct::HadCaloObj> &mainStream
				){
	l1ct::HadCaloObj tempObj;
	ap_uint<576> word = link;
	const int eta_offset = region.hwEtaCenter;
	const int phi_offset = region.hwPhiCenter;

	for(loop j=0; j<N_PF_LINK; j++) {
		tempObj.hwPt = word.range(11,0);
		tempObj.hwEta = (word.range(16,12) - eta_offset);
		tempObj.hwPhi = (word.range(24,17) - phi_offset);
		word=word>>64;

		mainStream.write(tempObj);
	}
}

void fillExtra(	const ap_uint<576> &link_left,
				const ap_uint<576> &link_right,
				const l1ct::PFRegion &region,
				const int N_REGION,
				hls::stream<l1ct::HadCaloObj> &extraStream
				){
#pragma HLS DATAFLOW

	int region_l = (N_REGION==0) 				? (N_INPUT_LINKS-1):(N_REGION-1);
	int region_r = (N_REGION==N_INPUT_LINKS-1) 	? (0):(N_REGION+1);

	int phi_offset_l = 12*region_l + 5;
	int phi_offset_r = 12*region_r + 5;

	if (  N_REGION == 0              and    region_l == (N_INPUT_LINKS-1)) phi_offset_l+=72;
	if (  N_REGION == (N_SECTORS-1)  and    region_r == 0                ) phi_offset_r-=72;

	hls::stream<l1ct::HadCaloObj> leftStream;
	hls::stream<l1ct::HadCaloObj> rightStream;

	getStream(link_left, phi_offset_l, region, leftStream);
	getStream(link_right, phi_offset_r, region, rightStream);

	merge(leftStream, rightStream, extraStream);

}

void merge(	hls::stream<l1ct::HadCaloObj> &leftStream,
			hls::stream<l1ct::HadCaloObj> &rightStream,
			hls::stream<l1ct::HadCaloObj> &extraStream
			){
	l1ct::HadCaloObj left;
	l1ct::HadCaloObj right;
	l1ct::HadCaloObj clear;

	left = leftStream.read();
	right = rightStream.read();

	for(loop idx=0; idx<N_PF_LINK*2; idx++) {
		if (idx < (NCALO - N_PF_LINK)){
			if (left.hwPt == 0 && right.hwPt == 0) {
				break;
			}
			else if(left.hwPt < right.hwPt){
				extraStream.write(right);
				right = rightStream.read();
			}
			else {
				extraStream.write(left);
				left = leftStream.read();
			}
		}
		else {
			if(!rightStream.empty()){
				clear = rightStream.read();
			}
			if(!leftStream.empty()){
				clear = leftStream.read();
			}
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
				hls::stream<l1ct::HadCaloObj> &outstream
				){

	l1ct::HadCaloObj tempObj;
	ap_uint<576> word = link;

	for(loop j=0; j<N_PF_LINK; j++) {
		ap_uint<8>  ETA = (word.range(16,12) - 12);
		ap_int<9>   PHI = (word.range(24,17) - phi_offset);

		bool isInside = region.isInside(ETA, PHI);
		if(isInside){
			tempObj.hwPt = word.range(11,0);
			tempObj.hwEta = ETA;
			tempObj.hwPhi = PHI;
			outstream.write(tempObj);
		}

		word=word>>64;
	}

	tempObj.hwPt = 0;
	tempObj.hwEta = 0;
	tempObj.hwPhi = 0;
	outstream.write(tempObj);
}




void algo_topIP1(
    ap_uint<576> link_in[N_INPUT_LINKS],
    ap_uint<576> link_out[N_OUTPUT_LINKS]
    ) {
#pragma HLS PIPELINE II=8

#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS INTERFACE ap_ctrl_hs port=return

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
    }

}


