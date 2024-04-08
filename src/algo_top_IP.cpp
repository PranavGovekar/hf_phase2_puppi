#ifndef ALGO_TOPIP1_CPP
#define ALGO_TOPIP1_CPP

//#define DEBUGME

#ifdef DEBUGME
#include <bitset>
#endif

#include "linpuppi.h"
#include "algo_topIP1.h"
#include <cstdint>

void compute(const ap_uint<576> link_center,
			const ap_uint<576> link_left,
			const ap_uint<576> link_right,
			const ap_uint<3> sector,
			ap_uint<576> &link_out
				){
	#pragma HLS PIPELINE

	l1ct::PFRegion region;
    region.hwEtaCenter = l1ct::glbeta_t(12);
    region.hwEtaHalfWidth = l1ct::eta_t(12);
    region.hwEtaExtra = l1ct::eta_t(0);
    region.hwPhiExtra = l1ct::phi_t(2);
    region.hwPhiHalfWidth = l1ct::phi_t(6);
    region.hwPhiCenter = l1ct::glbphi_t((12*sector)+5);

	l1ct::HadCaloObj puppiIn[NCALO];
//#pragma HLS ARRAY_PARTITION dim=1 type=complete variable=puppiIn
	l1ct::PuppiObj pfselne[NNEUTRALS];
//#pragma HLS ARRAY_PARTITION dim=1 type=complete variable=pfselne

	fillCenterLink(link_center, region, puppiIn);
	fillExtra(link_left, link_right, region, sector, puppiIn);

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
	ap_uint<576> temp_link;

    for(loop idx=0; idx<NNEUTRALS; idx++) {
    	temp_link.range(start+(BW-1),start) = pfselne[idx].pack();

        start = start + BW;
    }

    link_out = temp_link;
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
				const ap_uint<3> N_REGION,
				l1ct::HadCaloObj puppiIn[NCALO]
				){
	ap_uint<3> region_l;
	ap_uint<3> region_r;
	if (N_REGION==0){
		region_l = 5;
	}
	else {
		region_l = N_REGION-1;
	}
	if (N_REGION==5){
		region_r = 0;
	}
	else {
		region_r = N_REGION+1;
	}

	int phi_offset_l = 12*N_REGION + 5;
	int phi_offset_r = 12*N_REGION + 5;

	if (  N_REGION == 0              and    region_l == (5)) phi_offset_l+=72;
	if (  N_REGION == (5)  and    region_r == 0                ) phi_offset_r-=72;

	l1ct::HadCaloObj leftStream[N_EXTRA];
	l1ct::HadCaloObj rightStream[N_EXTRA];

	for(loop idx=0; idx<N_EXTRA; idx++) {
		leftStream[idx].clear();
		rightStream[idx].clear();
	}

	getInside(link_left, phi_offset_l, region, leftStream);
	getInside(link_right, phi_offset_r, region, rightStream);

	mergeSort(leftStream, rightStream, puppiIn);

}

void mergeSort(	l1ct::HadCaloObj leftStream[N_EXTRA],
				l1ct::HadCaloObj rightStream[N_EXTRA],
				l1ct::HadCaloObj puppiIn[NCALO]
	){

	l1ct::HadCaloObj left;
	l1ct::HadCaloObj right;

	int idx_left = 0;
	int idx_right = 0;

	left = leftStream[idx_left];
	right = rightStream[idx_right];

	for(loop idx=0; idx<N_EXTRA; idx++) {
		if(left.hwPt < right.hwPt){
			puppiIn[NNEUTRALS+idx] = right;
			idx_right = idx_right + 1;
			right = rightStream[idx_right];
		}
		else {
			puppiIn[NNEUTRALS+idx] = left;
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

void getInside(const ap_uint<576> &link,
				const int &phi_offset,
				const l1ct::PFRegion &region,
				l1ct::HadCaloObj outstream[N_EXTRA]
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
		if (count == N_EXTRA) break;
		word=word>>64;
	}
}




void algo_topIP1(
    ap_uint<576> link_in[N_INPUT_LINKS],
    ap_uint<576> link_out[N_OUTPUT_LINKS]
    ) {
#pragma HLS PIPELINE
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS INTERFACE ap_ctrl_hs port=return

	ap_uint<576> link_in_1[N_INPUT_LINKS];
#pragma HLS ARRAY_PARTITION variable=link_in_1 complete dim=0
	ap_uint<576> link_in_2[N_INPUT_LINKS];
#pragma HLS ARRAY_PARTITION variable=link_in_2 complete dim=0
	for(int idx=0; idx < N_INPUT_LINKS; idx++) {
#pragma HLS UNROLL
		link_in_1[idx] = link_in[idx];
		link_in_2[idx] = link_in[idx];
	}

	ap_uint<3> left;
	ap_uint<3> right;

    for(int idx=0; idx < N_SECTORS; idx++) {

    	if (idx==0){
    		left = 5;
    	}
    	else {
    		left = idx-1;
    	}
    	if (idx==5){
    		right = 0;
    	}
    	else {
    		right = idx+1;
    	}

		compute(link_in[idx], link_in_1[left], link_in_2[right], idx, link_out[idx]);
	}
}

void ReadWrite( ap_uint<64> in[54], 
                ap_uint<64> out[54]
            ){


void stream_in(
    hls::stream<axi_stream> &inputStream,
    ap_uint<576> *link_in
    ) { 
    
    axi_stream buffer;
    ap_uint<576> tempLink;
    ap_uint<BIT_WIDTH> tempArray[N_WORDS];
    ap_uint<12> start = 0;

    for(loop i=0; i<N_WORDS; i++) {
        buffer = inputStream.read();
        tempLink.range(start+(BIT_WIDTH-1),start) = buffer.data;

        #ifdef DEBUGME

        std::cout << "in_stream buff :" << buffer.data << endl;
        std::cout << "in_stream buff :" << tempLink.range(start+(BIT_WIDTH-1),start) << endl;
        std::cout << "tempLink inside stream_in : " << std::bitset<64>(tempLink.range(start+(BIT_WIDTH-1),start)) << endl;
        std::cout << "start+B , start :" << start+(BIT_WIDTH-1) << "," << start << endl;
        #endif
        
        start = start + BIT_WIDTH;

        if(i == N_WORDS-1){
            start = 0;
        }
    }

    #ifdef DEBUGME
    std::cout << "tempLink inside stream_in : " << std::bitset<576>(tempLink) << endl;
    #endif

    *link_in = tempLink;
}

void stream_out(
    hls::stream<axi_stream> &outputStream,
    ap_uint<576> *link_out
    ) {

    axi_stream buffer;
    ap_uint<576> tempLink;
    ap_uint<12> start = 0;

    tempLink = *link_out;

    #ifdef DEBUGME
    std::cout << "*link_out : " << std::bitset<576>(*link_out) << endl;
    std::cout << "tempLink : " << std::bitset<576>(tempLink) << endl;
    #endif

    for(loop i=0; i<N_WORDS; i++) {
        #ifdef DEBUGME
        std::cout << "start+B , start :" << start+(BIT_WIDTH-1) << "," << start << endl;
        std::cout << "tempLink ranged: " << std::bitset<64>(tempLink.range(start+(BIT_WIDTH-1),start)) << endl;
        #endif
        
        buffer.data = tempLink.range(start+(BIT_WIDTH-1),start);
        start = start + BIT_WIDTH;

        #ifdef DEBUGME
        std::cout << "out_stream :" << buffer.data << endl;
        std::cout << "i insdie the stream_out for loop :" << i << endl;
        #endif

        if (i == N_WORDS-1){
            #ifdef DEBUGME
            cout << "i insdie the stream_out for loop :" << i << endl;            
            #endif
            
            buffer.last = 1;
            start = 0;
        }
        else {
            buffer.last = 0;
        }

        outputStream.write(buffer);
    } 
}

void AXIStream_wrapper(
    hls::stream<axi_stream> &inputStream0,
    hls::stream<axi_stream> &inputStream1,
    hls::stream<axi_stream> &inputStream2,
    hls::stream<axi_stream> &inputStream3,
    hls::stream<axi_stream> &inputStream4,
    hls::stream<axi_stream> &inputStream5,
    // if it works, it works.
    hls::stream<axi_stream> &outputStream0,
    hls::stream<axi_stream> &outputStream1,
    hls::stream<axi_stream> &outputStream2,
    hls::stream<axi_stream> &outputStream3,
    hls::stream<axi_stream> &outputStream4,
    hls::stream<axi_stream> &outputStream5
    ) {

#pragma HLS INTERFACE mode=ap_ctrl_hs port=return
#pragma HLS DATAFLOW
#pragma HLS INTERFACE mode=axis register_mode=both port=inputStream0 register
#pragma HLS INTERFACE mode=axis register_mode=both port=inputStream1 register
#pragma HLS INTERFACE mode=axis register_mode=both port=inputStream2 register
#pragma HLS INTERFACE mode=axis register_mode=both port=inputStream3 register
#pragma HLS INTERFACE mode=axis register_mode=both port=inputStream4 register
#pragma HLS INTERFACE mode=axis register_mode=both port=inputStream5 register
#pragma HLS INTERFACE mode=axis register_mode=both port=outputStream0 register
#pragma HLS INTERFACE mode=axis register_mode=both port=outputStream1 register
#pragma HLS INTERFACE mode=axis register_mode=both port=outputStream2 register
#pragma HLS INTERFACE mode=axis register_mode=both port=outputStream3 register
#pragma HLS INTERFACE mode=axis register_mode=both port=outputStream4 register
#pragma HLS INTERFACE mode=axis register_mode=both port=outputStream5 register
    
    ap_uint<576> link_in[N_INPUT_LINKS];
    ap_uint<576> link_out[N_OUTPUT_LINKS];

    stream_in(inputStream0, &link_in[0]);
    stream_in(inputStream1, &link_in[1]);
    stream_in(inputStream2, &link_in[2]);
    stream_in(inputStream3, &link_in[3]);
    stream_in(inputStream4, &link_in[4]);
    stream_in(inputStream5, &link_in[5]);

     algo_topIP1(link_in, link_out);

//   for(loop i=0; i<6; i++) {
//   	link_out[i] = link_in[i];
//   }

    stream_out(outputStream0, &link_out[0]);
    stream_out(outputStream1, &link_out[1]);
    stream_out(outputStream2, &link_out[2]);
    stream_out(outputStream3, &link_out[3]);
    stream_out(outputStream4, &link_out[4]);
    stream_out(outputStream5, &link_out[5]);

}

#endif
