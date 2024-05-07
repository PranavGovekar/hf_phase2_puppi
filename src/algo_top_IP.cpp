#include "algo_topIP1.h"


void processInputLinks(ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], hftower HFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI]){
	#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
  #pragma HLS ARRAY_PARTITION variable=HFTowers complete dim=0

for(loop j=0; j<TOWERS_IN_PHI; j=j+2){
    #pragma HLS UNROLL
    for(loop i=0; i<TOWERS_IN_ETA-1; i++){
    #pragma HLS UNROLL
      ap_uint<10> start   = i*10;
      ap_uint<10> end = start + 9;
      HFTowers[i][j] = hftower(link_in[j/2].range(end, start));
  }
    for(loop i=0; i<TOWERS_IN_ETA-1; i++){
    #pragma HLS UNROLL
      ap_uint<10> start   = i*10+110;
      ap_uint<10> end = start + 9;
      HFTowers[i][j+1] = hftower(link_in[j/2].range(end, start));
  }
}

    for(loop j=0; j<TOWERS_IN_PHI; j=j+2){
    #pragma HLS UNROLL
    hftower A10 = HFTowers[TOWERS_IN_ETA-2][j] ;
    hftower B10 = HFTowers[TOWERS_IN_ETA-2][j+1] ;


    ap_uint<8> halfA = A10.energy >> 1 ;
    ap_uint<8> halfB = B10.energy >> 1 ;


    A10.energy = halfA ;
    B10.energy = halfB ;


    HFTowers[TOWERS_IN_ETA-2][j] = A10;
    HFTowers[TOWERS_IN_ETA-2][j+1] = A10;

    HFTowers[TOWERS_IN_ETA-1][j] = B10;
    HFTowers[TOWERS_IN_ETA-1][j+1] = B10;
  }
}

void algo_topIP1(ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], ap_uint<LINK_WIDTH> link_out[N_OUTPUT_LINKS]){
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS PIPELINE II=9
#pragma HLS INTERFACE ap_ctrl_hs port=return

        hftower HFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI] ;
        #pragma HLS ARRAY_PARTITION variable=HFTowers complete dim=0

	processInputLinks(link_in, HFTowers) ;


        ap_uint<10> start ;
        ap_uint<10> end ;

        for(loop j=0; j<TOWERS_IN_PHI; j++)
          {
          	for(loop i=0; i<TOWERS_IN_ETA; i++)
		{
		start=i*10 ; end=start+9;
		link_out[j].range(end, start) = HFTowers[i][j].gettower() ;
	  }}
}

