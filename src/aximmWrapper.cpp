#include "aximmWrapper.h"

void ReadWrite( ap_uint<64> in[54],
                ap_uint<64> out[54]
              ) {

    #pragma HLS INTERFACE m_axi port=in bundle=aximm1
    #pragma HLS INTERFACE m_axi port=out bundle=aximm2

    const ap_uint<8> BW = 64;
    const ap_uint<8> N_WORDS = 9;

    ap_uint<576> link_in[N_INPUT_LINKS];
    ap_uint<576> link_out[N_OUTPUT_LINKS];
    ap_uint<576> tempLink;
    ap_uint<12> start = 0;

    for(loop i=0; i<N_INPUT_LINKS; i++) {
        for(loop j=0; j<N_WORDS; j++) {
            #pragma HLS PIPELINE
#ifdef __SYNTHESIS__
            std::cout << "\nelement : " << i*N_WORDS+j;
            std::cout << "\n" << start << endl;
            std::cout << in[(i*N_WORDS)+j] <<endl;
            std::cout << std::bitset<64>(in[(i*N_WORDS)+j]) << endl;
#endif
            tempLink.range(start+(BW-1),start) = in[(i*N_WORDS)+j];
            start = start + BW;
        }
#ifdef __SYNTHESIS__
        std::cout << "\ntempLink " << i << " : " <<  tempLink;
        std::cout << "\ntempLink binary " << i << " : " <<  std::bitset<576>(tempLink);
#endif

        link_in[i] = tempLink;
        start = 0;
    }

#ifdef __SYNTHESIS__
    std::cout << "HERE+++++++++++++++++++++++++++" << endl;
    for(loop i=0; i<N_OUTPUT_LINKS; i++) {
        std::cout << std::bitset<576>(link_in[i]) << endl;
    }
    std::cout << "HERE+++++++++++++++++++++++++++" << endl;
#endif

    algo_topIP1(link_in,link_out);

#ifdef __SYNTHESIS__
    std::cout << "HERE+++++++++++++++++++++++++++" << endl;
    for(loop i=0; i<N_OUTPUT_LINKS; i++) {
        std::cout << std::bitset<576>(link_out[i]) << endl;
    }
    std::cout << "HERE+++++++++++++++++++++++++++" << endl;
#endif

    for(loop i=0; i<N_OUTPUT_LINKS; i++) {
        tempLink = link_out[i];
        for(loop j=0; j<N_WORDS; j++) {
            #pragma HLS PIPELINE
            out[(i*N_WORDS)+j] = tempLink.range(start+(BW-1),start);
            start = start + BW;
        }
        start = 0;
    }
}
