#include "axisWrapper.h"

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

        #ifdef __SYNTHESIS__

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

    #ifdef __SYNTHESIS__
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

    #ifdef __SYNTHESIS__
    std::cout << "*link_out : " << std::bitset<576>(*link_out) << endl;
    std::cout << "tempLink : " << std::bitset<576>(tempLink) << endl;
    #endif

    for(loop i=0; i<N_WORDS; i++) {
        #ifdef __SYNTHESIS__
        std::cout << "start+B , start :" << start+(BIT_WIDTH-1) << "," << start << endl;
        std::cout << "tempLink ranged: " << std::bitset<64>(tempLink.range(start+(BIT_WIDTH-1),start)) << endl;
        #endif

        buffer.data = tempLink.range(start+(BIT_WIDTH-1),start);
        start = start + BIT_WIDTH;

        #ifdef __SYNTHESIS__
        std::cout << "out_stream :" << buffer.data << endl;
        std::cout << "i insdie the stream_out for loop :" << i << endl;
        #endif

        if (i == N_WORDS-1){
            #ifdef __SYNTHESIS__
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
