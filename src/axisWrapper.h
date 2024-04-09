#ifndef AXIS_WRAPPER_H
#define AXIS_WRAPPER_H

#include "algo_topIP1.h"
#include "hls_stream.h"
#include "ap_axi_sdata.h"

#define N_WORDS 9
#define BIT_WIDTH 64
typedef ap_axis <BIT_WIDTH,0,0,1> axi_stream;

void stream_in(
    hls::stream<axi_stream> &inputStream,
    ap_uint<576> *link_in);

void stream_out(
    hls::stream<axi_stream> &outputStream,
    ap_uint<576> *link_out);

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
    hls::stream<axi_stream> &outputStream5);

#endif
