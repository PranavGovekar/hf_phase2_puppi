#include "algo_topIP1.h"
#include <cstdlib>
#include <fstream>
#include <stdlib.h>
#include <linpuppi.h>


int main(){
	loop start, end ;
    ap_uint<64> inputLink[] = {0x14a01f,0x0d001d,0x14401c,0x124006,0x16e006,0x172006,0x66005,0x086005,0x0,0x1ce022,0x19201e,0x2e6009,0x18e006,0x1ae006,0x2ee006,0x270006,0x290006,0x0,0x30602a,0x3ce025,0x46600a,0x3ae007,0x472007,0x3e8006,0x408006,0x36a006,0x0,0x548035,0x486030,0x590020,0x49201d,0x528009,0x52a008,0x54a008,0x4a8007,0x0,0x748031,0x650025,0x728009,0x6e6008,0x706008,0x766008,0x628008,0x648008,0x0,0x0886029,0x084e028,0x0812022,0x786008,0x082e008,0x0866007,0x08e6006,0x7ec006,0x0};
    ap_uint<64> outputLink[54];

    ReadWrite(inputLink, outputLink);

    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";
    std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++";

    for(loop i=0; i<54; i++) {
        std::cout << outputLink[i] << endl;
    }

    const ap_uint<8> BW = 64;
    for(loop i=0; i<N_OUTPUT_LINKS; i++) {
        start = 0;
        end = start+BW ;
        for(loop k=0; k<N_PUPPI_LINK; k++) {
            l1ct::PuppiObj puppiOutObj = l1ct::PuppiObj::unpack(outputLink[(i*N_PUPPI_LINK)+k]);
            start = start + BW ;
            end = start + BW ;
            std::cout << "OUT : " << puppiOutObj.hwPt << ","
                      << puppiOutObj.hwEta<< "," << puppiOutObj.hwPhi << "\n";
        }
    }
}
