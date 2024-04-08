#include "algo_topIP1.h"
#include "axisWrapper.h"
#include <cstdlib>
#include <fstream>
#include <stdlib.h>
#include <linpuppi.h>

#include "hls_stream.h"  
#include "ap_axi_sdata.h"

#define N_TEST_FILES 1
#define N_WORDS 9
#define BIT_WIDTH 64

typedef ap_axis<BIT_WIDTH,0,0,1> axi_stream;

int main() {

    hls::stream<axi_stream> inputStream0;
    hls::stream<axi_stream> inputStream1;
    hls::stream<axi_stream> inputStream2;
    hls::stream<axi_stream> inputStream3;
    hls::stream<axi_stream> inputStream4;
    hls::stream<axi_stream> inputStream5;

    hls::stream<axi_stream> outputStream0;
    hls::stream<axi_stream> outputStream1;
    hls::stream<axi_stream> outputStream2;
    hls::stream<axi_stream> outputStream3;
    hls::stream<axi_stream> outputStream4;
    hls::stream<axi_stream> outputStream5;

    axi_stream buffer;

    srand((unsigned)time(0));

    loop start, end ;

    ap_uint<64> inputLink[6][9];
    ap_uint<64> outputLink[6][9];

    PFcluster InPFcluster[N_PF];

    for (int fileNum = 0; fileNum < N_TEST_FILES; fileNum++) {
        string filename = "../../../../../data/evt_" + to_string(fileNum) + "_clipped.csv";
        string save = "../../../../../outputs/evt_" + to_string(fileNum) + "_puppi.csv";

        std::fstream read (filename);
        if (!read.is_open()) {
            std::cout << "error: file open failed | probably you don't have input files ಠ_ಠ" << "\n";
            return 1;
        }

        string data;
        int lineNum = 0;
        int num = 0;
        float temp;
        while (getline(read, data)) {
            if (lineNum == 0) {
                lineNum++;
                continue;
            }
            if (lineNum > 49) {
                break;
            }

            stringstream ss(data);


            while (ss.good()) {
                string substr;
                getline(ss, substr, ',');
                temp = std::atof(substr.data());
                if (num == 0) {
                    InPFcluster[lineNum-1].ET = temp;
                    //				std::cout<< InPFcluster[lineNum-1].ET << "\n";
                    num++;
                }
                else if (num == 1) {
                    InPFcluster[lineNum-1].Eta = temp;
                    //               	std::cout<< Eta << "\n";
                    num++;
                }
                else if (num == 2) {
                    InPFcluster[lineNum-1].Phi = temp;
                    //               	std::cout<< Phi << "\n";
                    num = 0;
                }
            }
            lineNum++;
        }

        read.close();

        for(loop i=0; i<N_INPUT_LINKS; i++) {
            start = 0;
            end = start+63;
            for(loop k=0; k<N_PF_LINK; k++) {
                std::cout << "Et : " << InPFcluster[i*N_PF_LINK+k].ET << " Eta : " << InPFcluster[i*N_PF_LINK+k].Eta << " Phi : " << InPFcluster[i*N_PF_LINK+k].Phi << endl;
                inputLink[i][k] = InPFcluster[i*N_PF_LINK+k].data()  ;
                start = start + 64 ;
                end = start + 63 ;
            }
        }

        for(loop k=0; k<N_PF_LINK; k++) {
        	buffer.data = inputLink[0][k]; buffer.last = 0; inputStream0.write(buffer);
        }
        buffer.data = 0; buffer.last = 1; inputStream0.write(buffer);
        for(loop k=0; k<N_PF_LINK; k++) {
        	buffer.data = inputLink[1][k]; buffer.last = 0; inputStream1.write(buffer);
        }
        buffer.data = 0; buffer.last = 1; inputStream1.write(buffer);
        for(loop k=0; k<N_PF_LINK; k++) {
        	buffer.data = inputLink[2][k]; buffer.last = 0; inputStream2.write(buffer);
        }
        buffer.data = 0; buffer.last = 1; inputStream2.write(buffer);
        for(loop k=0; k<N_PF_LINK; k++) {
        	buffer.data = inputLink[3][k]; buffer.last = 0; inputStream3.write(buffer);
        }
        buffer.data = 0; buffer.last = 1; inputStream3.write(buffer);
        for(loop k=0; k<N_PF_LINK; k++) {
        	buffer.data = inputLink[4][k]; buffer.last = 0; inputStream4.write(buffer);
        }
        buffer.data = 0; buffer.last = 1; inputStream4.write(buffer);
        for(loop k=0; k<N_PF_LINK; k++) {
        	buffer.data = inputLink[5][k]; buffer.last = 0; inputStream5.write(buffer);
        }
        buffer.data = 0; buffer.last = 1; inputStream5.write(buffer);

//        algo_topIP1(link_in, link_out) ;
        AXIStream_wrapper(
            inputStream0,
            inputStream1,
            inputStream2,
            inputStream3,
            inputStream4,
            inputStream5,

            outputStream0,
            outputStream1,
            outputStream2,
            outputStream3,
            outputStream4,
            outputStream5
            );

        std::ofstream write (save);
        write << "#et,eta,phi\n";

        const ap_uint<8> BW = 64;

        loop i = 0;
		start = 0;
		end = start+BW ;
		for(loop k=0; k<N_PUPPI_LINK; k++) {
			buffer = outputStream0.read();
			l1ct::PuppiObj puppiOutObj = l1ct::PuppiObj::unpack(buffer.data);
			start = start + BW ;
			end = start + BW ;
			write << puppiOutObj.hwPt << "," << puppiOutObj.hwEta<< "," << puppiOutObj.hwPhi << "\n";
			std::cout << "Link : " << i << " cluster : " << k << " IN : " << InPFcluster[i*N_PF_LINK+k].ET << ","
					  << InPFcluster[i*N_PF_LINK+k].Eta << "," << InPFcluster[i*N_PF_LINK+k].Phi << ",  "<< " OUT : " << puppiOutObj.hwPt << ","
					  << puppiOutObj.hwEta<< "," << puppiOutObj.hwPhi << "\n";
		}

		i = 1;
		start = 0;
		end = start+BW ;
		for(loop k=0; k<N_PUPPI_LINK; k++) {
			buffer = outputStream1.read();
			l1ct::PuppiObj puppiOutObj = l1ct::PuppiObj::unpack(buffer.data);
			start = start + BW ;
			end = start + BW ;
			write << puppiOutObj.hwPt << "," << puppiOutObj.hwEta<< "," << puppiOutObj.hwPhi << "\n";
			std::cout << "Link : " << i << " cluster : " << k << " IN : " << InPFcluster[i*N_PF_LINK+k].ET << ","
					  << InPFcluster[i*N_PF_LINK+k].Eta << "," << InPFcluster[i*N_PF_LINK+k].Phi << ",  "<< " OUT : " << puppiOutObj.hwPt << ","
					  << puppiOutObj.hwEta<< "," << puppiOutObj.hwPhi << "\n";
		}

		i = 2;
		start = 0;
		end = start+BW ;
		for(loop k=0; k<N_PUPPI_LINK; k++) {
			buffer = outputStream2.read();
			l1ct::PuppiObj puppiOutObj = l1ct::PuppiObj::unpack(buffer.data);
			start = start + BW ;
			end = start + BW ;
			write << puppiOutObj.hwPt << "," << puppiOutObj.hwEta<< "," << puppiOutObj.hwPhi << "\n";
			std::cout << "Link : " << i << " cluster : " << k << " IN : " << InPFcluster[i*N_PF_LINK+k].ET << ","
					  << InPFcluster[i*N_PF_LINK+k].Eta << "," << InPFcluster[i*N_PF_LINK+k].Phi << ",  "<< " OUT : " << puppiOutObj.hwPt << ","
					  << puppiOutObj.hwEta<< "," << puppiOutObj.hwPhi << "\n";
		}

		i = 3;
		start = 0;
		end = start+BW ;
		for(loop k=0; k<N_PUPPI_LINK; k++) {
			buffer = outputStream3.read();
			l1ct::PuppiObj puppiOutObj = l1ct::PuppiObj::unpack(buffer.data);
			start = start + BW ;
			end = start + BW ;
			write << puppiOutObj.hwPt << "," << puppiOutObj.hwEta<< "," << puppiOutObj.hwPhi << "\n";
			std::cout << "Link : " << i << " cluster : " << k << " IN : " << InPFcluster[i*N_PF_LINK+k].ET << ","
					  << InPFcluster[i*N_PF_LINK+k].Eta << "," << InPFcluster[i*N_PF_LINK+k].Phi << ",  "<< " OUT : " << puppiOutObj.hwPt << ","
					  << puppiOutObj.hwEta<< "," << puppiOutObj.hwPhi << "\n";
		}

		i = 4;
		start = 0;
		end = start+BW ;
		for(loop k=0; k<N_PUPPI_LINK; k++) {
			buffer = outputStream4.read();
			l1ct::PuppiObj puppiOutObj = l1ct::PuppiObj::unpack(buffer.data);
			start = start + BW ;
			end = start + BW ;
			write << puppiOutObj.hwPt << "," << puppiOutObj.hwEta<< "," << puppiOutObj.hwPhi << "\n";
			std::cout << "Link : " << i << " cluster : " << k << " IN : " << InPFcluster[i*N_PF_LINK+k].ET << ","
					  << InPFcluster[i*N_PF_LINK+k].Eta << "," << InPFcluster[i*N_PF_LINK+k].Phi << ",  "<< " OUT : " << puppiOutObj.hwPt << ","
					  << puppiOutObj.hwEta<< "," << puppiOutObj.hwPhi << "\n";
		}

		i = 5;
		start = 0;
		end = start+BW ;
		for(loop k=0; k<N_PUPPI_LINK; k++) {
			buffer = outputStream5.read();
			l1ct::PuppiObj puppiOutObj = l1ct::PuppiObj::unpack(buffer.data);
			start = start + BW ;
			end = start + BW ;
			write << puppiOutObj.hwPt << "," << puppiOutObj.hwEta<< "," << puppiOutObj.hwPhi << "\n";
			std::cout << "Link : " << i << " cluster : " << k << " IN : " << InPFcluster[i*N_PF_LINK+k].ET << ","
					  << InPFcluster[i*N_PF_LINK+k].Eta << "," << InPFcluster[i*N_PF_LINK+k].Phi << ",  "<< " OUT : " << puppiOutObj.hwPt << ","
					  << puppiOutObj.hwEta<< "," << puppiOutObj.hwPhi << "\n";
		}


        cout << "writing to evt_" + to_string(fileNum) + "_puppi.csv" << endl;
        write.close();
    }
    return 0;
}
