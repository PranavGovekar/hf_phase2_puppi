#include "algo_topIP1.h"
#include <cstdlib>
#include <fstream>
#include <stdlib.h>
#include <linpuppi.h>

#define N_TEST_FILES 1

int main() {
    srand((unsigned)time(0));

    loop start, end ;

    ap_uint<576> link_in[N_INPUT_LINKS] ;
    ap_uint<576> link_out[N_OUTPUT_LINKS] ;

    PFcluster InPFcluster[N_PF];
    string i_pattern = "../../../../../outputs/inputPattern.mem";
    string o_pattern = "../../../../../outputs/outputPattern.mem";

    std::ofstream outFile (o_pattern);
    std::ofstream inFile (i_pattern);

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
                link_in[i].range(end, start) = InPFcluster[i*N_PF_LINK+k].data()  ;
                start = start + 64 ;
                end = start + 63 ;
            }
			for(loop p=0; p<576; p++) {
				inFile << link_in[i][575 - p];
//				std::cout << link_in[i][575 - p];
			}
			inFile << "\n";
//			std::cout << "\n";
//			std::cout << link_in[i];
//			std::cout << "\n";
        }

        algo_topIP1(link_in, link_out) ;

        std::ofstream write (save);
        write << "#et,eta,phi\n";



        const ap_uint<8> BW = 64;
        for(loop i=0; i<N_OUTPUT_LINKS; i++) {
            start = 0;
            end = start+BW ;
            for(loop k=0; k<N_PUPPI_LINK; k++) {
                l1ct::PuppiObj puppiOutObj = l1ct::PuppiObj::unpack( link_out[i].range(end, start)  );
                start = start + BW ;
                end = start + BW ;
                write << puppiOutObj.hwPt << "," << puppiOutObj.hwEta<< "," << puppiOutObj.hwPhi << "\n";
                std::cout << "Link : " << i << " cluster : " << k << " IN : " << InPFcluster[i*N_PF_LINK+k].ET << ","
                          << InPFcluster[i*N_PF_LINK+k].Eta << "," << InPFcluster[i*N_PF_LINK+k].Phi << ",  "<< " OUT : " << puppiOutObj.hwPt << ","
                          << puppiOutObj.hwEta<< "," << puppiOutObj.hwPhi << "\n";
            }
			for(loop p=0; p<576; p++) {
				outFile << link_out[i][575 - p];
//				std::cout << link_out[i][575 - p];
			}
			outFile << "\n";
//			std::cout << "\n";
//			std::cout << link_out[i];
//			std::cout << "\n";
        }

        cout << "writing to evt_" + to_string(fileNum) + "_puppi.csv" << endl;
        write.close();
    }
    outFile.close();
    inFile.close();

    return 0;
}
