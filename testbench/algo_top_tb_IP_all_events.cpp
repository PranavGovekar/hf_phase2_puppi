#include "algo_topIP1.h"
#include <cstdlib>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <bitset>

int main(){
	srand((unsigned)time(0));

	loop start, end ;

	ap_uint<576> link_in[N_INPUT_LINKS] ;
	ap_uint<576> link_out[N_OUTPUT_LINKS] ;

    PFcluster InPFcluster[N_PF];

    for (int fileNum = 0; fileNum < 10; fileNum++){
    	std::cout<<"## EVENT_BEG "<<fileNum<<"\n";
    	string filename = "../../../../../inputs/evt_" + to_string(fileNum) + "_clipped.csv";
    	string save = "../../../../../outputs/evt_" + to_string(fileNum) + "_puppi.csv";

		std::fstream read (filename);
		if (!read.is_open()) {
			std::cout << "error: file open failed " << ".\n";
			return 1;
		}

		string data;
		int lineNum = 0;
		int num = 0;
		float temp,eta,eta3;
		uint8_t eta_discrete, eta_discrete_masked ;
		while (getline(read, data)) {
			if (lineNum == 0){
				lineNum++;
				continue;
			}
			if (lineNum > 49){
				break;
			}

			stringstream ss(data);

			while (ss.good()) {
				string substr;
				getline(ss, substr, ',');
				temp = std::atof(substr.data());
				if (num == 0){
					InPFcluster[lineNum-1].ET = temp;
					num++;
				}
				else if (num == 1){
        			eta_discrete =  abs(temp) - 30;
					if (temp < 0)
					{
           				eta_discrete_masked = 1<<4 | eta_discrete;
					
					}
					else
						{
							eta_discrete_masked = 	eta_discrete ;
						}
					InPFcluster[lineNum-1].Eta= eta_discrete_masked;

					num++;
					}
			   
				else if (num == 2){
					InPFcluster[lineNum-1].Phi = temp;
					num = 0;
			   }
			}
			lineNum++;
		}

		read.close();

		for(loop i=0; i<N_INPUT_LINKS; i++){
		start = 0; end = start+63 ;
				for(loop k=0; k<N_PF_LINK; k++){
				  link_in[i].range(end, start) = InPFcluster[i*N_PF_LINK+k].data()  ;
				  start = start + 64 ; end = start + 63 ;
		}}

		algo_topIP1(link_in, link_out) ;

		std::ofstream write (save);
		write << "#et,eta,phi\n";

			for(loop i=0; i<N_OUTPUT_LINKS; i++){
			start = 0; end = start+63 ;
					for(loop k=0; k<N_PF_LINK; k++){
					  InPFcluster[i*N_PF_LINK+k].getPFcluster(link_out[i].range(end, start)) ;
					  start = start + 64 ; end = start + 63 ;
					  write << InPFcluster[i*8+k].ET << "," << InPFcluster[i*8+k].Eta << "," << InPFcluster[i*8+k].Phi << "\n";			  
					}
			}
		 cout << "writing to evt_" + to_string(fileNum) + "_puppi.csv" << endl;
		 write.close();
		 std::cout<<"## EVENT_END "<<fileNum<<"\n";
    }
  return 0;
}
