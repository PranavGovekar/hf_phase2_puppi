#include "algo_topIP1.h"
#include <cstdlib>

int main(){

srand((unsigned)time(0));

loop start, end ;
ap_uint<576> link_in[N_INPUT_LINKS] ;
ap_uint<576> link_out[N_OUTPUT_LINKS] ;
    PFcluster InPFcluster[N_PF];

        for(loop i=0; i<N_INPUT_LINKS; i++){
                link_in[i] = 0 ;
        }
        for(loop i=0; i<N_PF; i++){
        	if(i==0) {InPFcluster[i].ET = 90.0; InPFcluster[i].Eta = 5; InPFcluster[i].Phi = 58.0; }  // eta = 35.0
        	if(i==1) {InPFcluster[i].ET = 70.0; InPFcluster[i].Eta = 4; InPFcluster[i].Phi = 40.0; }  // eta = 34.0
        	if(i==2) {InPFcluster[i].ET = 68.0; InPFcluster[i].Eta = 1<<4 | 5; InPFcluster[i].Phi = 34.0; }  // eta = -35.0
        	if(i==3) {InPFcluster[i].ET = 65.0; InPFcluster[i].Eta = 1<<4 | 7; InPFcluster[i].Phi = 42.0; }  // eta = -37.0
        	if(i==4) {InPFcluster[i].ET = 62.0; InPFcluster[i].Eta = 5; InPFcluster[i].Phi = 22.0; }  // eta = 35.0
        	if(i==5) {InPFcluster[i].ET = 57.0; InPFcluster[i].Eta = 1<<4 | 6; InPFcluster[i].Phi = 24.0; }  // eta = -36.0
        	if(i==6) {InPFcluster[i].ET = 57.0; InPFcluster[i].Eta = 3; InPFcluster[i].Phi = 54.0; }  // eta = 33.0
        	if(i==7) {InPFcluster[i].ET = 55.0; InPFcluster[i].Eta = 1<<4 | 6; InPFcluster[i].Phi = 2.0; }  // eta = -36.0
        	if(i==8) {InPFcluster[i].ET = 54.0; InPFcluster[i].Eta = 6; InPFcluster[i].Phi = 44.0; }  // eta = 36.0
        	if(i==9) {InPFcluster[i].ET = 53.0; InPFcluster[i].Eta = 1<<4 | 4; InPFcluster[i].Phi = 12.0; }  // eta = -34.0
        	if(i==10) {InPFcluster[i].ET = 52.0; InPFcluster[i].Eta = 1<<4 | 5; InPFcluster[i].Phi = 42.0; }  // eta = -35.0
        	if(i==11) {InPFcluster[i].ET = 51.0; InPFcluster[i].Eta = 1<<4 | 3; InPFcluster[i].Phi = 28.0; }  // eta = -33.0
        	if(i==12) {InPFcluster[i].ET = 50.0; InPFcluster[i].Eta = 3; InPFcluster[i].Phi = 8.0; }  // eta = 33.0
        	if(i==13) {InPFcluster[i].ET = 49.0; InPFcluster[i].Eta = 1<<4 | 4; InPFcluster[i].Phi = 48.0; }  // eta = -34.0
        	if(i==14) {InPFcluster[i].ET = 49.0; InPFcluster[i].Eta = 9; InPFcluster[i].Phi = 4.0; }  // eta = 39.0
        	if(i==15) {InPFcluster[i].ET = 48.0; InPFcluster[i].Eta = 1<<4 | 7; InPFcluster[i].Phi = 18.0; }  // eta = -37.0
        	if(i==16) {InPFcluster[i].ET = 48.0; InPFcluster[i].Eta = 1<<4 | 3; InPFcluster[i].Phi = 72.0; }  // eta = -33.0
        	if(i==17) {InPFcluster[i].ET = 47.0; InPFcluster[i].Eta = 5; InPFcluster[i].Phi = 30.0; }  // eta = 35.0
        	if(i==18) {InPFcluster[i].ET = 45.0; InPFcluster[i].Eta = 3; InPFcluster[i].Phi = 70.0; }  // eta = 33.0
        	if(i==19) {InPFcluster[i].ET = 44.0; InPFcluster[i].Eta = 1<<4 | 11; InPFcluster[i].Phi = 62.0; }  // eta = -41.0
        	if(i==20) {InPFcluster[i].ET = 44.0; InPFcluster[i].Eta = 1<<4 | 9; InPFcluster[i].Phi = 70.0; }  // eta = -39.0
        	if(i==21) {InPFcluster[i].ET = 44.0; InPFcluster[i].Eta = 6; InPFcluster[i].Phi = 10.0; }  // eta = 36.0
        	if(i==22) {InPFcluster[i].ET = 44.0; InPFcluster[i].Eta = 8; InPFcluster[i].Phi = 34.0; }  // eta = 38.0
        	if(i==23) {InPFcluster[i].ET = 43.0; InPFcluster[i].Eta = 8; InPFcluster[i].Phi = 66.0; }  // eta = 38.0
        	if(i==24) {InPFcluster[i].ET = 41.0; InPFcluster[i].Eta = 1<<4 | 4; InPFcluster[i].Phi = 18.0; }  // eta = -34.0
        	if(i==25) {InPFcluster[i].ET = 41.0; InPFcluster[i].Eta = 1<<4 | 3; InPFcluster[i].Phi = 36.0; }  // eta = -33.0
        	if(i==26) {InPFcluster[i].ET = 41.0; InPFcluster[i].Eta = 5; InPFcluster[i].Phi = 72.0; }  // eta = 35.0
        	if(i==27) {InPFcluster[i].ET = 40.0; InPFcluster[i].Eta = 1<<4 | 5; InPFcluster[i].Phi = 64.0; }  // eta = -35.0
        	if(i==28) {InPFcluster[i].ET = 40.0; InPFcluster[i].Eta = 1<<4 | 2; InPFcluster[i].Phi = 8.0; }  // eta = -32.0
        	if(i==29) {InPFcluster[i].ET = 39.0; InPFcluster[i].Eta = 1<<4 | 9; InPFcluster[i].Phi = 52.0; }  // eta = -39.0
        	if(i==30) {InPFcluster[i].ET = 39.0; InPFcluster[i].Eta = 1<<4 | 8; InPFcluster[i].Phi = 30.0; }  // eta = -38.0
        	if(i==31) {InPFcluster[i].ET = 38.0; InPFcluster[i].Eta = 8; InPFcluster[i].Phi = 56.0; }  // eta = 38.0
        	if(i==32) {InPFcluster[i].ET = 38.0; InPFcluster[i].Eta = 9; InPFcluster[i].Phi = 12.0; }  // eta = 39.0
        	if(i==33) {InPFcluster[i].ET = 37.0; InPFcluster[i].Eta = 1<<4 | 11; InPFcluster[i].Phi = 2.0; }  // eta = -41.0
        	if(i==34) {InPFcluster[i].ET = 37.0; InPFcluster[i].Eta = 1<<4 | 11; InPFcluster[i].Phi = 42.0; }  // eta = -41.0
        	if(i==35) {InPFcluster[i].ET = 37.0; InPFcluster[i].Eta = 1<<4 | 9; InPFcluster[i].Phi = 14.0; }  // eta = -39.0
        	if(i==36) {InPFcluster[i].ET = 34.0; InPFcluster[i].Eta = 1<<4 | 11; InPFcluster[i].Phi = 30.0; }  // eta = -41.0
        	if(i==37) {InPFcluster[i].ET = 28.0; InPFcluster[i].Eta = 1<<4 | 5; InPFcluster[i].Phi = 58.0; }  // eta = -35.0
        	if(i==38) {InPFcluster[i].ET = 23.0; InPFcluster[i].Eta = 5; InPFcluster[i].Phi = 57.0; }  // eta = 35.0
        	if(i==39) {InPFcluster[i].ET = 13.0; InPFcluster[i].Eta = 1<<4 | 7; InPFcluster[i].Phi = 41.0; }  // eta = -37.0
        	if(i==40) {InPFcluster[i].ET = 12.0; InPFcluster[i].Eta = 1<<4 | 5; InPFcluster[i].Phi = 33.0; }  // eta = -35.0
        	if(i==41) {InPFcluster[i].ET = 12.0; InPFcluster[i].Eta = 1<<4 | 3; InPFcluster[i].Phi = 27.0; }  // eta = -33.0
        	if(i==42) {InPFcluster[i].ET = 12.0; InPFcluster[i].Eta = 4; InPFcluster[i].Phi = 39.0; }  // eta = 34.0
        	if(i==43) {InPFcluster[i].ET = 12.0; InPFcluster[i].Eta = 4; InPFcluster[i].Phi = 59.0; }  // eta = 34.0
        	if(i==44) {InPFcluster[i].ET = 12.0; InPFcluster[i].Eta = 4; InPFcluster[i].Phi = 60.0; }  // eta = 34.0
        	if(i==45) {InPFcluster[i].ET = 11.0; InPFcluster[i].Eta = 1<<4 | 6; InPFcluster[i].Phi = 23.0; }  // eta = -36.0
        	if(i==46) {InPFcluster[i].ET = 10.0; InPFcluster[i].Eta = 1<<4 | 6; InPFcluster[i].Phi = 1.0; }  // eta = -36.0
        	if(i==47) {InPFcluster[i].ET = 10.0; InPFcluster[i].Eta = 1<<4 | 5; InPFcluster[i].Phi = 41.0; }  // eta = -35.0
        	cout << "IN cluster " << i << " ET "<< InPFcluster[i].ET << " Eta " << InPFcluster[i].Eta<< " Phi " << InPFcluster[i].Phi<<"\n";
        }
        for(loop i=0; i<N_INPUT_LINKS; i++){
        start = 0; end = start+63 ;
                for(loop k=0; k<N_PF_LINK; k++){
                  link_in[i].range(end, start) = InPFcluster[i*N_PF_LINK+k].data()  ;
                  start = start + 64 ; end = start + 63 ;
        }
    }

      algo_topIP1(link_in, link_out) ;

        for(loop i=0; i<N_OUTPUT_LINKS; i++){
        start = 0; end = start+63 ;
        	for(loop k=0; k<N_PF_LINK; k++){
        		InPFcluster[i*N_PF_LINK+k].getPFcluster(link_out[i].range(end, start)) ;
                start = start + 64 ; end = start + 63 ;
                cout << "OUT cluster " << i*N_PF_LINK+k << " ET " << InPFcluster[i*N_PF_LINK+k].ET << " Eta " << InPFcluster[i*N_PF_LINK+k].Eta << " Phi " << InPFcluster[i*N_PF_LINK+k].Phi<<"\n";
        }
     }
  return 0;
}
