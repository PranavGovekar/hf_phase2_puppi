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
                if(i==0) {InPFcluster[i].ET = 10 ; InPFcluster[i].Eta = 5; InPFcluster[i].Phi = 3; }
                if(i==40) {InPFcluster[i].ET = 110 ; InPFcluster[i].Eta = 7; InPFcluster[i].Phi = 31; }
                if(i==41) {InPFcluster[i].ET = 90  ; InPFcluster[i].Eta = 6; InPFcluster[i].Phi = 28; }
                cout << "IN cluster " << i << " ET "
                	  << InPFcluster[i].ET << " Eta " << InPFcluster[i].Eta
					  << " Phi " << InPFcluster[i].Phi << endl ;
        }

        for(loop i=0; i<N_INPUT_LINKS; i++){
        start = 0; end = start+63 ;
                for(loop k=0; k<N_PF_LINK; k++){
                  link_in[i].range(end, start) = InPFcluster[i*N_PF_LINK+k].data()  ;
                  start = start + 64 ; end = start + 63 ;
        }}

      algo_topIP1(link_in, link_out) ;


        for(loop i=0; i<N_OUTPUT_LINKS; i++){
        start = 0; end = start+63 ;
                for(loop k=0; k<N_PF_LINK; k++){
                  InPFcluster[i*N_PF_LINK+k].getPFcluster(link_out[i].range(end, start)) ;
                  start = start + 64 ; end = start + 63 ;
cout << "OUT cluster " << i*8+k << " ET " << InPFcluster[i*8+k].ET << " Eta " << InPFcluster[i*8+k].Eta << " Phi " << InPFcluster[i*8+k].Phi << endl ;
        }}


  return 0;
}
