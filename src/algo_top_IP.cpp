#include "algo_topIP1.h"
void algo_topIP1(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]){
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS PIPELINE II=9
#pragma HLS INTERFACE ap_ctrl_hs port=return
PFcluster APFcluster[N_PF];
PFcluster PFclusterOutput[N_PF];
#pragma HLS ARRAY_PARTITION variable=APFcluster complete dim=0
#pragma HLS ARRAY_PARTITION variable=PFclusterOutput complete dim=0
ap_uint<12> start, end;

        for(loop i=0; i<N_INPUT_LINKS; i++){
        start = 0; end = start+63 ;
                for(loop k=0; k<N_PF_LINK; k++){
                  APFcluster[i*N_PF_LINK+k].getPFcluster(link_in[i].range(end, start)) ;
                  start = start + 64 ; end = start + 63 ;
        }}
/* --------------
 * put your code here
 * and recalculate each PF cluster
 ------------------ */
        for(loop i=0 ; i <N_PF ; i++)
        {
        	PFclusterOutput[i] = APFcluster[i] ;
        	PFclusterOutput[i].ET = int(0.8 * PFclusterOutput[i].ET );
			PFclusterOutput[i].Eta = int(1 +  PFclusterOutput[i].Eta );
        }

        for(loop i=0; i<N_OUTPUT_LINKS; i++){
        start = 0; end = start+63 ;
                for(loop k=0; k<N_PF_LINK; k++){
                  link_out[i].range(end, start) = PFclusterOutput[i*N_PF_LINK+k].data()  ;
                  start = start + 64 ; end = start + 63 ;
        }}
}
