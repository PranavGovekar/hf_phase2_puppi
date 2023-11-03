#include "linpuppi.h"
#include "algo_topIP1.h"


void convertPFClustersToHadCaloObjs(PFcluster* pfclusters, l1ct::HadCaloObj* hadcalobjs) {
	for (int i = 0; i < N_PF; i++) {
		hadcalobjs[i].clear();

		hadcalobjs[i].hwPt = pfclusters[i].ET;
		hadcalobjs[i].hwEta = pfclusters[i].Eta;
		hadcalobjs[i].hwPhi = pfclusters[i].Phi;
		}
}

void convertPuppiObjsToPFClusters(l1ct::PuppiObj* puppiobjs, PFcluster* pfclusters) {
	for (int i = 0; i < N_PF; i++) {
		pfclusters[i].ET = puppiobjs[i].hwPt;
		pfclusters[i].Eta = puppiobjs[i].hwEta;
		pfclusters[i].Phi = puppiobjs[i].hwPhi;

		pfclusters[i].Spare = 0;
		pfclusters[i].all = 0;
		}
}


void algo_topIP1(ap_uint<576> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS]){

	PFcluster APFcluster[N_PF];
	PFcluster PFclusterOutput[N_PF];

	#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
	#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
	#pragma HLS PIPELINE II=9
	#pragma HLS INTERFACE ap_ctrl_hs port=return
	#pragma HLS ARRAY_PARTITION variable=APFcluster complete dim=0
	#pragma HLS ARRAY_PARTITION variable=PFclusterOutput complete dim=0

	ap_uint<12> start, end;

	// empty array of PuppiObj 
	l1ct::PuppiObj pfselne[N_PF];

	//dummy PF Region Object
	l1ct::PFRegion region;
	region.hwEtaCenter = l1ct::glbeta_t(0);
	region.hwPhiCenter = l1ct::glbphi_t(0);
	region.hwEtaHalfWidth = l1ct::eta_t(100);
	region.hwPhiHalfWidth = l1ct::phi_t(100);
	region.hwEtaExtra = l1ct::eta_t(0);
	region.hwPhiExtra = l1ct::phi_t(0);

	// Convert Links to PFClusters
	for(loop i=0; i<N_INPUT_LINKS; i++){
		start = 0; end = start+63 ;
				for(loop k=0; k<N_PF_LINK; k++){
					APFcluster[i*N_PF_LINK+k].getPFcluster(link_in[i].range(end, start)) ;
					start = start + 64 ; end = start + 63 ;
		}
	}

	//Created HadCaloObj object and converted APFcluster
	l1ct::HadCaloObj H_in[N_PF];
	convertPFClustersToHadCaloObjs(APFcluster, H_in);

	//PUPPI algo
	//void fwdlinpuppi( const l1ct::PFRegion & region, const l1ct::HadCaloObj caloin[NCALO], l1ct::PuppiObj& pfselne[NNEUTRALS]);
	fwdlinpuppi(region, H_in, pfselne);

	//Converted PuppiObj to PFClusters
	convertPuppiObjsToPFClusters(pfselne, PFclusterOutput);

	// Covert PFclusterOutput back to link
	for(loop i=0; i<N_OUTPUT_LINKS; i++){
	start = 0; end = start+63 ;
			for(loop k=0; k<N_PF_LINK; k++){
				link_out[i].range(end, start) = PFclusterOutput[i*N_PF_LINK+k].data()  ;
				start = start + 64 ; end = start + 63 ;
	 }
	}

}
