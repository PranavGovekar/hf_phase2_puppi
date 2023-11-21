#include "linpuppi.h"
#include "algo_topIP1.h"

void convertPFClustersToHadCaloObjs(PFcluster* pfclusters, l1ct::HadCaloObj* hadcalobjs) {
	for (int i = 0; i < N_PF; i++) {
		hadcalobjs[i].clear(); //to clear up garbage data from previous instance
		hadcalobjs[i].hwPt = pfclusters[i].ET;
		hadcalobjs[i].hwEta =  pfclusters[i].Eta & (0<<4 | 1<<3 | 1<<2 | 1<<1 | 1) ;
		if(  pfclusters[i].Eta & (1<<4))
			hadcalobjs[i].hwEta = -1 * hadcalobjs[i].hwEta ;

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

void algo_topIP1(ap_uint<576> link_in[N_INPUT_LINKS],ap_uint<576> link_out[N_OUTPUT_LINKS]){

	PFcluster APFcluster[N_PF];
	PFcluster PFclusterOutput[N_PF];

	#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
	#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
	#pragma HLS INTERFACE ap_ctrl_hs port=return
	#pragma HLS ARRAY_PARTITION variable=APFcluster complete dim=0
	#pragma HLS ARRAY_PARTITION variable=PFclusterOutput complete dim=N_SECTORS

	ap_uint<12> start, end;

	// empty array of PuppiObj 
	l1ct::PuppiObj pfselne[N_SECTORS][NNEUTRALS];

	//dummy PF Region Object
	l1ct::PFRegion region[N_SECTORS];
	region[0].hwEtaCenter = l1ct::glbeta_t(6);
	region[0].hwPhiCenter = l1ct::glbphi_t(17);
	region[0].hwEtaHalfWidth = l1ct::eta_t(5);
	region[0].hwPhiHalfWidth = l1ct::phi_t(18);
	region[0].hwEtaExtra = l1ct::eta_t(0);
	region[0].hwPhiExtra = l1ct::phi_t(2);

	region[1].hwEtaCenter = l1ct::glbeta_t(6);
	region[1].hwPhiCenter = l1ct::glbphi_t(54);
	region[1].hwEtaHalfWidth = l1ct::eta_t(5);
	region[1].hwPhiHalfWidth = l1ct::phi_t(18);
	region[1].hwEtaExtra = l1ct::eta_t(0);
	region[1].hwPhiExtra = l1ct::phi_t(2);

	region[2].hwEtaCenter = l1ct::glbeta_t(-6);
	region[2].hwPhiCenter = l1ct::glbphi_t(17);
	region[2].hwEtaHalfWidth = l1ct::eta_t(5);
	region[2].hwPhiHalfWidth = l1ct::phi_t(18);
	region[2].hwEtaExtra = l1ct::eta_t(0);
	region[2].hwPhiExtra = l1ct::phi_t(2);

	region[3].hwEtaCenter = l1ct::glbeta_t(-6);
	region[3].hwPhiCenter = l1ct::glbphi_t(54);
	region[3].hwEtaHalfWidth = l1ct::eta_t(5);
	region[3].hwPhiHalfWidth = l1ct::phi_t(18);
	region[3].hwEtaExtra = l1ct::eta_t(0);
	region[3].hwPhiExtra = l1ct::phi_t(2);

	// Convert Links to PFClusters
	for(loop i=0; i<N_INPUT_LINKS; i++){
		start = 0; end = start+63 ;
				for(loop k=0; k<N_PF_LINK; k++){
					APFcluster[i*N_PF_LINK+k].getPFcluster(link_in[i].range(end, start)) ;
					start = start + 64 ; end = start + 63 ;
		}
	}

	//Create HadCaloObj object and convert APFcluster
	l1ct::HadCaloObj H_in[N_PF];
	convertPFClustersToHadCaloObjs(APFcluster, H_in);

	l1ct::HadCaloObj H_in_regionized[N_SECTORS][NCALO];

	for(loop i=0; i<N_SECTORS; i++) {
		for(loop j=0; j<NCALO; j++) {
			H_in_regionized[i][j].clear();
		}
	}

	int nClusterInRegion[4]={0,0,0,0};
	for(loop i=0; i<N_PF; i++) {
		for(loop j=0; j<N_SECTORS; j++) {
			if(region[j].isInside(H_in[i].hwEta-region[j].hwEtaCenter, H_in[i].hwPhi-region[j].hwPhiCenter) )  {
				if( nClusterInRegion[j] < NCALO ) {
					H_in_regionized[j][nClusterInRegion[j]].hwEta = H_in[i].hwEta-region[j].hwEtaCenter;
					H_in_regionized[j][nClusterInRegion[j]].hwPhi = H_in[i].hwPhi-region[j].hwPhiCenter;
					H_in_regionized[j][nClusterInRegion[j]].hwPt = H_in[i].hwPt;
					nClusterInRegion[j]++;
				}
			}
		}
	}

	//PUPPI algorithm
	// Loop
	for(loop i=0; i<N_SECTORS; i++) {
		fwdlinpuppi(region[i], H_in_regionized[i], pfselne[i]);
	}


	//Convert PuppiObj to PFClusters
	for(loop i=0; i<N_SECTORS; i++) {
	convertPuppiObjsToPFClusters(pfselne[i], &PFclusterOutput[i*NNEUTRALS]);
	}
	// Convert PFclusterOutput back to link
	for(loop i=0; i<N_OUTPUT_LINKS; i++){
	start = 0; end = start+63 ;
			for(loop k=0; k<N_PF_LINK; k++){
				link_out[i].range(end, start) = PFclusterOutput[i*N_PF_LINK+k].data();
				start = start + 64 ; end = start + 63 ;
			}
	}
}
