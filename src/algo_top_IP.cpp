#include "algo_topIP1.h"

void processInputLinks(ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], hftower HFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI]){
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=HFTowers complete dim=0

	hftower tmepHFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI/2];

	for(loop j=0; j<TOWERS_IN_PHI/2; j=j+2){
		for(loop i=0; i<TOWERS_IN_ETA-1; i++){
		  ap_uint<10> start   = i*10;
		  ap_uint<10> end = start + 9;
		  tmepHFTowers[i][j] = hftower(link_in[j/2].range(end, start));
	  }
		for(loop i=0; i<TOWERS_IN_ETA-1; i++){
		  ap_uint<10> start   = i*10+110;
		  ap_uint<10> end = start + 9;
		  tmepHFTowers[i][j+1] = hftower(link_in[j/2].range(end, start));
	  }
	}

    for(loop j=0; j<TOWERS_IN_PHI/2; j=j+2){
    hftower A10 = tmepHFTowers[TOWERS_IN_ETA-2][j] ;
    hftower B10 = tmepHFTowers[TOWERS_IN_ETA-2][j+1] ;

    ap_uint<8> halfA = A10.energy >> 1 ;
    ap_uint<8> halfB = B10.energy >> 1 ;


    A10.energy = halfA ;
    B10.energy = halfB ;


    tmepHFTowers[TOWERS_IN_ETA-2][j] = A10;
    tmepHFTowers[TOWERS_IN_ETA-2][j+1] = A10;

    tmepHFTowers[TOWERS_IN_ETA-1][j] = B10;
    tmepHFTowers[TOWERS_IN_ETA-1][j+1] = B10;
  }

    for(loop i=0; i<TOWERS_IN_ETA; i++){
		for(loop j=0; j<TOWERS_IN_PHI/2; j++){

			tmepHFTowers[i][j].energy = tmepHFTowers[i][j].energy >> 1;

			HFTowers[i][j*2] = tmepHFTowers[i][j];
			HFTowers[i][(j*2)+1] = tmepHFTowers[i][j];
		}
    }
}


void fillRegions (	const hftower HFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI],
					hftower HFRegions[N_SECTORS][TOWERS_IN_ETA + EXTRA_IN_ETA*2][(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2]){
#pragma HLS ARRAY_PARTITION variable=HFTowers complete dim=0
#pragma HLS ARRAY_PARTITION variable=HFRegions complete dim=0

	for(loop sector=0; sector<N_SECTORS; sector++){
		for(loop ele=0; ele<(TOWERS_IN_PHI/N_SECTORS) + TOWERS_IN_PHI; ele++){
			HFRegions[sector][0][ele] = hftower();
		}

		if(sector == 0){
			const ap_uint<7> phi_start = TOWERS_IN_PHI - EXTRA_IN_PHI;
			for(loop eta=1; eta<TOWERS_IN_ETA; eta++){
				ap_uint<7> phi_region = 0;
				for(loop phi=phi_start; phi<TOWERS_IN_PHI; phi++){
					HFRegions[sector][eta][phi_region] = HFTowers[eta][phi];
					phi_region++;
				}
			}

			const ap_uint<7> phi_start_ = ((TOWERS_IN_PHI/N_SECTORS)*sector);
			for(loop eta=1; eta<TOWERS_IN_ETA; eta++){
				ap_uint<7> phi_region = 2;
				for(loop phi=phi_start_; phi<phi_start_ + (TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI; phi++){
					HFRegions[sector][eta][phi_region] = HFTowers[eta][phi];
					phi_region++;
				}
			}

		}

		else if(sector == N_SECTORS-1){
			const ap_uint<7> phi_start = ((TOWERS_IN_PHI/N_SECTORS)*sector) - EXTRA_IN_PHI;
			for(loop eta=1; eta<TOWERS_IN_ETA; eta++){
				ap_uint<7> phi_region = 0;
				for(loop phi=phi_start; phi<phi_start + (TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI; phi++){
					HFRegions[sector][eta][phi_region] = HFTowers[eta][phi];
					phi_region++;
				}
			}

			const ap_uint<7> phi_start_ = 0;
			for(loop eta=1; eta<TOWERS_IN_ETA; eta++){
				ap_uint<7> phi_region = 14;
				for(loop phi=phi_start_; phi<EXTRA_IN_PHI; phi++){
					HFRegions[sector][eta][phi_region] = HFTowers[eta][phi];
					phi_region++;
				}
			}

		}

		else{
			const ap_uint<7> phi_start = ((TOWERS_IN_PHI/N_SECTORS)*sector) - EXTRA_IN_PHI;
			for(loop eta=1; eta<TOWERS_IN_ETA; eta++){
				ap_uint<7> phi_region = 0;
				for(loop phi=phi_start; phi<phi_start + (TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2; phi++){
					HFRegions[sector][eta][phi_region] = HFTowers[eta][phi];
					phi_region++;
				}
			}
		}

		for(loop ele=0; ele<(TOWERS_IN_PHI/N_SECTORS) + TOWERS_IN_PHI; ele++){
			HFRegions[sector][TOWERS_IN_ETA + EXTRA_IN_ETA][ele] = hftower();
		}
	}

#ifndef __SYNTHESIS__
	for(loop sector=0; sector<N_SECTORS; sector++){
		for(loop eta=0; eta<TOWERS_IN_ETA + EXTRA_IN_ETA*2; eta++){
			for(loop phi=0;phi<(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2; phi++) {
				std::cout<< "\t" << HFRegions[sector][eta][phi].energy;
			}
			std::cout << endl;
		}
		std::cout << endl << endl;
	}
#endif
}

void clustrizer (	hftower HFRegion[TOWERS_IN_ETA + EXTRA_IN_ETA*2][(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2],
					PFcluster Clusters[N_PF_CLUSTERS]){

#pragma HLS ARRAY_PARTITION variable=HFRegion complete dim=0
#pragma HLS ARRAY_PARTITION variable=Clusters complete dim=0

#ifndef __SYNTHESIS__
	for(loop phiI=0;phiI<(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2; phiI++){
		std::cout<<std::setw(5)<<phiI;
	}
	std::cout<<"  <--Phi"<<"\n";
	for(loop etaI=0;etaI<TOWERS_IN_ETA + EXTRA_IN_ETA*2; etaI++){
		for(loop phiI=0;phiI<(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2; phiI++){
			std::cout<<std::setw(5)<< HFRegion[etaI][phiI].energy;
		}
		std::cout<<"  <--"<<etaI<<"\n";
	}
	std::cout<<endl;
#endif

	ap_uint<8> etmax = 0;
	ap_uint<7> phiC = 0;
	ap_uint<7> etaC = 0;
	ap_uint<10> etaSum = 0;
	ap_uint<4> clusterIdx = 0;

	for(loop cluster=0; cluster<N_PF_CLUSTERS; cluster++){
		for(loop eta=1;eta<(TOWERS_IN_ETA + EXTRA_IN_ETA*2)-1; eta++){
			for(loop phi=1;phi<((TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2)-1; phi++){
				if(etmax < HFRegion[eta][phi].energy){
					etmax = HFRegion[eta][phi].energy;
					etaC = eta;
					phiC = phi;
				}
			}
		}

		if (etmax <= MIN_CLUSTER_SEED_ENERGY){
			continue;
		}
		else{
			etaSum += HFRegion[etaC][phiC].energy;
			HFRegion[etaC][phiC].energy = 0;

			etaSum += HFRegion[etaC+1][phiC].energy;
			HFRegion[etaC+1][phiC].energy = 0;

			etaSum += HFRegion[etaC-1][phiC].energy;
			HFRegion[etaC-1][phiC].energy = 0;

			etaSum += HFRegion[etaC][phiC+1].energy;
			HFRegion[etaC][phiC+1].energy = 0;

			etaSum += HFRegion[etaC][phiC-1].energy;
			HFRegion[etaC][phiC-1].energy = 0;


			etaSum += HFRegion[etaC+1][phiC+1].energy;
			HFRegion[etaC+1][phiC+1].energy = 0;

			etaSum += HFRegion[etaC-1][phiC+1].energy;
			HFRegion[etaC-1][phiC+1].energy = 0;

			etaSum += HFRegion[etaC+1][phiC-1].energy;
			HFRegion[etaC+1][phiC-1].energy = 0;

			etaSum += HFRegion[etaC-1][phiC-1].energy;
			HFRegion[etaC-1][phiC-1].energy = 0;
			
			if(phiC > 1 && phiC < 14){
				Clusters[clusterIdx].Eta = etaC;
				Clusters[clusterIdx].Phi = phiC;
				Clusters[clusterIdx].ET = etaSum;

				etaSum = 0;
				clusterIdx++;
			}
			etmax = 0;
		}
	}

#ifndef __SYNTHESIS__
	for(loop cluster=0; cluster<N_PF_CLUSTERS; cluster++){
		std::cout<<"Cluster " << cluster << " E = "<< Clusters[cluster].ET
				<<", center = ("<< Clusters[cluster].Eta <<","<< Clusters[cluster].Phi << ")\n";
	}
	std::cout<<endl;
#endif
}

void algo_topIP1(ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], ap_uint<576> link_out[6]){
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS INTERFACE ap_ctrl_hs port=return

	hftower HFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI];
	hftower regions [N_SECTORS][TOWERS_IN_ETA + EXTRA_IN_ETA*2][(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2];

#pragma HLS ARRAY_PARTITION variable=HFTowers complete dim=0
#pragma HLS ARRAY_PARTITION variable=regions complete dim=0

	processInputLinks(link_in, HFTowers);

	fillRegions(HFTowers, regions);
	
	PFcluster Clusters[N_SECTORS][N_PF_CLUSTERS];

	for(loop sector=0; sector<N_SECTORS; sector++){
		clustrizer(regions[sector],Clusters[sector]);
		ap_uint<10> start = 0;
		ap_uint<10> end = start+63;

		for(loop cluster=0; cluster<9; cluster++){
				Clusters[sector][cluster].Eta += -EXTRA_IN_ETA;
				Clusters[sector][cluster].Phi += sector*TOWERS_IN_PHI;

				link_out[sector].range(end, start) = Clusters[sector][cluster].data();
			}

            start = start + 64 ;
            end = start + 63 ;

	}
}

