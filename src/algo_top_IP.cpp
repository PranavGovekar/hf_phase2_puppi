#include "algo_topIP1.h"

//18 regions

// This function takes one link and makes 12x4 HFTowers Array
void processInputLink( ap_uint<LINK_WIDTH> link_in,
						hftower HFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI/N_INPUT_LINKS]){

#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=HFTowers complete dim=0

	for(loop j=0; j<(TOWERS_IN_PHI/N_INPUT_LINKS)/2; j=j+2){
		for(loop i=0; i<TOWERS_IN_ETA-1; i++){
		  hftower halfEnergy = hftower(link_in.range((i*10)+9, i*10));
		  halfEnergy.energy = halfEnergy.energy >> 1;

		  HFTowers[i][j] = halfEnergy;
		  HFTowers[i][j+1] = halfEnergy;

		  hftower halfEnergy_0 = hftower(link_in.range((i*10+110)+9, i*10+110));
		  halfEnergy_0.energy = halfEnergy_0.energy >> 1;

		  HFTowers[i][j+2] = halfEnergy_0;
		  HFTowers[i][j+3] = halfEnergy_0;
	  }
	}

    hftower A10 = HFTowers[TOWERS_IN_ETA-2][0] ;
    hftower B10 = HFTowers[TOWERS_IN_ETA-2][2] ;

    A10.energy = A10.energy >> 1 ;
    B10.energy = B10.energy >> 1 ;

    HFTowers[TOWERS_IN_ETA-2][0] = A10;
    HFTowers[TOWERS_IN_ETA-2][1] = A10;
    HFTowers[TOWERS_IN_ETA-2][2] = A10;
    HFTowers[TOWERS_IN_ETA-2][3] = A10;

    HFTowers[TOWERS_IN_ETA-1][0] = B10;
    HFTowers[TOWERS_IN_ETA-1][1] = B10;
    HFTowers[TOWERS_IN_ETA-1][2] = B10;
    HFTowers[TOWERS_IN_ETA-1][3] = B10;


#ifndef __SYNTHESIS__
	for(loop eta=0; eta<TOWERS_IN_ETA; eta++){
		for(loop phi=0;phi<TOWERS_IN_PHI/N_INPUT_LINKS; phi++) {
			std::cout<< "\t" << HFTowers[eta][phi].energy;
		}
		std::cout << endl;
	}
	std::cout << endl;
#endif
}


// this function takes in 3 links and makes a region of 3 wedges (1 + 2 extra) i.e eta:(12+2) and phi:(4+8)
// this does not affect the seed searching because only the required area is searched. search region is eta:(12) and phi:(4+2)
void makeRegion(const ap_uint<LINK_WIDTH> link_in[LINKS_PER_REGION],
	hftower HFRegions[TOWERS_IN_ETA + EXTRA_IN_ETA*2][(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2]){

#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=HFRegions complete dim=0

	for(loop link=0; link<LINKS_PER_REGION; link++){

		hftower HFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI/N_INPUT_LINKS];

		#pragma HLS ARRAY_PARTITION variable=HFTowers complete dim=0
		processInputLink(link_in[link], HFTowers);

		for(loop eta=0; eta<TOWERS_IN_ETA; eta++){
			for(loop phi=0;phi<TOWERS_IN_PHI/N_INPUT_LINKS; phi++) {
				HFRegions[eta+1][(link*(TOWERS_IN_PHI/N_INPUT_LINKS))+phi] = HFTowers[eta][phi];
			}
		}
	}

#ifndef __SYNTHESIS__
	for(loop eta=0; eta<TOWERS_IN_ETA + EXTRA_IN_ETA*2; eta++){
		for(loop phi=0;phi<(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2; phi++) {
			std::cout<< "\t" << HFRegions[eta][phi].energy;
		}
		std::cout << endl;
	}
	std::cout << endl;
#endif
}

void findMaxEnergyTowerInEta(const hftower EtaTowers[TOWERS_IN_ETA], ap_uint<5>& etaC){
#pragma HLS PIPELINE II=1
	hftower tempArray[16];
	ap_uint<5> index[16];

#pragma HLS ARRAY_PARTITION variable=tempArray complete dim=0
#pragma HLS ARRAY_PARTITION variable=index complete dim=0

	for(loop i=0; i<12; i++){
		tempArray[i] = EtaTowers[i];
		index[i] = i;
	}


	for(loop i=16; i>1; i=(i/2)){
		for(loop j=0; j < i/2; j++){
			if(tempArray[index[j*2]].energy < tempArray[index[(j*2) + 1]].energy){
				index[j] = index[(j*2)+1];
			}
			else{
				index[j] = index[j*2];
			}
		}
	}

	etaC = index[0];
}

void findMaxEnergyTowerInPhi(const hftower HFRegion[TOWERS_IN_ETA + EXTRA_IN_ETA*2][(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2],
								const ap_uint<5> etaCenters[6], ap_uint<8>& phiC){
#pragma HLS PIPELINE II=1
	hftower tempArray[8];
	ap_uint<4> index[8];

#pragma HLS ARRAY_PARTITION variable=tempArray complete dim=0
#pragma HLS ARRAY_PARTITION variable=index complete dim=0

	for(loop i=0; i<6; i++){
		tempArray[i] = HFRegion[etaCenters[i]][i+3];
		index[i] = i;
	}

	for(loop i=8; i>1; i=(i/2)){
		for(loop j=0; j < i/2; j++){
			if(tempArray[index[j*2]].energy < tempArray[index[(j*2) + 1]].energy){
				index[j] = index[(j*2)+1];
			}
			else{
				index[j] = index[j*2];
			}
		}
	}

	phiC = index[0];
}

void findMaxEnergyTower(hftower HFRegion[TOWERS_IN_ETA + EXTRA_IN_ETA*2][(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2],
				ap_uint<5>& etaC,
				ap_uint<8>& phiC){
//3,4,5,6,7,8
#pragma HLS ARRAY_PARTITION variable=HFRegion complete dim=0

	ap_uint<5> towersPhi[6];
#pragma HLS ARRAY_PARTITION variable=towersPhi complete dim=0
	ap_uint<8> tempPhi;
	for(ap_uint<5> phi = 3; phi < 9; phi++) {
		hftower towersEta[12];
#pragma HLS ARRAY_PARTITION variable=towersEta complete dim=0
		for(loop eta = 0; eta < 12; eta++) {
			towersEta[eta] = HFRegion[eta][phi];
		}
		findMaxEnergyTowerInEta(towersEta, towersPhi[phi-3]);
	}

	findMaxEnergyTowerInPhi(HFRegion, towersPhi, tempPhi);

	etaC = towersPhi[tempPhi];
	phiC = tempPhi + 3;
}

// Function to form the cluster and zero out the tower energies
void formClusterAndZeroOut(hftower HFRegion[TOWERS_IN_ETA + EXTRA_IN_ETA*2][(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2],
		ap_uint<5> etaC,
		ap_uint<8> phiC,
		ap_uint<12>& etaSum){

#pragma HLS ARRAY_PARTITION variable=HFRegion complete dim=0

    etaSum = HFRegion[etaC][phiC].energy +
             HFRegion[etaC+1][phiC].energy +
             HFRegion[etaC-1][phiC].energy +
             HFRegion[etaC][phiC+1].energy +
             HFRegion[etaC][phiC-1].energy +
             HFRegion[etaC+1][phiC+1].energy +
             HFRegion[etaC-1][phiC+1].energy +
             HFRegion[etaC+1][phiC-1].energy +
             HFRegion[etaC-1][phiC-1].energy;

    HFRegion[etaC][phiC].energy = 0;
    HFRegion[etaC+1][phiC].energy = 0;
    HFRegion[etaC-1][phiC].energy = 0;
    HFRegion[etaC][phiC+1].energy = 0;
    HFRegion[etaC][phiC-1].energy = 0;
    HFRegion[etaC+1][phiC+1].energy = 0;
    HFRegion[etaC-1][phiC+1].energy = 0;
    HFRegion[etaC+1][phiC-1].energy = 0;
    HFRegion[etaC-1][phiC-1].energy = 0;
}

// Main clustrizer function
void clustrizer (const ap_uint<LINK_WIDTH> regionLinks[3],
		PFcluster Clusters[N_PF_CLUSTERS]){

#pragma HLS PIPELINE
#pragma HLS ARRAY_PARTITION variable=Clusters complete dim=0

	hftower HFRegion[TOWERS_IN_ETA + EXTRA_IN_ETA*2][(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2];
	#pragma HLS ARRAY_PARTITION variable=HFRegion complete dim=0

	makeRegion(regionLinks, HFRegion);

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

    for(loop cluster = 0; cluster < N_PF_CLUSTERS; cluster++) {
        ap_uint<8> etmax = 0;
        ap_uint<8> phiC = 1;
        ap_uint<5> etaC = 1;
        ap_uint<12> etaSum = 0;

        findMaxEnergyTower(HFRegion, etaC, phiC);

        if(HFRegion[etaC][phiC].energy > MIN_CLUSTER_SEED_ENERGY){
			formClusterAndZeroOut(HFRegion, etaC, phiC, etaSum);

			if(phiC > 3 && phiC < 8){
				Clusters[cluster].Eta = etaC;
				Clusters[cluster].Phi = phiC;
				Clusters[cluster].ET = etaSum;
			}
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

void getTopPFClusters(PFcluster Clusters[N_PF_CLUSTERS]) {

#pragma HLS ARRAY_PARTITION variable=Clusters complete dim=0
    for (int i = 0; i < 8; i++) {
        for (int j = i + 1; j < N_PF_CLUSTERS*3; j++) {
            if (Clusters[j].ET > Clusters[i].ET) {
                PFcluster temp = Clusters[i];
                Clusters[i] = Clusters[j];
                Clusters[j] = temp;
            }
        }
    }
}

void packer(PFcluster Clusters[N_PF_CLUSTERS], const ap_uint<576>& link_out_sector, const ap_uint<7> sector) {

#pragma HLS ARRAY_PARTITION variable=Clusters complete dim=0

    ap_uint<10> start = 0;
    ap_uint<10> end = start + 63;
    PFcluster zero;
    packer_loop:
    for (int cluster = 0; cluster < 8; cluster++) {
        #pragma HLS UNROLL
        if (Clusters[cluster].ET > 0) {
//            Clusters[cluster].Eta += -EXTRA_IN_ETA;
//            Clusters[cluster].Phi = (TOWERS_IN_PHI / N_SECTORS) * sector + (Clusters[cluster].Phi - EXTRA_IN_PHI);

            link_out_sector.range(end, start) = Clusters[cluster].data();
        } else {
            link_out_sector.range(end, start) = 0 ;//zero.data();
        }
        start += 64;
        end = start + 63;
    }
}

void algo_topIP1(ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], ap_uint<576> link_out[6]){

#pragma HLS PIPELINE II=9
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS INTERFACE ap_ctrl_hs port=return

	const ap_uint<LINK_WIDTH> regionLinks[N_SECTORS][3] = {{link_in[17], link_in[0], link_in[1]},
															{link_in[0], link_in[1], link_in[2]},
															{link_in[1], link_in[2], link_in[3]},
															{link_in[2], link_in[3], link_in[4]},
															{link_in[3], link_in[4], link_in[5]},
															{link_in[4], link_in[5], link_in[6]},
															{link_in[5], link_in[6], link_in[7]},
															{link_in[6], link_in[7], link_in[8]},
															{link_in[7], link_in[8], link_in[9]},
															{link_in[8], link_in[9], link_in[10]},
															{link_in[9], link_in[10], link_in[11]},
															{link_in[10], link_in[11], link_in[12]},
															{link_in[11], link_in[12], link_in[13]},
															{link_in[12], link_in[13], link_in[14]},
															{link_in[13], link_in[14], link_in[15]},
															{link_in[14], link_in[15], link_in[16]},
															{link_in[15], link_in[16], link_in[17]},
															{link_in[16], link_in[17], link_in[0]},
	};

//	PFcluster Clusters[6][N_PF_CLUSTERS*LINKS_PER_REGION];
//	#pragma HLS ARRAY_PARTITION variable=Clusters complete dim=0
//
//	for(loop region=0; region<N_SECTORS; region++){
//#pragma HLS UNROLL factor=2
//		clustrizer(regionLinks[region], &Clusters[(region-2)/3][(region%3)*4]);
//
//		if(region%3 == 2){
//			getTopPFClusters(Clusters[(region-2)/3]);
//
//			#ifndef __SYNTHESIS__
//				std::cout<< "After Sorting : " << endl;
//				for(loop cluster=0; cluster<N_PF_CLUSTERS*3; cluster++){
//					std::cout<<"Cluster " << cluster << " E = "<< Clusters[(region-2)/3][cluster].ET
//							<<", center = ("<< Clusters[(region-2)/3][cluster].Eta <<","<< Clusters[(region-2)/3][cluster].Phi << ")\n";
//				}
//				std::cout<<endl;
//			#endif
//
//			packer(Clusters[(region-2)/3], link_out[(region-2)/3], (region-2)/3);
//		}
//	}


	PFcluster Clusters[6][N_PF_CLUSTERS*LINKS_PER_REGION];
	#pragma HLS ARRAY_PARTITION variable=Clusters complete dim=0

	link_unroll_loop:
	for(loop link=0; link<6; link++){

		sector_unroll_loop:
		for(loop sector=0; sector<3; sector++){

			PFcluster tempClusters[N_PF_CLUSTERS];
			#pragma HLS ARRAY_PARTITION variable=tempClusters complete dim=0

			clustrizer(regionLinks[(link*3) + sector], tempClusters);

			for(loop i=0; i<N_PF_CLUSTERS; i++){
				tempClusters[i].Eta += -EXTRA_IN_ETA;
				tempClusters[i].Phi = (tempClusters[i].Phi - 4) + ((link*3)+sector)*4;
				Clusters[link][(sector*N_PF_CLUSTERS)+i] = tempClusters[i];
			}
		}

		getTopPFClusters(Clusters[link]);

#ifndef __SYNTHESIS__
	std::cout<< "After Sorting : " << endl;
	for(loop cluster=0; cluster<N_PF_CLUSTERS*3; cluster++){
		std::cout<<"Cluster " << cluster << " E = "<< Clusters[link][cluster].ET
				<<", center = ("<< Clusters[link][cluster].Eta <<","<< Clusters[link][cluster].Phi << ")\n";
	}
	std::cout<<endl;
#endif

        packer(Clusters[link], link_out[link], link);
	}

}

