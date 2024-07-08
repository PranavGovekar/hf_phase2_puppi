#include "algo_topIP1.h"

//6 regions

// This function takes one link and makes 12x4 HFTowers Array
void processInputLink(const ap_uint<LINK_WIDTH> link_in,
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


// this function takes in 5 links and makes a region of 5 wedges (3 + 2 extra) i.e eta:(12+2) and phi:(12+8)
// this does not affect the seed searching because only the required area is searched. search region is eta:(12) and phi:(12+2)
void makeRegion(const ap_uint<LINK_WIDTH> link_in[LINKS_PER_REGION],
	hftower HFRegions[TOWERS_IN_ETA + EXTRA_IN_ETA*2][(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2]){

	//zero padding
	for(loop ele=0; ele<(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2; ele++){
		HFRegions[0][ele] = hftower();
		HFRegions[(TOWERS_IN_ETA + EXTRA_IN_ETA*2)-1][ele] = hftower();
	}

	for(loop link=0; link<LINKS_PER_REGION; link++){

		hftower HFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI/N_INPUT_LINKS];
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

// linear search
void findMaxEnergyTower(hftower HFRegion[TOWERS_IN_ETA + EXTRA_IN_ETA*2][(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2],
				ap_uint<8>& etmax,
				ap_uint<7>& etaC,
				ap_uint<7>& phiC){

    etmax = MIN_CLUSTER_SEED_ENERGY;

    for(loop eta = 1; eta < NTOWER_IN_ETA_PER_SECTOR - 1; eta++) {
        for(loop phi = 3; phi < NTOWER_IN_PHI_PER_SECTOR - 4; phi++) {
            if(etmax < HFRegion[eta][phi].energy) {
                etmax = HFRegion[eta][phi].energy;
                etaC = eta;
                phiC = phi;
            }
        }
    }
}

// Function to form the cluster and zero out the tower energies
void formClusterAndZeroOut(hftower HFRegion[TOWERS_IN_ETA + EXTRA_IN_ETA*2][(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2],
		ap_uint<7> etaC,
		ap_uint<7> phiC,
		ap_uint<10>& etaSum){

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

// Function to check if the cluster is inside the region and update the Clusters array
void isInside(PFcluster Clusters[N_PF_CLUSTERS],
		ap_uint<4>& clusterIdx,
		ap_uint<7> etaC,
		ap_uint<7> phiC,
		ap_uint<10> etaSum) {

    if(phiC > 3 && phiC < 16) {
        Clusters[clusterIdx].Eta = etaC;
        Clusters[clusterIdx].Phi = phiC;
        Clusters[clusterIdx].ET = etaSum;
        clusterIdx++;
    }
}

// Main clustrizer function
void clustrizer (hftower HFRegion[TOWERS_IN_ETA + EXTRA_IN_ETA*2][(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2],
		PFcluster Clusters[N_PF_CLUSTERS]){

#pragma HLS PIPELINE

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
    ap_uint<7> phiC = 1;
    ap_uint<7> etaC = 1;
    ap_uint<10> etaSum = 0;
    ap_uint<4> clusterIdx = 0;

    for(loop cluster = 0; cluster < N_PF_CLUSTERS; cluster++) {
        findMaxEnergyTower(HFRegion, etmax, etaC, phiC);
        formClusterAndZeroOut(HFRegion, etaC, phiC, etaSum);
        isInside(Clusters, clusterIdx, etaC, phiC, etaSum);
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
    for (int i = 0; i < 8; i++) {
        for (int j = i + 1; j < N_PF_CLUSTERS; j++) {
            if (Clusters[j].ET > Clusters[i].ET) {
                PFcluster temp = Clusters[i];
                Clusters[i] = Clusters[j];
                Clusters[j] = temp;
            }
        }
    }
}

void packer(PFcluster Clusters[N_PF_CLUSTERS], ap_uint<576>& link_out_sector, ap_uint<7> sector) {
    ap_uint<10> start = 0;
    ap_uint<10> end = start + 63;
    PFcluster zero;

    for (int cluster = 0; cluster < 8; cluster++) {
        #pragma HLS UNROLL
        if (Clusters[cluster].ET > 0) {
            Clusters[cluster].Eta += -EXTRA_IN_ETA;
            Clusters[cluster].Phi = (TOWERS_IN_PHI / N_SECTORS) * sector + (Clusters[cluster].Phi - EXTRA_IN_PHI);

            link_out_sector.range(end, start) = Clusters[cluster].data();
        } else {
            link_out_sector.range(end, start) = zero.data();
        }
        start += 64;
        end = start + 63;
    }
}

void algo_topIP1(ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], ap_uint<576> link_out[6]){

#pragma HLS PIPELINE
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS INTERFACE ap_ctrl_hs port=return

	const ap_uint<LINK_WIDTH> regionLinks[N_SECTORS][5] = {{link_in[17], link_in[0], link_in[1], link_in[2], link_in[3]},
															{link_in[2], link_in[3], link_in[4], link_in[5], link_in[6]},
															{link_in[5], link_in[6], link_in[7], link_in[8], link_in[9]},
															{link_in[8], link_in[9], link_in[10], link_in[11], link_in[12]},
															{link_in[11], link_in[12], link_in[13], link_in[14], link_in[15]},
															{link_in[14], link_in[15], link_in[16], link_in[17], link_in[0]}};

	hftower HFRegion[TOWERS_IN_ETA + EXTRA_IN_ETA*2][(TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2];
	PFcluster Clusters[N_PF_CLUSTERS];

#pragma HLS ARRAY_PARTITION variable=HFRegion complete dim=0
#pragma HLS ARRAY_PARTITION variable=Clusters complete dim=0


	for(loop sector=0; sector<N_SECTORS; sector++){

		makeRegion(regionLinks[sector], HFRegion);

		clustrizer(HFRegion, Clusters);

		getTopPFClusters(Clusters);

        packer(Clusters, link_out[sector], sector);
	}

}

