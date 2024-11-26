#include "pfClusterering.h"

// this function takes in 3 links and makes a region of 3 wedges (1 + 2 extra) i.e eta:(12+2) and phi:(4+8)
// this does not affect the seed searching because only the required area is searched. search region is eta:(12) and phi:(4+2)
void makeRegion(const ap_uint<LINK_WIDTH> link_in[LINKS_PER_REGION],
                hftower HFRegions[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR])
{

//    #pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
//    #pragma HLS ARRAY_PARTITION variable=HFRegions complete dim=0

	makeRegion_1_1:
    for(loop link=0; link<LINKS_PER_REGION; link++)
    {

        hftower HFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI/N_INPUT_LINKS];

//        #pragma HLS ARRAY_PARTITION variable=HFTowers complete dim=0
        processInputLink(link_in[link], HFTowers);

        makeRegion_1_2:
        for(loop eta=0; eta<TOWERS_IN_ETA; eta++)
        {
        	makeRegion_1_3:
            for(loop phi=0; phi<TOWERS_IN_PHI/N_INPUT_LINKS; phi++)
            {
                HFRegions[eta+1][(link*(TOWERS_IN_PHI/N_INPUT_LINKS))+phi] = HFTowers[eta][phi];

                HFRegions[eta+1][(link*(TOWERS_IN_PHI/N_INPUT_LINKS))+phi].eta = eta+1;
                HFRegions[eta+1][(link*(TOWERS_IN_PHI/N_INPUT_LINKS))+phi].phi = 4*link + phi;
            }
        }
    }

#ifndef __SYNTHESIS__

    if (DEBUG_LEVEL > 9)
    {
        for(loop eta=0; eta<NTOWER_IN_ETA_PER_SECTOR; eta++)
        {
            for(loop phi=0; phi<NTOWER_IN_PHI_PER_SECTOR; phi++)
            {
                std::cout<< "\t" << HFRegions[eta][phi].energy;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
#endif
}

void findMaxEnergyTowerInPhi(const hftower EtaTowers[6],
							hftower& phiC) {
#pragma HLS INLINE off
    hftower tempArray[8];
//    #pragma HLS ARRAY_PARTITION variable=tempArray complete dim=0

    findMaxEnergyTowerInPhi_1_1:
    for(loop i=0; i<6; i++){
        tempArray[i] = EtaTowers[i];
    }

    findMaxEnergyTowerInPhi_2_1:
    for(loop i=8; i>1; i=(i/2))
    {
    	findMaxEnergyTowerInPhi_2_2:
        for(loop j=0; j < i/2; j++)
        {
        	tempArray[j] = bestOf2(tempArray[j*2], tempArray[(j*2) + 1]);
        }
    }

    phiC = tempArray[0];
}


void findMaxEnergyTowerInEta(const hftower EtaTowers[TOWERS_IN_ETA], hftower& etaC) {
#pragma HLS INLINE off
	hftower tempArray[16];
//    #pragma HLS ARRAY_PARTITION variable=tempArray complete dim=0

	findMaxEnergyTowerInEta_1_1:
    for(loop i=0; i<12; i++)
    {
        tempArray[i] = EtaTowers[i];
    }

    findMaxEnergyTowerInEta_2_1:
    for(loop i=16; i>1; i=(i/2))
    {
    	findMaxEnergyTowerInEta_2_2:
        for(loop j=0; j < i/2; j++)
        {
        	tempArray[j] = bestOf2(tempArray[j*2], tempArray[(j*2) + 1]);
        }
    }

    etaC = tempArray[0];

}


void findMaxEnergyTower(const hftower HFRegion[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR],
                        ap_uint<5>& etaC,
                        ap_uint<8>& phiC)
{
#pragma HLS PIPELINE
//    #pragma HLS ARRAY_PARTITION variable=HFRegion complete dim=0

	hftower towersEta[12];
//    #pragma HLS ARRAY_PARTITION variable=towersPhi complete dim=0
	hftower maxTower;

	findMaxEnergyTower_1_1:
	for(loop eta = 0; eta < 12; eta++){
		hftower towersPhi[6];
		findMaxEnergyTower_1_2:
		for(ap_uint<5> phi = 0; phi < 6; phi++){
			towersPhi[phi] = HFRegion[eta+1][phi+3];
		}
		findMaxEnergyTowerInPhi(towersPhi, towersEta[eta]);
	}
	findMaxEnergyTowerInEta(towersEta, maxTower);

    etaC = maxTower.eta;
    phiC = maxTower.phi;


#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL >7 )
    {
        std::cout<<"  Center found at : "<<etaC<<" , "<<phiC<<"\n" ;
    }
#endif
}


void formClusterAndZeroOut(hftower HFRegion[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR],
                           ap_uint<5> etaC,
                           ap_uint<8> phiC,
                           ap_uint<12>& etSum) {

	ap_uint<12> etSum__ = 0;
//    #pragma HLS ARRAY_PARTITION variable=HFRegion complete dim=0

    ap_uint<1> mask[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR] = {0};
    ap_uint<1> zero[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR] = {0};

    formClusterAndZeroOut_1_1:
    for (ap_uint<5> i = 0; i < NTOWER_IN_ETA_PER_SECTOR; i++) {
        if (i + 1 >= etaC && i <= etaC + 1) {
        	formClusterAndZeroOut_1_2:
            for (ap_uint<8> j = 0; j < NTOWER_IN_PHI_PER_SECTOR; j++) {
                if (j + 1 >= phiC && j <= phiC + 1) {
                    mask[i][j] = 1;
                }
            }
        }
    }

    formClusterAndZeroOut_2_1:
    for (ap_uint<5> i = 0; i < NTOWER_IN_ETA_PER_SECTOR; i++) {
    	formClusterAndZeroOut_2_2:
        for (ap_uint<8> j = 0; j < NTOWER_IN_PHI_PER_SECTOR; j++) {
        	etSum__ += HFRegion[i][j].energy * mask[i][j];
        }
    }

    formClusterAndZeroOut_3_1:
    for (ap_uint<5> i = 0; i < NTOWER_IN_ETA_PER_SECTOR; i++) {
    	formClusterAndZeroOut_3_2:
        for (ap_uint<8> j = 0; j < NTOWER_IN_PHI_PER_SECTOR; j++) {
            zero[i][j] = (ap_uint<1>)(1 - mask[i][j]);
        }
    }

    formClusterAndZeroOut_4_1:
    for (ap_uint<5> i = 0; i < NTOWER_IN_ETA_PER_SECTOR; i++) {
    	formClusterAndZeroOut_4_2:
        for (ap_uint<8> j = 0; j < NTOWER_IN_PHI_PER_SECTOR; j++) {
            HFRegion[i][j].energy *= zero[i][j];
        }
    }

    etSum = etSum__;
}

void swap_1(const PFcluster Clusters_in[N_PF_CLUSTERS], PFcluster Clusters_out[N_PF_CLUSTERS]){
#pragma HLS INLINE off
    GreaterSmaller res;

    res = AscendDescend(Clusters_in[0], Clusters_in[1]);
    Clusters_out[0] = res.greater;
    Clusters_out[1] = res.smaller;

    res = AscendDescend(Clusters_in[2], Clusters_in[3]);
    Clusters_out[2] = res.greater;
    Clusters_out[3] = res.smaller;
}

void swap_2(const PFcluster Clusters_in[N_PF_CLUSTERS], PFcluster Clusters_out[N_PF_CLUSTERS]){
#pragma HLS INLINE off
    GreaterSmaller res;

    res = AscendDescend(Clusters_in[0], Clusters_in[3]);
    Clusters_out[0] = res.greater;
    Clusters_out[3] = res.smaller;

    res = AscendDescend(Clusters_in[1], Clusters_in[2]);
    Clusters_out[1] = res.greater;
    Clusters_out[2] = res.smaller;
}

void sortPFcluster_4(const PFcluster Clusters_in[N_PF_CLUSTERS], PFcluster Clusters_out[N_PF_CLUSTERS]){
#pragma HLS INLINE off
	PFcluster temp_1[N_PF_CLUSTERS];
	PFcluster temp_2[N_PF_CLUSTERS];

	swap_1(Clusters_in, temp_1);
	swap_2(temp_1, temp_2);
	swap_1(temp_2, Clusters_out);

}

void findMaxEnergyPFcluster_32(const PFcluster Clusters[18], ap_uint<5>& maxIdx) {
		PFcluster tempArray[32];
		ap_uint<6> index[32];

//	    #pragma HLS ARRAY_PARTITION variable=tempArray complete dim=0
//	    #pragma HLS ARRAY_PARTITION variable=index complete dim=0

		findMaxEnergyPFCluster_1_1:
	    for(loop i=0; i<18; i++)
	    {
	        tempArray[i] = Clusters[i];
	        index[i] = i;
	    }

	    findMaxEnergyPFCluster_3_1:
	    for(loop i=32; i>1; i=(i/2))
	    {
	    	findMaxEnergyPFCluster_3_2:
	        for(loop j=0; j < i/2; j++)
	        {
	            if(tempArray[index[j*2]].ET < tempArray[index[(j*2) + 1]].ET)
	            {
	                index[j] = index[(j*2)+1];
	            }
	            else
	            {
	                index[j] = index[j*2];
	            }
	        }
	    }

	    maxIdx = index[0];
	}


// Main clustrizer function
void makeCaloClusters (const ap_uint<LINK_WIDTH> regionLinks[LINKS_PER_REGION],
		PFcluster Clusters[N_PF_CLUSTERS],
		PFcluster EGClusters[N_PF_CLUSTERS],
		const ap_uint<10> sector)
{

//    #pragma HLS ARRAY_PARTITION variable=Clusters complete dim=0

    hftower HFRegion[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR];
//    #pragma HLS ARRAY_PARTITION variable=HFRegion complete dim=0

    PFcluster Clusters_in [N_PF_CLUSTERS];
    PFcluster EGClusters_in [N_PF_CLUSTERS];

    makeRegion(regionLinks, HFRegion);

#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 7)
    {
        for(loop phiI=0; phiI<NTOWER_IN_PHI_PER_SECTOR; phiI++)
        {
            std::cout<<std::setw(7)<<phiI;
        }
        std::cout<<"  <--Phi"<<"\n";
        for(loop etaI=0; etaI<NTOWER_IN_ETA_PER_SECTOR; etaI++)
        {
            for(loop phiI=0; phiI<NTOWER_IN_PHI_PER_SECTOR; phiI++)
            {
                std::cout<<std::setw(5)<< HFRegion[etaI][phiI].energy<<" | ";
            }
            std::cout<<"  <--"<<etaI-1<<"\n";
        }
        std::cout<<std::endl;
    }
#endif

    makeCaloClusters_1_1:
    for(loop cluster = 0; cluster < 6; cluster++)
    {
        ap_uint<8> phiC=0;
        ap_uint<5> etaC=0;
        ap_uint<12> etSum=0;
        ap_uint<10> seedET=0;

        findMaxEnergyTower(HFRegion, etaC, phiC);

        if(HFRegion[etaC][phiC].energy > MIN_CLUSTER_SEED_ENERGY) {

        	seedET = HFRegion[etaC][phiC].energy;
            formClusterAndZeroOut(HFRegion, etaC, phiC, etSum);

#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 7)
    {
        for(loop phiI=0; phiI<NTOWER_IN_PHI_PER_SECTOR; phiI++)
        {
            std::cout<<std::setw(7)<<phiI;
        }
        std::cout<<"  <--Phi"<<"\n";
        for(loop etaI=0; etaI<NTOWER_IN_ETA_PER_SECTOR; etaI++)
        {
            for(loop phiI=0; phiI<NTOWER_IN_PHI_PER_SECTOR; phiI++)
            {
                std::cout<<std::setw(5)<< HFRegion[etaI][phiI].energy<<" | ";
            }
            std::cout<<"  <--"<<etaI-1<<"\n";
        }
        std::cout<<std::endl;
    }
#endif

            if(phiC > 3 && phiC < 8) {
            	Clusters_in[cluster].Eta = etaC-1;
            	Clusters_in[cluster].Phi = phiC + sector*4;
            	Clusters_in[cluster].ET = etSum;

                if(etSum == 2*seedET) {
                	EGClusters_in[cluster].Eta = etaC-1;
                	EGClusters_in[cluster].Phi = phiC + sector*4;
                	EGClusters_in[cluster].ET = etSum;
                }
            }
        }
    }

    sortPFcluster_4(EGClusters_in, EGClusters);
    sortPFcluster_4(Clusters_in, Clusters);

#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 7)
    {
        for(loop cluster=0; cluster<N_PF_CLUSTERS; cluster++)
        {
            std::cout<<"Cluster " << cluster << " E = "<< Clusters_in[cluster].ET
                     <<", center = ("<< Clusters_in[cluster].Eta <<","<< Clusters_in[cluster].Phi << ")\n";
        }
        std::cout<<std::endl;
        for(loop cluster=0; cluster<N_PF_CLUSTERS; cluster++)
        {
            std::cout<<"Cluster " << cluster << " E = "<< Clusters[cluster].ET
                     <<", center = ("<< Clusters[cluster].Eta <<","<< Clusters[cluster].Phi << ")\n";
        }
        std::cout<<std::endl;
    }
#endif
}

void findMaxEnergyPFCluster(const ap_uint<12> EtaTowers[72], ap_uint<8>& etaC)
{
	ap_uint<12> tempArray[128];
	ap_uint<8> index[128];

//    #pragma HLS ARRAY_PARTITION variable=tempArray complete dim=0
//    #pragma HLS ARRAY_PARTITION variable=index complete dim=0

	findMaxEnergyPFCluster_1_1:
    for(loop i=0; i<72; i++)
    {
        tempArray[i] = EtaTowers[i];
        index[i] = i;
    }
    findMaxEnergyPFCluster_2_1:
    for(loop i=72; i<128; i++)
    {
        tempArray[i] = 0;
        index[i] = 127;
    }


    findMaxEnergyPFCluster_3_1:
    for(loop i=128; i>1; i=(i/2))
    {
    	findMaxEnergyPFCluster_3_2:
        for(loop j=0; j < i/2; j++)
        {
            if(tempArray[index[j*2]] < tempArray[index[(j*2) + 1]])
            {
                index[j] = index[(j*2)+1];
            }
            else
            {
                index[j] = index[j*2];
            }
        }
    }

    etaC = index[0];
}


void selectEGClusters(const PFcluster caloClusters[N_SECTORS][4] ,l1ct::HadCaloObj egClusters[9] )
{
    PFcluster  egCandidateClusters[N_SECTORS];
    ap_uint<3> cluaterIdx[N_SECTORS];
//#pragma HLS ARRAY_PARTITION variable=egCandidateClusters complete dim=0
//#pragma HLS ARRAY_PARTITION variable=cluaterIdx complete dim=0
    
    for(int sector=0 ; sector < N_SECTORS ; sector++){
    	egCandidateClusters[sector] = caloClusters[sector][0];
    	cluaterIdx[sector] = 1;
    }
    
    for(int cluster=0 ; cluster < 9 ; cluster++){
    	ap_uint<5> maxIdx;

    	findMaxEnergyPFcluster_32(egCandidateClusters, maxIdx);

    	egClusters[cluster].hwEta = egCandidateClusters[maxIdx].Eta;
    	egClusters[cluster].hwPhi = egCandidateClusters[maxIdx].Phi;
    	egClusters[cluster].hwPt = egCandidateClusters[maxIdx].ET;
    	egClusters[cluster].hwPt >>= 2;

    	egCandidateClusters[maxIdx] = caloClusters[maxIdx][cluaterIdx[maxIdx]];
    	cluaterIdx[maxIdx]++;
    }
}


void pfExy(const l1ct::PuppiObj pfselne[N_SECTORS_PF][NNEUTRALS],
		ap_fixed<32,16>& Ex, ap_fixed<32,16>& Ey){

#pragma HLS INLINE off

    ap_fixed<32,16> Ex_temp = 0;
    ap_fixed<32,16> Ey_temp = 0;

	ap_fixed<32,16> sin_LUT[72] = {0.0436, 0.1305, 0.2164, 0.3007, 0.3827, 0.4617, 0.5373, 0.6088,
            0.6756, 0.7373, 0.7934, 0.8434, 0.8870, 0.9239, 0.9537, 0.9763,
            0.9914, 0.9990, 0.9990, 0.9914, 0.9763, 0.9537, 0.9239, 0.8870,
            0.8434, 0.7934, 0.7373, 0.6756, 0.6088, 0.5373, 0.4617, 0.3827,
            0.3007, 0.2164, 0.1305, 0.0436, -0.0436, -0.1305, -0.2164, -0.3007,
            -0.3827, -0.4617, -0.5373, -0.6088, -0.6756, -0.7373, -0.7934, -0.8434,
            -0.8870, -0.9239, -0.9537, -0.9763, -0.9914, -0.9990, -0.9990, -0.9914,
            -0.9763, -0.9537, -0.9239, -0.8870, -0.8434, -0.7934, -0.7373, -0.6756,
            -0.6088, -0.5373, -0.4617, -0.3827, -0.3007, -0.2164, -0.1305, -0.0436};

//#pragma HLS ARRAY_PARTITION variable=sin_LUT complete dim=0

    for(int i=0; i < N_SECTORS_PF; i++ ){
        for(int j =0 ; j< NNEUTRALS ; j++){

//			#pragma HLS RESOURCE variable=pfselne[i][j].hwPt * sin_LUT[(pfselne[i][j].hwPhi + 18) % 72] core=Mul_LUT
            Ex_temp += pfselne[i][j].hwPt * sin_LUT[(pfselne[i][j].hwPhi + 18)%72];

//			#pragma HLS RESOURCE variable=pfselne[i][j].hwPt * sin_LUT[pfselne[i][j].hwPhi] core=Mul_LUT
            Ey_temp += pfselne[i][j].hwPt * sin_LUT[pfselne[i][j].hwPhi];
        }
    }

//#ifndef __SYNTHESIS__
//std::cout << "\ntemp Ex: " << Ex_temp << "\nEy: " << Ey_temp << std::endl;
//#endif

    Ex = Ex_temp;
    Ey = Ey_temp;
}

void doPFClustringChain( const ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS],
		l1ct::HadCaloObj egClusters[8],
		l1ct::PuppiObj pfSelectedNutrals[N_SECTORS_PF][NNEUTRALS],
		ap_fixed<32,16>& Ex,
		ap_fixed<32,16>& Ey)
{
#pragma HLS PIPELINE II=9

//=== define the  sector boundaries and overlaps ========================================
    l1ct::PFRegion region[N_SECTORS_PF];
    doPFClustringChain_2_1:
    for(int i=0 ; i < N_SECTORS_PF ; i++)
        {
            region[i].hwEtaCenter = l1ct::glbeta_t(6);
            region[i].hwPhiCenter = l1ct::glbphi_t(12*i+5);

            region[i].hwEtaExtra = l1ct::eta_t(0);
            region[i].hwPhiExtra = l1ct::phi_t(2);

            region[i].hwEtaHalfWidth = l1ct::eta_t(7);
            region[i].hwPhiHalfWidth = l1ct::phi_t(6);
            
            #ifndef __SYNTHESIS__
            if(DEBUG_LEVEL > 9)
                {
                    std::cout<<"Sector "<<i<<" | eta center : "<<region[i].hwEtaCenter<<" , phi center : "<<region[i].hwPhiCenter
                         <<" [ eta HW : "<<region[i].hwEtaHalfWidth<<" phi HW : "<< region[i].hwPhiHalfWidth  <<" ] "   
                         <<"\n";
                }
            #endif
        }
//=======================================================================================


    PFcluster caloClusters[N_SECTORS][4];
    PFcluster EGClusters[N_SECTORS][4];
//#pragma HLS ARRAY_PARTITION variable=caloClusters complete dim=0
//#pragma HLS ARRAY_PARTITION variable=EGClusters complete dim=0

    ap_uint<LINK_WIDTH> linksInSector[N_SECTORS][3] ;

    linksInSector[0][0] = link_in[17];
    linksInSector[0][1] = link_in[ 0];
    linksInSector[0][2] = link_in[ 1];
    doPFClustringChain_1_1:
    for(int i=1; i<N_SECTORS-1; i++) {
        linksInSector[i][0] = link_in[i-1];
        linksInSector[i][1] = link_in[i  ];
        linksInSector[i][2] = link_in[i+1];
    }
    linksInSector[N_SECTORS -1][0] = link_in[N_SECTORS -2 ];
    linksInSector[N_SECTORS -1][1] = link_in[N_SECTORS -1];
    linksInSector[N_SECTORS -1][2] = link_in[ 0];

    doPFClustringChain_3_1:
    for(loop sector=0; sector<N_SECTORS; sector++) {
		makeCaloClusters(linksInSector[sector], caloClusters[sector],EGClusters[sector], sector);
    }
    
#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 0)
    {
    for(loop sector=0; sector<N_SECTORS_PF; sector++)
        {
            std::cout<<"@@CALOCLUSTER | "<< sector<<" | ";
            for(int wedge =0 ; wedge<LINKS_PER_REGION ; wedge++)
            {
                for(int i =0 ; i <4 ;i++)
                 {
               //     std::cout<<caloClusters[sector*LINKS_PER_REGION+wedge][i].ET<<","
               //              <<caloClusters[sector*LINKS_PER_REGION+wedge][i].Eta<<"," 
               //              <<caloClusters[sector*LINKS_PER_REGION+wedge][i].Phi<<"" 
               //              <<" | " ;
                 std::cout<<caloClusters[sector*LINKS_PER_REGION+wedge][i].data()<<" | " ;
                }
            }
            std::cout<<"\n";
        }

    for(loop sector=0; sector<N_SECTORS_PF; sector++)
        {
            std::cout<<"CALOCLUSTER | "<< sector<<"\n";
            for(int wedge =0 ; wedge<LINKS_PER_REGION ; wedge++)
            {
                for(int i =0 ; i <4 ;i++)
                 {
                    std::cout<<caloClusters[sector*LINKS_PER_REGION+wedge][i].ET<<","
                             <<caloClusters[sector*LINKS_PER_REGION+wedge][i].Eta<<","
                             <<caloClusters[sector*LINKS_PER_REGION+wedge][i].Phi<<"\n";
                }
            }
            std::cout<<"\n";
        }
    }
#endif

    l1ct::HadCaloObj egCandidates[9];
//#pragma HLS ARRAY_PARTITION variable=egCandidates complete dim=0

    l1ct::PuppiObj pfselne[N_SECTORS_PF][NNEUTRALS];
    l1ct::HadCaloObj pfHadronicClusters[N_SECTORS_PF][NCALO];
//#pragma HLS ARRAY_PARTITION variable=pfselne complete dim=0
//#pragma HLS ARRAY_PARTITION variable=pfHadronicClusters complete dim=0

    for(loop sector=0; sector<N_SECTORS_PF; sector++){
    	for(loop wedge=0; wedge<LINKS_PER_REGION; wedge++){
    		for(loop cluster=0; cluster<N_PF_CLUSTERS; cluster++){
    			pfHadronicClusters[sector][wedge*N_PF_CLUSTERS + cluster].hwEta =
    					caloClusters[sector*LINKS_PER_REGION + wedge][cluster].Eta - region[sector].hwEtaCenter;

    			pfHadronicClusters[sector][wedge*N_PF_CLUSTERS + cluster].hwPhi =
    			    					caloClusters[sector*LINKS_PER_REGION + wedge][cluster].Phi - region[sector].hwPhiCenter;

    			pfHadronicClusters[sector][wedge*N_PF_CLUSTERS + cluster].hwPt =
    			    					caloClusters[sector*LINKS_PER_REGION + wedge][cluster].ET;

    			pfHadronicClusters[sector][wedge*N_PF_CLUSTERS + cluster].hwPt >>= 2;
    		}
    	}
    }


    for(loop sector=0; sector<N_SECTORS_PF; sector++){

    	int phi_offset_l;
    	int phi_offset_r;

    	if(sector == 0){
    		phi_offset_l = -72;
    		phi_offset_r = 0;
    	}
    	else if(sector == 5){
    		phi_offset_l = 0;
			phi_offset_r = 72;
    	}
    	else{
    		phi_offset_l = 0;
			phi_offset_r = 0;
    	}


    	PFcluster left;
    	PFcluster right;

    	int idx_left = 0;
    	int idx_right = 0;

    	left = caloClusters[(sector*3)-1][idx_left];
    	right = caloClusters[(sector*3)+3][idx_right];

    	for(loop idx=0; idx<4; idx++) {
    		if(left.ET < right.ET){
    			pfHadronicClusters[sector][12+idx].hwEta =
    					right.Eta - region[sector].hwEtaCenter;
    			pfHadronicClusters[sector][12+idx].hwPhi =
    					right.Phi - region[sector].hwPhiCenter + phi_offset_r;
    			pfHadronicClusters[sector][12+idx].hwPt =
    					right.ET;
    			pfHadronicClusters[sector][12+idx].hwPt >>= 2;

    			idx_right++;
    			right = caloClusters[(sector*3)+3][idx_right];
    		}
    		else {
    			pfHadronicClusters[sector][12+idx].hwEta =
						left.Eta - region[sector].hwEtaCenter;
				pfHadronicClusters[sector][12+idx].hwPhi =
						left.Phi - region[sector].hwPhiCenter + phi_offset_l;
				pfHadronicClusters[sector][12+idx].hwPt =
						left.ET;
				pfHadronicClusters[sector][12+idx].hwPt >>= 2;

    			idx_left++;
    			left = caloClusters[(sector*3)-1][idx_left];
    		}
    	}

    }

#ifndef __SYNTHESIS__
    for(loop sector=0; sector<N_SECTORS_PF; sector++) {
            std::cout<<"HADCALO | "<< sector<<"\n";
            for(int cluster =0 ; cluster<16 ; cluster++)
            {
				std::cout<<pfHadronicClusters[sector][cluster].hwPt<<","
						 <<pfHadronicClusters[sector][cluster].hwEta<<","
						 <<pfHadronicClusters[sector][cluster].hwPhi<<"\n";
            }
            std::cout<<"\n";
        }
#endif


    selectEGClusters(EGClusters,egCandidates);
    

    doPFClustringChain_5_1:
    for(loop i=0; i<N_SECTORS_PF; i++) {
        fwdlinpuppi(region[i], pfHadronicClusters[i], pfselne[i]);
    }



    ap_fixed<32,16> Ex_temp = 0;
    ap_fixed<32,16> Ey_temp = 0;

    pfExy(pfselne, Ex_temp, Ey_temp);

	doPFClustringChain_6_1:
    for(int i=0;i < 8 ; i++){
            egClusters[i]=egCandidates[i];
    }
    
    doPFClustringChain_7_1:
    for(int i=0; i < N_SECTORS_PF; i++ ) {
    	doPFClustringChain_7_2:
        for(int j =0 ; j< NNEUTRALS ; j++) {
            pfSelectedNutrals[i][j]=pfselne[i][j];
        }
    }

    Ex = Ex_temp;
    Ey = Ey_temp;

	#ifndef __SYNTHESIS__
	std::cout << "\nEx: " << Ex << "\nEy: " << Ey << std::endl;
	#endif

}






