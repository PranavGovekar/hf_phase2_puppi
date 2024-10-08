#include "pfClusterering.h"

// this function takes in 3 links and makes a region of 3 wedges (1 + 2 extra) i.e eta:(12+2) and phi:(4+8)
// this does not affect the seed searching because only the required area is searched. search region is eta:(12) and phi:(4+2)
void makeRegion(const ap_uint<LINK_WIDTH> link_in[LINKS_PER_REGION],
                hftower HFRegions[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR])
{

    #pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
    #pragma HLS ARRAY_PARTITION variable=HFRegions complete dim=0

    for(loop link=0; link<LINKS_PER_REGION; link++)
    {

        hftower HFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI/N_INPUT_LINKS];

        #pragma HLS ARRAY_PARTITION variable=HFTowers complete dim=0
        processInputLink(link_in[link], HFTowers);

        for(loop eta=0; eta<TOWERS_IN_ETA; eta++)
        {
            for(loop phi=0; phi<TOWERS_IN_PHI/N_INPUT_LINKS; phi++)
            {
                HFRegions[eta+1][(link*(TOWERS_IN_PHI/N_INPUT_LINKS))+phi] = HFTowers[eta][phi];
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
            std::cout << endl;
        }
        std::cout << endl;
    }
#endif
}

void findMaxEnergyTowerInPhi(const hftower HFRegion[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR],
                             const ap_uint<5> etaCenters[6], ap_uint<8>& phiC)
{
    #pragma HLS PIPELINE II=1
    hftower tempArray[8];
    ap_uint<4> index[8];

    #pragma HLS ARRAY_PARTITION variable=tempArray complete dim=0
    #pragma HLS ARRAY_PARTITION variable=index complete dim=0

    for(loop i=0; i<6; i++)
    {
        tempArray[i] = HFRegion[etaCenters[i]][i+3];
        index[i] = i;
    }

    for(loop i=8; i>1; i=(i/2))
    {
        for(loop j=0; j < i/2; j++)
        {
            if(tempArray[index[j*2]].energy < tempArray[index[(j*2) + 1]].energy)
            {
                index[j] = index[(j*2)+1];
            }
            else
            {
                index[j] = index[j*2];
            }
        }
    }

    phiC = index[0];
}


void findMaxEnergyTowerInEta(const hftower EtaTowers[TOWERS_IN_ETA], ap_uint<5>& etaC)
{
    #pragma HLS PIPELINE II=1
    hftower tempArray[16];
    ap_uint<5> index[16];

    #pragma HLS ARRAY_PARTITION variable=tempArray complete dim=0
    #pragma HLS ARRAY_PARTITION variable=index complete dim=0

    for(loop i=0; i<12; i++)
    {
        tempArray[i] = EtaTowers[i];
        index[i] = i;
    }


    for(loop i=16; i>1; i=(i/2))
    {
        for(loop j=0; j < i/2; j++)
        {
            if(tempArray[index[j*2]].energy < tempArray[index[(j*2) + 1]].energy)
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


void findMaxEnergyTower(const hftower HFRegion[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR],
                        ap_uint<5>& etaC,
                        ap_uint<8>& phiC)
{
    #pragma HLS ARRAY_PARTITION variable=HFRegion complete dim=0

    ap_uint<5> towersPhi[6];
    #pragma HLS ARRAY_PARTITION variable=towersPhi complete dim=0
    ap_uint<8> tempPhi;
    for(ap_uint<5> phi = 3; phi < 9; phi++)
    {
        hftower towersEta[12];
        #pragma HLS ARRAY_PARTITION variable=towersEta complete dim=0
        for(loop eta = 0; eta < 12; eta++)
        {
            towersEta[eta] = HFRegion[eta][phi];
        }
        findMaxEnergyTowerInEta(towersEta, towersPhi[phi-3]);
    }

    findMaxEnergyTowerInPhi(HFRegion, towersPhi, tempPhi);

    etaC = towersPhi[tempPhi];
    phiC = tempPhi + 3;
#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL >7 )
    {
        std::cout<<"  Center found at : "<<etaC<<" , "<<phiC<<"\n" ;
    }
#endif
}

// Function to form the cluster and zero out the tower energies
void formClusterAndZeroOut(hftower HFRegion[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR],
                           ap_uint<5> etaC,
                           ap_uint<8> phiC,
                           ap_uint<12>& etaSum)
{

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
void makeCaloClusters (const ap_uint<LINK_WIDTH> regionLinks[LINKS_PER_REGION],PFcluster Clusters[N_PF_CLUSTERS])
{

    #pragma HLS PIPELINE
    #pragma HLS ARRAY_PARTITION variable=Clusters complete dim=0

    hftower HFRegion[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR];
    #pragma HLS ARRAY_PARTITION variable=HFRegion complete dim=0

    makeRegion(regionLinks, HFRegion);

#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 7)
    {
        for(loop phiI=0; phiI<NTOWER_IN_PHI_PER_SECTOR; phiI++)
        {
            std::cout<<std::setw(3)<<phiI;
        }
        std::cout<<"  <--Phi"<<"\n";
        for(loop etaI=0; etaI<NTOWER_IN_ETA_PER_SECTOR; etaI++)
        {
            for(loop phiI=0; phiI<NTOWER_IN_PHI_PER_SECTOR; phiI++)
            {
                std::cout<<std::setw(5)<< HFRegion[etaI][phiI].energy<<" | ";
            }
            std::cout<<"  <--"<<etaI<<"\n";
        }
        std::cout<<endl;
    }
#endif

    for(loop cluster = 0; cluster < N_PF_CLUSTERS; cluster++)
    {
        ap_uint<8> etmax = 0;
        ap_uint<8> phiC = 1;
        ap_uint<5> etaC = 1;
        ap_uint<12> etaSum = 0;

        findMaxEnergyTower(HFRegion, etaC, phiC);

        if(HFRegion[etaC][phiC].energy > MIN_CLUSTER_SEED_ENERGY)
        {
            formClusterAndZeroOut(HFRegion, etaC, phiC, etaSum);

            if(phiC > 3 && phiC < 8)
            {
                Clusters[cluster].Eta = etaC;
                Clusters[cluster].Phi = phiC;
                Clusters[cluster].ET = etaSum;
            }
        }
    }

#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 2)
    {
        for(loop cluster=0; cluster<N_PF_CLUSTERS; cluster++)
        {
            std::cout<<"Cluster " << cluster << " E = "<< Clusters[cluster].ET
                     <<", center = ("<< Clusters[cluster].Eta <<","<< Clusters[cluster].Phi << ")\n";
        }
        std::cout<<endl;
    }
#endif
}

void packer(PFcluster Clusters[N_PF_CLUSTERS], const ap_uint<576>& link_out_sector, const ap_uint<7> sector)
{

    #pragma HLS ARRAY_PARTITION variable=Clusters complete dim=0

    ap_uint<10> start = 0;
    ap_uint<10> end = start + 63;
    PFcluster zero;
packer_loop:
    for (int cluster = 0; cluster < 8; cluster++)
    {
        #pragma HLS UNROLL
        if (Clusters[cluster].ET > 0)
        {
            link_out_sector.range(end, start) = Clusters[cluster].data();
        }
        else
        {
            link_out_sector.range(end, start) = 0 ;//zero.data();
        }
        start += 64;
        end = start + 63;
    }
}

void makePFClusters( const ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], PFcluster ClustersOut[N_OUTPUT_LINKS][N_SORT_ELEMENTS]    )
{
    // TODO : FOR LOOP TO DO THIS PLEASE !! :)
    ap_uint<LINK_WIDTH> linksInSector[N_SECTORS][3] ;

    linksInSector[0][0] = link_in[17];
    linksInSector[0][1] = link_in[ 0];
    linksInSector[0][2] = link_in[ 1];
    for(int i=1; i<N_SECTORS-1; i++)
    {
        linksInSector[i][0] = link_in[i-1];
        linksInSector[i][1] = link_in[i  ];
        linksInSector[i][2] = link_in[i+1];
    }
    linksInSector[N_SECTORS -1][0] = link_in[N_SECTORS -2 ];
    linksInSector[N_SECTORS -1][1] = link_in[N_SECTORS -1];
    linksInSector[N_SECTORS -1][2] = link_in[ 0];

    // Too many '#define ed' elements ?
    PFcluster caloClusters[N_OUTPUT_LINKS][N_SORT_ELEMENTS];
    #pragma HLS ARRAY_PARTITION variable=Clusters complete dim=0

link_unroll_loop:
    for(loop link=0; link<N_OUTPUT_LINKS; link++)
    {
sector_unroll_loop:
        for(loop sector=0; sector<LINKS_PER_REGION; sector++)
        {
            PFcluster tempClusters[N_PF_CLUSTERS];
            #pragma HLS ARRAY_PARTITION variable=tempClusters complete dim=0

            makeCaloClusters(linksInSector[(link*LINKS_PER_REGION) + sector], tempClusters);

            for(loop i=0; i<N_PF_CLUSTERS; i++)
            {
                if(tempClusters[i].ET > 0 )
                {
                    caloClusters[link][(sector*N_PF_CLUSTERS)+i].ET  = tempClusters[i].ET;
                    caloClusters[link][(sector*N_PF_CLUSTERS)+i].Eta = tempClusters[i].Eta - EXTRA_IN_ETA;
                    caloClusters[link][(sector*N_PF_CLUSTERS)+i].Phi = tempClusters[i].Phi - 4  + ( link*3 + sector )*4  ; // NEED TO ADD LOGIC AS EXPLANATION
                }
            }
        }

#ifndef __SYNTHESIS__
        if (DEBUG_LEVEL > 5)
        {
            std::cout<< "Before Sorting : " << endl;
            for(loop cluster=0; cluster<16; cluster++)
            {
                std::cout<<"Cluster " << cluster << " E = "<< caloClusters[link][cluster].ET
                         <<", center = ("<< caloClusters[link][cluster].Eta <<","<< caloClusters[link][cluster].Phi << ")\n";
            }
            std::cout<<endl;
        }
#endif

        PFcluster sortedClusters[N_SORT_ELEMENTS];
        #pragma HLS ARRAY_PARTITION variable=sortedClusters complete dim=0
        bitonicSort16(caloClusters[link], sortedClusters);

#ifndef __SYNTHESIS__
        if (DEBUG_LEVEL > 2)
        {
            std::cout<< "Calo clusters in out link :  "<< link << endl;
            for(loop cluster=0; cluster<16; cluster++)
            {
                std::cout<<"  Cluster " << cluster << " E = "<< sortedClusters[cluster].ET
                         <<", center = ("<< sortedClusters[cluster].Eta <<","<< sortedClusters[cluster].Phi << ")\n";
            }
            std::cout<<endl;
        }
#endif
        for(int i =0 ; i< N_SORT_ELEMENTS ; i++)
        {
            ClustersOut[link][i]= sortedClusters[i];
        }
    }
}







