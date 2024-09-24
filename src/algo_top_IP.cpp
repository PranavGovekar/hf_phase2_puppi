#include "algo_topIP1.h"

#define DEBUG_LEVEL 1
// Unpacker function for a wedge
// This function takes one link and makes 12x4 grid of hftower Array
void processInputLink( ap_uint<LINK_WIDTH> link_in,
                       hftower towerGrid[TOWERS_IN_ETA][TOWERS_IN_PHI/N_INPUT_LINKS])
{

    #pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
    #pragma HLS ARRAY_PARTITION variable=towerGrid complete dim=0

    for(loop j=0; j<(TOWERS_IN_PHI/N_INPUT_LINKS)/2; j=j+2)
    {
        for(loop i=0; i<TOWERS_IN_ETA-1; i++)
        {
            hftower halfEnergy = hftower(link_in.range((i*10)+9, i*10));
            halfEnergy.energy = halfEnergy.energy >> 1;

            towerGrid[i][j] = halfEnergy;
            towerGrid[i][j+1] = halfEnergy;

            hftower halfEnergy_0 = hftower(link_in.range((i*10+110)+9, i*10+110));
            halfEnergy_0.energy = halfEnergy_0.energy >> 1;

            towerGrid[i][j+2] = halfEnergy_0;
            towerGrid[i][j+3] = halfEnergy_0;
        }
    }

    hftower A10 = towerGrid[TOWERS_IN_ETA-2][0] ;
    hftower B10 = towerGrid[TOWERS_IN_ETA-2][2] ;

    A10.energy = A10.energy >> 1 ;
    B10.energy = B10.energy >> 1 ;

    towerGrid[TOWERS_IN_ETA-2][0] = A10;
    towerGrid[TOWERS_IN_ETA-2][1] = A10;
    towerGrid[TOWERS_IN_ETA-2][2] = A10;
    towerGrid[TOWERS_IN_ETA-2][3] = A10;

    towerGrid[TOWERS_IN_ETA-1][0] = B10;
    towerGrid[TOWERS_IN_ETA-1][1] = B10;
    towerGrid[TOWERS_IN_ETA-1][2] = B10;
    towerGrid[TOWERS_IN_ETA-1][3] = B10;


#ifndef __SYNTHESIS__

    if (DEBUG_LEVEL > 10)
    {
        for(loop eta=0; eta<TOWERS_IN_ETA; eta++)
        {
            for(loop phi=0; phi<TOWERS_IN_PHI/N_INPUT_LINKS; phi++)
            {
                std::cout<< "\t" << towerGrid[eta][phi].energy;
            }
            std::cout << endl;
        }
        std::cout << endl;
    }
#endif
}

// Unpacker function
// This function  unpacks the input set of links into fixed-grid Super Towers

void unpackToSuperTowers(const ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS],
                         hftower superTowers[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2])
{

    hftower HFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI];

    #pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
    #pragma HLS ARRAY_PARTITION variable=HFTowers complete dim=0

    for(loop link=0; link<N_INPUT_LINKS; link++)
    {

        hftower HFTowers_temp[TOWERS_IN_ETA][TOWERS_IN_PHI/N_INPUT_LINKS];
        #pragma HLS ARRAY_PARTITION variable=HFTowers complete dim=0
        processInputLink(link_in[link], HFTowers_temp);

        for(loop eta=0; eta<TOWERS_IN_ETA; eta++)
        {
            for(loop phi=0; phi<TOWERS_IN_PHI/N_INPUT_LINKS; phi++)
            {
                HFTowers[eta][(link*(TOWERS_IN_PHI/N_INPUT_LINKS))+phi] = HFTowers_temp[eta][phi];
            }
        }
    }

    for (loop i = 0; i < 4; ++i)
    {
        for (loop j = 0; j < 24; ++j)
        {

            for (loop m = 0; m < 3; ++m)
            {
                for (loop n = 0; n < 3; ++n)
                {
                    superTowers[i+1][j+1].energy += HFTowers[i * 3 + m][j * 3 + n].energy;
                }
            }
        }
    }

#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 5)
    {
        std::cout<<" In  makeSuperTower \n"<<"HF Towers";
        for(loop eta=0; eta<TOWERS_IN_ETA; eta++)
        {
            for(loop phi=0; phi<TOWERS_IN_PHI; phi++)
            {
                std::cout<< "\t [ "<<eta<<","<<phi <<" ] : "<< HFTowers[eta][phi].energy<<" | ";
            }
            std::cout << endl;
        }
        std::cout << endl;
    }
    if(DEBUG_LEVEL > 2)
    {
        std::cout<<"HF Super Towers : "<<"\n";
        for(loop eta=0; eta<(TOWERS_IN_ETA/3) + 2; eta++)
        {
            for(loop phi=0; phi<(TOWERS_IN_PHI/3) + 2; phi++)
            {
                std::cout<< "["<<eta<<","<<phi<<"] :" << superTowers[eta][phi].energy;
            }
            std::cout << endl;
        }
        std::cout << endl;
    }

    if(DEBUG_LEVEL > 0)
    {
        for(loop eta=0; eta<(TOWERS_IN_ETA/3) + 2; eta++)
        {
            std::cout<<"@@HFSuperTowers | "<<eta<<" | ";
            for(loop phi=0; phi<(TOWERS_IN_PHI/3) + 2; phi++)
            {
                std::cout<<setw(3)<<superTowers[eta][phi].energy<<" |";
            }
            std::cout << endl;
        }
    }
#endif
}


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





template <int SIZE, int N_POWER_2>
void findMaxEnergyTowerInArray(const hftower Towers[SIZE], ap_uint<10>& center)
{
    #pragma HLS PIPELINE
//
//    // Calculate the next power of two
//    int N_POWER_2 = SIZE;
//    N_POWER_2 -= 1;
//    N_POWER_2 |= N_POWER_2 >> 1;
//    N_POWER_2 |= N_POWER_2 >> 2;
//    N_POWER_2 |= N_POWER_2 >> 4;
//    N_POWER_2 += 1;

    hftower tempArray[N_POWER_2];
    ap_int<10> index[N_POWER_2];

    #pragma HLS ARRAY_PARTITION variable=tempArray complete dim=0
    #pragma HLS ARRAY_PARTITION variable=index complete dim=0

    for(ap_uint<10> i=0; i<SIZE; i++)
    {
        tempArray[i] = Towers[i];
        index[i] = i;
    }
    for(ap_uint<10> i=SIZE; i<N_POWER_2; i++)
    {
        index[i] = SIZE;
    }


    for(ap_uint<10> i=N_POWER_2; i>1; i=(i/2))
    {
        //std::cout<<" -- > for i = "<<i<<"\n";
        //for(ap_uint<10> j=0; j < i; j++)
        //    std::cout<<setw(3)<<tempArray[index[j]].energy<<" , ";
        //std::cout<<std::endl;
        //for(ap_uint<10> j=0; j < i; j++)
        //    std::cout<<setw(3)<<index[j]<<" , ";
        //std::cout<<std::endl;
        for(ap_uint<10> j=0; j < i/2; j++)
        {
            //if(tempArray[index[j*2]].energy < tempArray[index[(j*2) + 1]].energy)
            if(tempArray[index[j]].energy < tempArray[index[j + i/2]].energy)
            {
                index[j] = index[j+i/2];
            }
            else
            {
                index[j] = index[j];
            }
        }
    }

    center = index[0];
}


void findMaxEnergySuperTower(const hftower HFRegion[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2],
                             ap_uint<10>& etaC,
                             ap_uint<10>& phiC)
{

    #pragma HLS PIPELINE
    #pragma HLS ARRAY_PARTITION variable=HFRegion complete dim=0

    hftower towersPhi[4];
    ap_uint<10> PhiCenters[4];
    ap_uint<10> EtaCenter;

    #pragma HLS ARRAY_PARTITION variable=towersPhi complete dim=0
    #pragma HLS ARRAY_PARTITION variable=PhiCenters complete dim=0

    for(ap_uint<5> eta = 1; eta < 5; eta++)
    {
        //std::cout<<"   ETA SCAN !! "<<eta<<" \n";
        findMaxEnergyTowerInArray<24,32>(&HFRegion[eta][1], PhiCenters[eta-1]);
        PhiCenters[eta-1] += 1;
        towersPhi[eta-1] = HFRegion[eta][PhiCenters[eta-1]];
        //std::cout<<" eta = "<<eta<<" | ";
        //for(int i =0 ; i < (TOWERS_IN_PHI/3)+2  ;  i++)
        //{
        //    std::cout<<HFRegion[eta][i].energy<<",";
        //}
        //std::cout<<"   | [ "<<towersPhi[PhiCenters[eta-1]].energy<<" / "<< PhiCenters[eta-1] <<" ] \n";
    }

    findMaxEnergyTowerInArray<4,4>(towersPhi, EtaCenter);
    //std::cout<<"   Phi SCAN !! "<< EtaCenter <<"\n";

    phiC = PhiCenters[EtaCenter];
    etaC = EtaCenter+1;
}

void formJetsAndZeroOut(hftower superTowers[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2],
                        ap_uint<10> etaC,
                        ap_uint<10> phiC,
                        ap_uint<12>& etaSum)
{

    #pragma HLS ARRAY_PARTITION variable=superTowers complete dim=0

    etaSum = superTowers[etaC][phiC].energy +
             superTowers[etaC+1][phiC].energy +
             superTowers[etaC-1][phiC].energy +
             superTowers[etaC][phiC+1].energy +
             superTowers[etaC][phiC-1].energy +
             superTowers[etaC+1][phiC+1].energy +
             superTowers[etaC-1][phiC+1].energy +
             superTowers[etaC+1][phiC-1].energy +
             superTowers[etaC-1][phiC-1].energy;

    superTowers[etaC][phiC].energy = 0;
    superTowers[etaC+1][phiC].energy = 0;
    superTowers[etaC-1][phiC].energy = 0;
    superTowers[etaC][phiC+1].energy = 0;
    superTowers[etaC][phiC-1].energy = 0;
    superTowers[etaC+1][phiC+1].energy = 0;
    superTowers[etaC-1][phiC+1].energy = 0;
    superTowers[etaC+1][phiC-1].energy = 0;
    superTowers[etaC-1][phiC-1].energy = 0;
}

void selectTaus(const jets Jet[9], jets Taus[9])
{

    ap_uint<4> count = 0;
    for(loop idx = 0; idx < 9; idx++)
    {
        // todo : fix to multiplication
        if(((float)Jet[idx].seedET)/((float)Jet[idx].ET) >= 0.7f)
        {
            Taus[count] = Jet[idx];
            count++;
        }
    }

#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 2)
    {
        std::cout<< "Taus : " << endl;
        for(loop cluster=0; cluster<9; cluster++)
        {
            std::cout<<"Tau " << cluster << " E = "<< Taus[cluster].ET << " Seed = " << Taus[cluster].seedET
                     <<", center = ("<< Taus[cluster].Eta <<","<< Taus[cluster].Phi << ")\n";
        }
        std::cout<<endl;
    }
#endif

}

void makeJets(hftower superTowers[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2],jets Jet[9])
{

    for(loop idx = 0; idx < 9; idx++)
    {
        ap_uint<10> phiC = 1;
        ap_uint<10> etaC = 1;
        ap_uint<12> etSum = 0;

        findMaxEnergySuperTower(superTowers, etaC, phiC);

        Jet[idx].seedET = superTowers[etaC][phiC].energy;
        formJetsAndZeroOut(superTowers, etaC, phiC, etSum);

        if(etSum ==0 )
        {
            etaC=0 ;
            phiC=0;
        }
        Jet[idx].ET  = etSum;
        Jet[idx].Eta = etaC;
        Jet[idx].Phi = phiC;

    }

#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 2)
    {
        std::cout<< "Jets : " << endl;
        for(loop cluster=0; cluster<9; cluster++)
        {
            std::cout<<"Jet " << cluster << " E = "<< Jet[cluster].ET << " Seed = " << Jet[cluster].seedET
                     <<", center = ("<< Jet[cluster].Eta <<","<< Jet[cluster].Phi << ")\n";
        }
        std::cout<<endl;
    }
#endif
}

void makeCaloObjects(const ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS],
                     jets Jets[9], jets Taus[9] )
{

    #pragma HLS PIPELINE

    hftower superTowerGrid[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2];
    jets _Jets[9];
    jets _Taus[9];

    unpackToSuperTowers(link_in, superTowerGrid);

    makeJets(superTowerGrid, _Jets);

    selectTaus(_Jets, _Taus);

    for(loop idx = 0; idx < 9; idx++)
    {
        Jets[idx] = _Jets[idx];
    }

    for(loop idx = 0; idx < 9; idx++)
    {
        Taus[idx] = _Taus[idx];
    }

}
//
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

void compute(
    l1ct::PuppiObj pfselne[N_SECTORS][NNEUTRALS],
    l1ct::HadCaloObj H_in_regionized[N_SECTORS][NCALO],
    const l1ct::PFRegion region[N_SECTORS]
)
{
    # pragma HLS INLINE off
puppi:
    for(loop i=0; i<N_SECTORS; i++)
    {
        fwdlinpuppi(region[i], H_in_regionized[i], pfselne[i]);
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


void algo_topIP1(ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], ap_uint<576> link_out[N_OUTPUT_LINKS])
{

    #pragma HLS PIPELINE II=9
    #pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
    #pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
    #pragma HLS INTERFACE ap_ctrl_hs port=return

    jets Jets[9];
    jets Taus[9];
    #pragma HLS ARRAY_PARTITION variable=Jets complete dim=0
    #pragma HLS ARRAY_PARTITION variable=link_in complete dim=0


#ifndef __SYNTHESIS__

    if(DEBUG_LEVEL > 0)
    {
        std::cout<<"# @@HFTowers | phi | tower_et_phi...\n";
        std::cout<<"# @@HFSuperTowers | eta | tower_et_phi...\n";
        std::cout<<"# @@Jets | nJet | jet.data()...\n";
        std::cout<<"# @@Taus | nJet | jet.data()...\n";
        std::cout<<"# @@CALOCLUSTER | link | pf_cluster.data() ...\n";
    }

    if(DEBUG_LEVEL > 0)
    {
        hftower towerGrid[TOWERS_IN_ETA][TOWERS_IN_PHI/N_INPUT_LINKS];
        for(int link = 0 ; link < N_INPUT_LINKS ; link++)
        {
            processInputLink(link_in[link], towerGrid);
            for(loop phi=0; phi<TOWERS_IN_PHI/N_INPUT_LINKS; phi++)
            {
                std::cout<<"@@HFTowers | "<<link<<","<<phi<<" | ";
                for(loop eta=0; eta<TOWERS_IN_ETA; eta++)
                {
                    std::cout<<setw(3)<<towerGrid[eta][phi].energy<<" | ";
                }
                std::cout << endl;
            }
        }
    }
#endif

    makeCaloObjects(link_in, Jets,Taus);
#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 0)
    {
        std::cout<<"@@JETS | 9 | ";
        for(int i=0 ; i< 9 ; i++)
        {
            std::cout<<Jets[i].data()<<" | ";
        }
        std::cout<<"\n";

        std::cout<<"@@Taus | 9 | ";
        for(int i=0 ; i< 9 ; i++)
        {
            std::cout<<Taus[i].data()<<" | ";
        }
        std::cout<<"\n";
    }
#endif

    PFcluster caloClusters[N_OUTPUT_LINKS][N_SORT_ELEMENTS];
    makePFClusters(link_in,caloClusters );

#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 0)
    {
        for(int i=0; i < N_OUTPUT_LINKS; i++ )
        {
            std::cout<<"@@CALOCLUSTER | "<< i<<" | ";
            for(int j =0 ; j< N_SORT_ELEMENTS ; j++)
            {
                std::cout<<caloClusters[i][j].data()<<" | ";
            }
            std::cout<<"\n";
        }
    }
#endif


    // PUPPI Funtionality
    l1ct::PuppiObj pfselne[N_SECTORS][NNEUTRALS];
    l1ct::HadCaloObj H_in_regionized[N_SECTORS][NCALO];

    #pragma HLS ARRAY_PARTITION variable=H_in_regionized complete dim=0
    #pragma HLS ARRAY_PARTITION variable=pfselne complete dim=0

    ap_uint<12> start, end;

    // define the  sector boundaries and overlaps
    l1ct::PFRegion region[N_SECTORS];
regions_init:
    for(int i=0 ; i < N_SECTORS ; i++)
    {
        region[i].hwEtaCenter = l1ct::glbeta_t(12);
        region[i].hwEtaHalfWidth = l1ct::eta_t(12);
        region[i].hwEtaExtra = l1ct::eta_t(0);
        region[i].hwPhiExtra = l1ct::phi_t(2);
        region[i].hwPhiHalfWidth = l1ct::phi_t(6);
        region[i].hwPhiCenter = l1ct::glbphi_t(12*i+5);
    }



    //initialize the neutral hadron objects
    //hadcalo_init:
    //    for( loop  i= 0; i < N_SECTORS; i++) {
    //    #pragma HLS unroll
    //        for( loop j = 0; j < NCALO ; j++) {
    //        #pragma HLS unroll
    //            H_in_regionized[i][j].clear();
    //        }
    //    }
    //unpack and regionize the PF Clusters to the neutral hadron objects

    //unpackLinksToHadCalo(link_in, H_in_regionized, region);

    // invoke puppi for the regionized clusters

    //compute(pfselne, H_in_regionized, region);


    // pack the puppi outputs into links
    //packLinks(pfselne, link_out);





}

