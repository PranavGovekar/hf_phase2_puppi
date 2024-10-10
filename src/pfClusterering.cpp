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
                if(etaSum = 2*HFRegion[etaC][phiC].energy  )
                {
                    Clusters[cluster].isEG= 1;
                }
                else {
                    Clusters[cluster].isEG= 0;
                }
            }
        }
    }

#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 7)
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
    l1ct::PuppiObj pfselne[N_SECTORS_PF][NNEUTRALS],
    l1ct::HadCaloObj H_in_regionized[N_SECTORS_PF][NCALO],
    const l1ct::PFRegion region[N_SECTORS_PF]
)
{
    # pragma HLS INLINE off
puppi:
    for(loop i=0; i<N_SECTORS_PF; i++)
    {
        fwdlinpuppi(region[i], H_in_regionized[i], pfselne[i]);
    }
}


template <int SIZE, int N_POWER_2> void getTopItemIndexInArray(const ap_uint<12> arrayIn[SIZE],ap_uint<8> &index)
{
    ap_uint<12> array[N_POWER_2];
    ap_uint<8> indices[N_POWER_2];
    for(uint8_t i=0;i < SIZE;i++)
    {
        array[i]=arrayIn[i];
        indices[i]=i;
    }   
    for(uint8_t i=SIZE;i < N_POWER_2;i++)
    {
        array[i]=0;
        indices[i]=SIZE-1;
    }   

    for(loop i=N_POWER_2; i>1; i=(i/2))
    {
        for(uint8_t j=0; j < i/2; j++)
        {
            
            if(array[indices[j]] < array[indices[j + i/2]])
            {
                 indices[j] = indices[j+i/2];
            }
           // else {
           //       indices[j] = indices[j];
           //  }
        }
    }
    index=indices[0]; 
}
void selectEGClusters(const PFcluster caloClusters[N_INPUT_LINKS][N_SORT_ELEMENTS] ,PFcluster egClusters[8] ) // TODO , remove hardcoded
{

    PFcluster  egCandidateClusters[N_SECTORS_PF*12/*N_SORT_ELEMENTS*LINKS_PER_REGION*/];
    ap_uint<12> egCandidateEnergies[N_SECTORS_PF*12/*N_SORT_ELEMENTS*LINKS_PER_REGION*/];
    
    PFcluster  egSelectedClusters[8];
    
    /*
    PFcluster inputClusters[N_INPUT_LINKS][N_SORT_ELEMENTS]
    for(int i =0 ; i < N_INPUT_LINKS ; i++) {
        for(int j =0 ; j < N_SORT_ELEMENTS ; j++) {
            inputClusters[i][j]=caloClusters[i][j];
        }
    }
    */

    // TODO .: REMOVE THE HARDCODING !
    

    for(int sector=0 ; sector < 6 ; sector++)
    {
        for(int wedge=0 ; wedge < 3 ; wedge++)
        {
            // can take top 3 per wedge too ! reducing the complexy of memory for 128 --> 64
            for(int i =0 ; i < 4 ; i++) {        
                int idx = (sector*3 + wedge)*4 +i ;
                //egSelectsClusters[idx] = inputClusters[sector*3 + wedge][i] ;
                egCandidateClusters[idx]   = caloClusters[sector*3 + wedge][i] ;
                egCandidateEnergies[idx]   = caloClusters[sector*3 + wedge][i].ET ;
                // std::cout<<" idx : "<<int(idx)<<" [ "<<sector<<","<<wedge<<" ] : "<<egCandidateEnergies[idx]<<"\n";

      //          if( ! egCandidateClusters[idx].isEG ) egCandidateEnergies[idx]=0;
            }
        }
    }
    
    ap_uint<8> select_index;
    for(int i=0 ; i < 8 ; i++)
    {
        getTopItemIndexInArray<N_SECTORS_PF*12,128>(egCandidateEnergies,select_index);
        //std::cout<<"  eg --> "<<i<<" : "<<  egCandidateEnergies[select_index]<<" ["<<select_index <<"] \n";
        egSelectedClusters[i] = egCandidateClusters[select_index]; 
        egCandidateEnergies[select_index]=0;
    }

    for(int i=0 ; i < 8 ; i++)
    {
        egClusters[i] = egSelectedClusters[i];
    }


}



void doPFClustringChain( const ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], PFcluster egClusters[8] , l1ct::PuppiObj pfSelectedNutrals[N_SECTORS_PF][NNEUTRALS]   )
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
    // PFcluster caloClusters[N_OUTPUT_LINKS][N_SORT_ELEMENTS];
    PFcluster caloClusters[N_INPUT_LINKS][N_SORT_ELEMENTS];
    PFcluster egCandidates[8]; // TODO Remive HARDCODING
    #pragma HLS ARRAY_PARTITION variable=Clusters complete dim=0

    l1ct::PuppiObj pfselne[N_SECTORS_PF][NNEUTRALS];
    l1ct::HadCaloObj pfHadronicClusters[N_SECTORS_PF][NCALO];

    // define the  sector boundaries and overlaps
    l1ct::PFRegion region[N_SECTORS_PF];
    regions_init:
    for(int i=0 ; i < N_SECTORS_PF ; i++)
        {
            region[i].hwEtaCenter = l1ct::glbeta_t(6);
            region[i].hwPhiCenter = l1ct::glbphi_t(12*i+5);

            region[i].hwEtaExtra = l1ct::eta_t(1);
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
    
sector_unroll_loop:
    for(loop sector=0; sector<N_SECTORS_PF; sector++)
    {
wedgesPreSector_unroll_loop:
        for(loop wedge=0; wedge<LINKS_PER_REGION; wedge++)
        {
            PFcluster tempClusters[N_PF_CLUSTERS];
            
            #pragma HLS ARRAY_PARTITION variable=tempClusters complete dim=0
            makeCaloClusters(linksInSector[(sector*LINKS_PER_REGION) + wedge], tempClusters);
            

            for(loop i=0; i<N_PF_CLUSTERS; i++)
            {
                 // PROPER ASSIGNMENT
                 // Need to set the  Etas and Phis Realtive
                 // Need to make sure the Eta and Phi has proper decimal(uint) representation
                 // Need to be carefull of the data type conversions
                  caloClusters[sector*LINKS_PER_REGION+wedge][i] = tempClusters[i] ;
                  
                //pfHadronicClusters[sector][N_PF_CLUSTERS*wedge + i ] = tempClusters[i] ;
                  pfHadronicClusters[sector][N_PF_CLUSTERS*wedge + i ].hwPt  = tempClusters[i].ET;
                  pfHadronicClusters[sector][N_PF_CLUSTERS*wedge + i ].hwEta = tempClusters[i].Eta - EXTRA_IN_ETA;
                  pfHadronicClusters[sector][N_PF_CLUSTERS*wedge + i ].hwPhi = tempClusters[i].Phi - 4  + ( sector*3 + wedge )*4 ;
                  
                  //caloClusters[link][(sector*N_PF_CLUSTERS)+i].ET  = tempClusters[i].ET;
                  //caloClusters[link][(sector*N_PF_CLUSTERS)+i].Eta = tempClusters[i].Eta - EXTRA_IN_ETA;
                  //caloClusters[link][(sector*N_PF_CLUSTERS)+i].Phi = tempClusters[i].Phi - 4  + ( link*3 + sector )*4  ; 


            #ifndef __SYNTHESIS__
                if(DEBUG_LEVEL > 8)
                {
                  if(tempClusters[i].Eta > 0 )
                  std::cout<<" for sector : "<<sector<<", wedge : "<<wedge
                           <<" et/eta/phi :: "
                           <<"calo : "<<tempClusters[i].ET<<" / "<<tempClusters[i].Eta<<" / "<<tempClusters[i].Phi<<" | "
                           <<"hadc : "<<pfHadronicClusters[sector][N_PF_CLUSTERS*wedge + i].hwPt<<" / "
                                      <<pfHadronicClusters[sector][N_PF_CLUSTERS*wedge + i].hwEta<<" / "
                                      <<pfHadronicClusters[sector][N_PF_CLUSTERS*wedge + i].hwPhi<<"  ";
               }   
            #endif
                  pfHadronicClusters[sector][N_PF_CLUSTERS*wedge + i ].hwEta = pfHadronicClusters[sector][N_PF_CLUSTERS*wedge + i ].hwEta-region[sector].hwEtaCenter;
                  pfHadronicClusters[sector][N_PF_CLUSTERS*wedge + i ].hwPhi = pfHadronicClusters[sector][N_PF_CLUSTERS*wedge + i ].hwPhi-region[sector].hwPhiCenter;
                
            #ifndef __SYNTHESIS__
                if(DEBUG_LEVEL > 8)
                {
                  if(tempClusters[i].Eta > 0 )
                  std::cout<<"hadc offset : "<<pfHadronicClusters[sector][N_PF_CLUSTERS*wedge + i].hwPt<<" / "
                                      <<pfHadronicClusters[sector][N_PF_CLUSTERS*wedge + i].hwEta<<" / "
                                      <<pfHadronicClusters[sector][N_PF_CLUSTERS*wedge + i].hwPhi<<"  "
                            <<std::endl;
                        
               } 
            #endif
            }   
       }
    }
#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 0)
    {
    for(loop sector=0; sector<N_SECTORS_PF; sector++)
        for(int i=0; i < N_OUTPUT_LINKS; i++ )
        {
            std::cout<<"@@CALOCLUSTER | "<< i<<" | ";
            for(int wedge =0 ; wedge< 3 ; wedge++)
            {
                for(int i =0 ; i <4 ;i++)
                  std::cout<<caloClusters[sector*LINKS_PER_REGION+wedge][i].data()<<" | " ;
            }
            std::cout<<"\n";
        }
    }
#endif


   selectEGClusters(caloClusters,egCandidates);
    
    for(loop sector=0; sector<N_SECTORS_PF; sector++)
    {
        int8_t left_sector  =   sector-1;
        int8_t right_sector =   sector+1;
        if(left_sector  <  0 )         left_sector  = N_SECTORS-1;
        if(right_sector == N_SECTORS ) right_sector = 0;
        
        uint8_t hadcalo_offset  = N_PF_CLUSTERS*(LINKS_PER_REGION);
        
        uint8_t left_offset  = N_PF_CLUSTERS*(LINKS_PER_REGION-1);
        uint8_t right_offset = 0;
        
        for( loop i=0 ; i < 4 ; i++)
        {
            if(pfHadronicClusters[left_offset][left_offset].hwPt < pfHadronicClusters[right_sector][right_offset].hwPt )
            {   
                // UPDATE THE OFFSETS HERE
                
                pfHadronicClusters[sector][hadcalo_offset]    = pfHadronicClusters[right_sector][right_offset];
                pfHadronicClusters[sector][hadcalo_offset].hwPhi +=12;
                right_offset++;
            }
            else 
            {
                // UPDATE THE OFFSETS HERE
                
                pfHadronicClusters[sector][hadcalo_offset] = pfHadronicClusters[left_sector][left_offset];
                pfHadronicClusters[sector][hadcalo_offset].hwPhi -=12;
                left_offset++;

            }
            hadcalo_offset++;
        }
    }

#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 5)
    {
        for(int i=0; i < N_SECTORS_PF; i++ )
        {
            std::cout<<"@@HADCALO | "<< i<<" | ";
            for(int j =0 ; j< NCALO ; j++)
            {
              if( pfHadronicClusters[i][j].hwPt > 0  )
                std::cout<<" [ "<<region[i].isInside(pfHadronicClusters[i][j])<<" ] "
                         <<pfHadronicClusters[i][j].hwPt<<","
                         <<pfHadronicClusters[i][j].hwEta<<","
                         <<pfHadronicClusters[i][j].hwPhi<<" | ";
            }
            std::cout<<"\n";
        }
    }
#endif

    
    compute(pfselne, pfHadronicClusters, region);

    for(int i=0;i < 8 ; i++)
    {
        egClusters[i]=egCandidates[i];
    }
    
    for(int i=0; i < N_SECTORS_PF; i++ )
    {
        for(int j =0 ; j< NNEUTRALS ; j++)
        {
            pfSelectedNutrals[i][j]=pfselne[i][j];
        }
    }

}







