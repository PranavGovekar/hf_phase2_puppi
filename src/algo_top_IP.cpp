#include "algo_topIP1.h"

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

    PFcluster egClusters[8]; // TODO . rmove hardcoding
    l1ct::PuppiObj pfSelectedNutrals[N_SECTORS_PF][NNEUTRALS];
    
    doPFClustringChain(link_in, egClusters , pfSelectedNutrals );


    #ifndef __SYNTHESIS__
           if(DEBUG_LEVEL > 0)
           {
                std::cout<<"@@EGClusters | " ;
                for(int i=0;i < 8 ; i++)
                {
                    std::cout<<egClusters[i].ET<<","<<egClusters[i].Eta<<","<<egClusters[i].Phi<<" | ";
                    //std::cout<<egClusters[i].data()<<" | ";
                }
                std::cout<<"\n";
          } 
    #endif

   
#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 0)
    {
        for(int i=0; i < N_SECTORS_PF; i++ )
        {
            std::cout<<"@@PUPPI | "<< i<<" | ";
            for(int j =0 ; j< NNEUTRALS ; j++)
            {
                std::cout<<pfSelectedNutrals[i][j].hwPt<<","<<pfSelectedNutrals[i][j].hwEta<<","<<pfSelectedNutrals[i][j].hwPhi<<","<<pfSelectedNutrals[i][j].hwPuppiW()<<" | ";
                //std::cout<<pfSelectedNutrals[i][j].data()<<" | ";
            }
            std::cout<<"\n";
        }
    }
#endif

 


    // PUPPI Funtionality

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

