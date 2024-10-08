
#include "algo_topIP1.h"

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

