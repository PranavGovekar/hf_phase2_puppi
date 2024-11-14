#include "algo_topIP1.h"

void packPuppi(	l1ct::PuppiObj pfselne[NNEUTRALS],
			ap_uint<576> &link_out){
#pragma HLS INLINE
	const ap_uint<8> BW = 64;
	ap_uint<12> start = 0;
	ap_uint<576> temp_link;

    for(loop idx=0; idx<NNEUTRALS; idx++) {
    	temp_link.range(start+(BW-1),start) = pfselne[idx].pack();

        start = start + BW;
    }

    link_out = temp_link;
}

void packEG(l1ct::HadCaloObj pfselne[NNEUTRALS],
			ap_uint<576> &link_out){
	const ap_uint<8> BW = 64;
	ap_uint<12> start = 0;
	ap_uint<576> temp_link;

    for(loop idx=0; idx<NNEUTRALS; idx++) {
    	temp_link.range(start+(BW-1),start) = pfselne[idx].pack();

        start = start + BW;
    }

    link_out = temp_link;
}

void packJets(jets Jets[9], ap_uint<576> &link_out){
	const ap_uint<8> BW = 64;
	ap_uint<12> start = 0;
	ap_uint<576> temp_link;

    for(loop idx=0; idx<6; idx++) {
    	temp_link.range(start+(BW-1),start) = Jets[idx].data();

        start = start + BW;
    }

    link_out = temp_link;
}


void algo_topIP1(ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], ap_uint<576> link_out[10])
{

    #pragma HLS PIPELINE II=9
    #pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
    #pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
    #pragma HLS INTERFACE ap_ctrl_hs port=return

	ap_uint<LINK_WIDTH> __link_in[N_INPUT_LINKS];
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
	for(int link = 0 ; link < N_INPUT_LINKS ; link++){
		__link_in[link] = link_in[link];
	}


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

    ap_fixed<32,16> Ex;
	ap_fixed<32,16> Ey;
	ap_uint<12> HT;

    makeCaloObjects(__link_in, Jets,Taus, Ex, Ey, HT);

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

    l1ct::HadCaloObj egClusters[8]; // TODO . rmove hardcoding
    l1ct::PuppiObj pfSelectedNutrals[N_SECTORS_PF][NNEUTRALS];
    ap_fixed<32,16> Ex_pu;
	ap_fixed<32,16> Ey_pu;
    
    doPFClustringChain(link_in, egClusters , pfSelectedNutrals, Ex_pu, Ey_pu);


    #ifndef __SYNTHESIS__
           if(DEBUG_LEVEL > 0)
           {
                std::cout<<"@@EGClusters | " ;
                for(int i=0;i < 8 ; i++)
                {
                    std::cout<<egClusters[i].hwPt<<","<<egClusters[i].hwEta<<","<<egClusters[i].hwPhi<<" | ";
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


    for(int i=0; i < N_SECTORS_PF; i++ ){
    	packPuppi(pfSelectedNutrals[i], link_out[i]);
	}

    packEG(egClusters, link_out[6]);
    packJets(Jets, link_out[7]);
    packJets(Taus, link_out[8]);

	ap_uint<576> temp_link;
    temp_link.range(11,0) = ap_uint<12>(Ex);
    temp_link.range(23,12) = ap_uint<12>(Ey);
    temp_link.range(35,24) = ap_uint<12>(Ex_pu);
    temp_link.range(47,36) = ap_uint<12>(Ey_pu);
    temp_link.range(59,48) = HT;

    link_out[9] = temp_link;

}

