#include "algo_topIP1.h"

void packPuppi(	l1ct::PuppiObj pfselne[NNEUTRALS],
			ap_uint<576> &link_out){
#pragma HLS INLINE
	const ap_uint<8> BW = 64;
	ap_uint<12> start = 0;
	ap_uint<576> temp_link;

	packPuppi_1_1:
    for(loop idx=0; idx<NNEUTRALS; idx++) {
    	temp_link.range(start+(BW-1),start) = pfselne[idx].pack();

        start = start + BW;
    }

    temp_link.range(575, 512) = 0;
    link_out = temp_link;
}

void sendOut(ap_uint<576> temp_link[10], ap_uint<576> link_out[10]){
#pragma HLS INLINE off

	for(loop idx=0; idx<10; idx++) {
		link_out[idx] = temp_link[idx];
	}

}

void packEG(l1ct::HadCaloObj pfselne[NNEUTRALS],
			ap_uint<576> &link_out){
	const ap_uint<8> BW = 64;
	ap_uint<12> start = 0;
	ap_uint<576> temp_link;

	packEG_1_1:
    for(loop idx=0; idx<8; idx++) {
//    	std::cout << "EG pack : " << idx << " | " << pfselne[idx].pack() << std::endl;
    	temp_link.range(start+(BW-1),start) = pfselne[idx].pack();
//    	std::cout << "end : " << start+(BW-1) << " start : " << start << std::endl;
        start = start + BW;
    }
    temp_link.range(575, 512) = 0;
    link_out = temp_link;
}

void packJets(jets Jets[9], ap_uint<576> &link_out){
	const ap_uint<8> BW = 64;
	ap_uint<12> start = 0;
	ap_uint<576> temp_link;

	packJets_1_1:
    for(loop idx=0; idx<9; idx++) {
//    	std::cout << "JETS pack : " << idx << " | " << Jets[idx].data() << std::endl;
    	temp_link.range(start+(BW-1),start) = Jets[idx].data();

        start = start + BW;
    }
    link_out = temp_link;
}


void algo_topIP1(ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], ap_uint<576> link_out[10])
{
    #pragma HLS PIPELINE II=9
    #pragma HLS ARRAY_PARTITION variable=link_in type=complete
    #pragma HLS ARRAY_PARTITION variable=link_out type=complete
    #pragma HLS INTERFACE ap_ctrl_hs port=return

//	for(int link = 0 ; link < 10 ; link++){
//		link_out[link] = 0;
//	}


	ap_uint<LINK_WIDTH> __link_in[N_INPUT_LINKS];
#pragma HLS ARRAY_PARTITION variable=__link_in type=complete
	algo_topIP1_1_1:
	for(int link = 0 ; link < N_INPUT_LINKS ; link++){
		__link_in[link] = link_in[link];
	}


    jets Jets[9];
    jets Taus[9];
    #pragma HLS ARRAY_PARTITION variable=Jets type=complete
    #pragma HLS ARRAY_PARTITION variable=Taus type=complete


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
                    std::cout<<std::setw(3)<<towerGrid[eta][phi].energy<<" | ";
                }
                std::cout << std::endl;
            }
        }
    }
#endif

    ap_fixed<32,16> Ex;
	ap_fixed<32,16> Ey;
	ap_uint<12> HT;

    makeCaloObjects(__link_in, Jets,Taus, Ex, Ey, HT);

//    std::cout<<"JETS : \n";
//    for(int i=0 ; i< 9 ; i++)
//    {
//        std::cout<< "ET : " << Jets[i].ET << "Eta : " << Jets[i].Eta << "Phi : " << Jets[i].Phi << std::endl;
//    }
//    std::cout<<"\n";

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

    l1ct::HadCaloObj egClusters[9]; // TODO . rmove hardcoding
#pragma HLS ARRAY_PARTITION variable=egClusters type=complete

    l1ct::PuppiObj pfSelectedNutrals[N_SECTORS_PF][NNEUTRALS];
#pragma HLS ARRAY_PARTITION variable=pfSelectedNutrals type=complete

    ap_fixed<32,16> Ex_pu;
	ap_fixed<32,16> Ey_pu;
    
    doPFClustringChain(link_in, egClusters , pfSelectedNutrals, Ex_pu, Ey_pu);


    #ifndef __SYNTHESIS__
           if(DEBUG_LEVEL > 0)
           {
                std::cout<<"@@EGClusters | " ;
                for(int i=0;i < 9 ; i++)
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

    ap_uint<576> temp_link[10];
    for(int i=0; i < 10; i++ ){
    	temp_link[i] = 0;
    }

    algo_topIP1_2_1:
    for(int i=0; i < N_SECTORS_PF; i++ ){
    	packPuppi(pfSelectedNutrals[i], temp_link[i]);
	}

    packEG(egClusters, temp_link[6]);
    packJets(Jets, temp_link[7]);
    packJets(Taus, temp_link[8]);

    temp_link[9].range(11,0) = ap_uint<12>(Ex);
    temp_link[9].range(23,12) = ap_uint<12>(Ey);
    temp_link[9].range(35,24) = ap_uint<12>(Ex_pu);
    temp_link[9].range(47,36) = ap_uint<12>(Ey_pu);
    temp_link[9].range(59,48) = HT;

    sendOut(temp_link, link_out);


//    std::cout << std::endl << "INSIDE:" << std::endl;
//    for(loop i = 0; i < 10; i++){
//    	std::cout << std::hex << (link_out[i]) << std::dec << std::endl;
//    }
//    std::cout << std::endl;

}

