#include "algo_topIP1.h"


void processInputLinks(ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], hftower HFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI]){
	#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
  #pragma HLS ARRAY_PARTITION variable=HFTowers complete dim=0

for(loop j=0; j<TOWERS_IN_PHI; j=j+2){
    #pragma HLS UNROLL
    for(loop i=0; i<TOWERS_IN_ETA-1; i++){
    #pragma HLS UNROLL
      ap_uint<10> start   = i*10;
      ap_uint<10> end = start + 9;
      HFTowers[i][j] = hftower(link_in[j/2].range(end, start));
  }
    for(loop i=0; i<TOWERS_IN_ETA-1; i++){
    #pragma HLS UNROLL
      ap_uint<10> start   = i*10+110;
      ap_uint<10> end = start + 9;
      HFTowers[i][j+1] = hftower(link_in[j/2].range(end, start));
  }
}

    for(loop j=0; j<TOWERS_IN_PHI; j=j+2){
    #pragma HLS UNROLL
    hftower A10 = HFTowers[TOWERS_IN_ETA-2][j] ;
    hftower B10 = HFTowers[TOWERS_IN_ETA-2][j+1] ;


    ap_uint<8> halfA = A10.energy >> 1 ;
    ap_uint<8> halfB = B10.energy >> 1 ;


    A10.energy = halfA ;
    B10.energy = halfB ;


    HFTowers[TOWERS_IN_ETA-2][j] = A10;
    HFTowers[TOWERS_IN_ETA-2][j+1] = A10;

    HFTowers[TOWERS_IN_ETA-1][j] = B10;
    HFTowers[TOWERS_IN_ETA-1][j+1] = B10;
  }
}

void clustrize( hftower HFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI] , PFcluster  Clusters[N_PF_CLUSTERS]  )
{
	#ifndef __SYNTHESIS__
		for(loop phiI=0;phiI<TOWERS_IN_PHI; phiI++)
		{
			std::cout<<std::setw(5)<<phiI;	
		}
		std::cout<<"  <--Phi"<<"\n";
	 	for(loop etaI=0;etaI<TOWERS_IN_ETA ; etaI++)
	 	{
			for(loop phiI=0;phiI<TOWERS_IN_PHI; phiI++)
			{
				std::cout<<std::setw(5)<< HFTowers[etaI][phiI].energy;	
			}
			std::cout<<"  <--"<<etaI<<"\n";
		}
	#endif
	for(loop clusIdx=0 ; clusIdx < 9 ; clusIdx++)
	{
		uint16_t etaC=0,phiC=0,etmax=0;
	 	for(loop etaI=0;etaI<TOWERS_IN_ETA ; etaI++)
	 	{
			for(loop phiI=0;phiI<TOWERS_IN_PHI ; phiI++)
			{	
				if( (etmax < HFTowers[etaI][phiI].energy ) and ( HFTowers[etaI][phiI].energy > MIN_CLUSTER_SEED_ENERGY) )
				{
					etmax=HFTowers[etaI][phiI].energy;
					etaC=etaI; phiC=phiI;
				}
			}
		}
		if( etmax > 0)
		{
			Clusters[clusIdx].Eta= etaC;
			Clusters[clusIdx].Phi= phiC;
			uint16_t en=etmax;
			uint16_t etaL=etaC-1;
			uint16_t etaH=etaC+1;
			//std::cout<<" Cluster centers : "<<etaC<<" , "<<phiC<<"\n";
			if( etaC == 0  ) etaL = 200;
			if( etaC == TOWERS_IN_ETA-1  ) etaH = 200;
			uint16_t phiL=phiC-1;
			uint16_t phiH=phiC+1;
			if( phiC == 0  ) phiL = TOWERS_IN_PHI-1;
			if( phiC == TOWERS_IN_PHI-1  ) phiH =0 ;
			
			if(etaC>0)
			{
				en = en + HFTowers[etaL][phiL].energy ; 
				en = en + HFTowers[etaL][phiC].energy ;
				en = en + HFTowers[etaL][phiH].energy ;
			}
			if(etaC<TOWERS_IN_ETA-1)
			{
				en = en + HFTowers[etaH][phiL].energy ;
				en = en + HFTowers[etaH][phiC].energy ;
				en = en + HFTowers[etaH][phiH].energy ;
			}

			en = en + HFTowers[etaC][phiL].energy ;
			en = en + HFTowers[etaC][phiH].energy ;
			
			#ifndef __SYNTHESIS__	
			std::cout<<"Cluster Obtained with E = "<<en<<", center = ("<<etaC<<","<<phiC<<","<<etmax<<")\n";
			std::cout<<"     CELLS > "
				 <<" ("<<etaL<<","<<phiL<<","<<HFTowers[etaL][phiL].energy<<"),("<<etaL<<","<<phiC<<","<<HFTowers[etaL][phiC].energy<<"),("<<etaL<<","<<phiH<<","<<HFTowers[etaL][phiL].energy<<")"
                                 <<",("<<etaC<<","<<phiL<<","<<HFTowers[etaC][phiL].energy<<"),("<<etaC<<","<<phiH<<","<<HFTowers[etaC][phiC].energy<<")"                  
				 <<",("<<etaH<<","<<phiL<<","<<HFTowers[etaH][phiL].energy<<"),("<<etaH<<","<<phiC<<","<<HFTowers[etaH][phiC].energy<<"),("<<etaH<<","<<phiH<<","<<HFTowers[etaH][phiH].energy<<")"
				 <<"\n";
			#endif
			HFTowers[etaC][phiC].energy = 0 ;  
			HFTowers[etaC][phiL].energy = 0 ;
                	HFTowers[etaC][phiH].energy = 0 ;
			if(etaC>0)
			{
			HFTowers[etaL][phiL].energy = 0 ;  
                	HFTowers[etaL][phiC].energy = 0 ;
                	HFTowers[etaL][phiH].energy = 0 ;
			}
			if(etaC<TOWERS_IN_ETA-1)
			{
                	HFTowers[etaH][phiL].energy = 0 ;
                	HFTowers[etaH][phiC].energy = 0 ;
                	HFTowers[etaH][phiH].energy = 0 ;
                	}
		}
	}
}

void algo_topIP1(ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], ap_uint<LINK_WIDTH> link_out[N_OUTPUT_LINKS]){
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS PIPELINE II=9
#pragma HLS INTERFACE ap_ctrl_hs port=return

        hftower HFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI] ;
        #pragma HLS ARRAY_PARTITION variable=HFTowers complete dim=0

	processInputLinks(link_in, HFTowers) ;
	
	PFcluster Clusters[N_PF_CLUSTERS];
	clustrize(HFTowers,Clusters);

        ap_uint<10> start ;
        ap_uint<10> end ;

        for(loop j=0; j<TOWERS_IN_PHI; j++)
          {
          	for(loop i=0; i<TOWERS_IN_ETA; i++)
		{
		start=i*10 ; end=start+9;
		link_out[j].range(end, start) = HFTowers[i][j].gettower() ;
	  }}
}

