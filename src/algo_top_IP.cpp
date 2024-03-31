#include "linpuppi.h"
#include "algo_topIP1.h"
#include <cstdint>

//#define DEBUGME

void packLinks( 
    l1ct::PuppiObj puppiobjs[N_SECTORS][NNEUTRALS] , 
    ap_uint<576> link_out[N_OUTPUT_LINKS] 
    ) { // packs the puppi objects to the output links
#pragma HLS INLINE

    const ap_uint<8> BW = 64;
    ap_uint<576> temp_link[N_OUTPUT_LINKS] = {0};
    #pragma HLS ARRAY_PARTITION dim=1 type=complete variable=temp_link
pack:
    for(loop i=0; i<N_SECTORS; i++) {
    #pragma HLS unroll
        for(loop k=0; k<NNEUTRALS; k++) {
        #pragma HLS unroll
            ap_uint<6> i_link = (i*NNEUTRALS+k)/N_PUPPI_LINK;
            ap_uint<5> offset = (i*NNEUTRALS+k) % N_PUPPI_LINK;
            ap_uint<12> start = offset*BW;

            temp_link[i_link].range(start+BW,start) = puppiobjs[i][k].pack();
        }
    }
    
memcopy: 
    for(loop i=0; i<N_OUTPUT_LINKS; i++){
    #pragma HLS UNROLL

    	link_out[i] = temp_link[i];
	}
}

void unpackLinks( 
    const ap_uint<576> link_in[N_INPUT_LINKS] , 
    ap_uint<12>  ET[N_INPUT_LINKS][N_PF_LINK] ,
    ap_uint<8>  ETA[N_INPUT_LINKS][N_PF_LINK] , 
    ap_int<9>   PHI[N_INPUT_LINKS][N_PF_LINK] 
    ) { // unpacks the puppi objects from the input links
#pragma HLS INLINE
unpack_link_linkLoop:
    for(loop i=0; i<N_INPUT_LINKS; i++) {
    #pragma HLS unroll   
        ap_uint<576> word= link_in[i];
        unpack_link_objLoop : 
            for(loop j=0; j<N_PF_LINK; j++) {
            #pragma HLS PIPELINE II=1
                etW :     ET[i][j]     = word.range(11,0);
                etaW:    ETA[i][j]     = word.range(16,12);
                phiW:    PHI[i][j]     = word.range(24,17);
                bShift:    word=word>>64;
                #ifdef DEBUGME    
                std::cout<<"Unpacked object : i,j == "<<i<<","<<j<<" : "
                        << "et  "<< ET[i][j]<<" | "
                        << "eta "<<ETA[i][j]<<" | "
                        << "phi "<<PHI[i][j]<<" \n";
                #endif
        }
    }

}

void regionizeToHCALObjects(
    const ap_uint<12>  ET[N_INPUT_LINKS][N_PF_LINK] ,
    const ap_uint<8>  ETA[N_INPUT_LINKS][N_PF_LINK] ,
    const ap_int<9>   Phi[N_INPUT_LINKS][N_PF_LINK] ,
    const l1ct::PFRegion region[N_SECTORS],
    l1ct::HadCaloObj H_in_regionized[N_SECTORS][NCALO]
    ) { 
#pragma HLS INLINE off

const uint8_t NLINKS_TO_SCAN(3);

const ap_uint<4> REGION_TO_LINK_MAP[N_SECTORS][NLINKS_TO_SCAN]={
    {0,1,5},
    {1,2,0},
    {2,3,1},
    {3,4,2},
    {4,5,3},
    {5,0,4}
    };
#pragma HLS ARRAY_PARTITION variable=REGION_TO_LINK_MAP dim=0

int nClusterInRegion[N_SECTORS] ; // can be changed to uint4
#pragma HLS ARRAY_PARTITION variable=nClusterInRegion dim=0

    for( loop i=0 ;i < N_SECTORS ; i++) {
        #pragma HLS unroll
        nClusterInRegion[i]=N_PF_LINK;
    }

sector_loop: 
    for( loop i=0 ;i < N_SECTORS ; i++) {
    #pragma hls unroll
        #ifdef DEBUGME
        std::cout<<" > Setting region : "<<i
                <<" [ "<<region[i].hwPhiCenter<<"+/-"<< region[i].hwPhiHalfWidth <<"] "<<"\n";
		#endif
	
	const loop lnk1 = REGION_TO_LINK_MAP[i][0];
    central_link_loop_per_sector:
	for(loop k=0; k < N_PF_LINK ; k++)
	{
		localEta_def_subs :	ap_int<9>  local_eta = ETA[lnk1][k]-  region[i].hwEtaCenter;
		localPhi_def_sum  :	ap_int<9>  local_phi = Phi[lnk1][k]-  region[i].hwPhiCenter;
					ap_int<12> local_et  =  ET[lnk1][k];
		HCALObjWriteEta :       H_in_regionized[i][k].hwEta = local_eta;
		HCALObjWritePhi :       H_in_regionized[i][k].hwPhi = local_phi;
		HCALObjWriteEt  :       H_in_regionized[i][k].hwPt =  ET[lnk1][k];
	
			#ifdef DEBUGME
			std::cout<<"  > Adding the def. i,j : "<<i<<","<<k<<" : "
					<<"("<<Phi[lnk1][k]<<"->"<<local_phi<<", "<<ETA[lnk1][k]<<"->"<<local_eta<<"   et : "<<local_et
					<<" ) is inside region "<<i<<" : "<<region[i].isInside(local_eta, local_phi)
					<<"\n" ;
			#endif
	}
	
	const loop lnkL = REGION_TO_LINK_MAP[i][1];
	const loop lnkR = REGION_TO_LINK_MAP[i][2];
	 loop idxL(0),idxR(0);
	int phi_offsetL=region[i].hwPhiCenter;
	int phi_offsetR=region[i].hwPhiCenter;
	int hid;
	ap_int<9> local_etaR,local_phiR;		
	ap_int<9> local_etaL,local_phiL;		
        ap_uint<12> local_etR,local_etL;
	if     (  i==0   and    lnkL==(N_INPUT_LINKS-1) ) phi_offsetL+=72;
	else if(  i==(N_SECTORS-1) and    lnkR==0       ) phi_offsetR-=72;

     link_loop_per_sector:
      for(int kk=0;kk< N_PF_LINK*2; kk++)
	{
	//#pragma HLS loop_tripcount max=16
		// TODO SET THE MAX TRIP COUNT
		if ( (idxR >= N_PF_LINK) and ( idxL >= N_PF_LINK)) break;
		if(nClusterInRegion[i] >= NCALO ) break;
		bool isLCand=false;
		if (idxL < N_PF_LINK)
		{
			local_phiL = Phi[lnkL][idxL]-  phi_offsetL;
			local_etaL = ETA[lnkL][idxL]-region[i].hwEtaCenter;
			local_etL  = ET[lnkL][idxL] ;
			
			#ifdef DEBUGME
			std::cout<<"  > Checking L(r"<<i<<",l"<<lnkL<<",p"<<lnkL<<") : "
					<<"("<<Phi[lnkL][idxL]<<"->"<<local_phiL<<", "<<ETA[lnkL][idxL]<<"->"<<local_etaL<<"   et : "<<local_etL
					<<" ) is inside region "<<i<<" : "<<region[i].isInside(local_etaL, local_phiL)
					<<"\n" ;
			#endif

			if( region[i].isInside(local_etaL, local_phiL ) )
				isLCand=true;	
		}

		bool isRCand=false;
		if (isRCand < N_PF_LINK)
		{
			local_phiR = Phi[lnkR][idxR]-  phi_offsetR;
			local_etaR = ETA[lnkR][idxR]-region[i].hwEtaCenter;
			local_etR  = ET[lnkR][idxR] ;
			#ifdef DEBUGME
			std::cout<<"  > Checking R(r"<<i<<",l"<<lnkR<<",p"<<lnkR<<") : "
					<<"("<<Phi[lnkR][idxR]<<"->"<<local_phiR<<", "<<ETA[lnkR][idxR]<<"->"<<local_etaR<<"   et : "<<local_etR
					<<" ) is inside region "<<i<<" : "<<region[i].isInside(local_etaR, local_phiR)
					<<"\n" ;
			#endif
			if( region[i].isInside(local_etaR, local_phiR ) )
				isRCand=true;	
		}

		hid=nClusterInRegion[i];
		if (isLCand or isRCand)
			nClusterInRegion[i]++;
		if (isLCand and isRCand)
		{
			if ( ET[lnkL][idxL] < ET[lnkR][idxR]){
			  H_in_regionized[i][hid].hwEta = local_etaR;
			  H_in_regionized[i][hid].hwPhi = local_phiR;
			  H_in_regionized[i][hid].hwPt =  ET[lnkR][idxR];
			#ifdef DEBUGME
			std::cout<<"  >    adding L as "<<hid<<" \n";
			#endif
	                  idxR++;
			}
			else{
			#ifdef DEBUGME
			std::cout<<"  >    adding R as "<<hid<<" \n";
			#endif
			  H_in_regionized[i][hid].hwEta = local_etaL;
			  H_in_regionized[i][hid].hwPhi = local_phiL;
			  H_in_regionized[i][hid].hwPt =  ET[lnkL][idxL];
	                  idxL++;
		     }
		}
		else if (isLCand)
		{
			  H_in_regionized[i][hid].hwEta = local_etaL;
			  H_in_regionized[i][hid].hwPhi = local_phiL;
			  H_in_regionized[i][hid].hwPt =  ET[lnkL][idxL];
		          std::cout<<"  >    adding L as "<<hid<<" \n";
	                  idxL++;
	                  idxR++;
		}
		else if (isRCand)
		{
			  H_in_regionized[i][hid].hwEta = local_etaR;
			  H_in_regionized[i][hid].hwPhi = local_phiR;
		 	  H_in_regionized[i][hid].hwPt =  ET[lnkR][idxR];
		          std::cout<<"  >    adding R as "<<hid<<" \n";
	                  idxL++;
	                  idxR++;
		}
		else 
		{
			#ifdef DEBUGME
			std::cout<<"  >    updating both idxs from L"<<idxL<<",R"<<idxR<<"\n";
			#endif
	                  idxL++;
	                  idxR++;
		}
	}
    }

   #ifdef DEBUGME    
   for( loop i=0 ;i < N_SECTORS ; i++)
   {
       for( loop j=0 ;j < NCALO ; j++)
        {
            std::cout<<"HCALOobj "<<i<<" , "<<j
             <<" : phi "<<H_in_regionized[i][j].hwPhi
             <<" : eta "<<H_in_regionized[i][j].hwEta
             <<" : et  "<<H_in_regionized[i][j].hwPt
             <<"\n";

        }
   }
   #endif
	
}

void compute(
		l1ct::PuppiObj pfselne[N_SECTORS][NNEUTRALS],
		l1ct::HadCaloObj H_in_regionized[N_SECTORS][NCALO],
		const l1ct::PFRegion region[N_SECTORS]
		){
# pragma HLS INLINE off
	puppi:
	    for(loop i=0; i<N_SECTORS; i++) {
	         fwdlinpuppi(region[i], H_in_regionized[i], pfselne[i]);
	    }
   #ifdef DEBUGME    
   for( loop i=0 ;i < N_SECTORS ; i++)
   {
       for( loop j=0 ;j < NCALO ; j++)
        {
            std::cout<<"Puppi Obj : "<<i<<" , "<<j
             <<" : phi "<<pfselne[i][j].hwPhi
             <<" : eta "<<pfselne[i][j].hwEta
             <<" : et  "<<pfselne[i][j].hwPt
             <<"\n";

        }
   }
   #endif
	



}

// Unpacks the
void unpackLinksToHadCalo(
    const ap_uint<576> link_in[N_INPUT_LINKS] ,
    l1ct::HadCaloObj H_in_regionized[N_SECTORS][NCALO], 
    const l1ct::PFRegion region[N_SECTORS]
    ) {

#pragma HLS INLINE off
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=H_in_regionized complete dim=0

    ap_uint<12>  ET[N_INPUT_LINKS][N_PF_LINK] ;
    ap_uint<8>  ETA[N_INPUT_LINKS][N_PF_LINK] ;
    ap_int<9>   PHI[N_INPUT_LINKS][N_PF_LINK] ;

#pragma HLS ARRAY_PARTITION variable=ET dim=0
#pragma HLS ARRAY_PARTITION variable=ETA dim=0
#pragma HLS ARRAY_PARTITION variable=PHI dim=0
    fCall_unpackLink : unpackLinks(link_in , ET, ETA, PHI) ;
    fCall_regionizer : regionizeToHCALObjects( ET, ETA , PHI, region,H_in_regionized) ;

}

void algo_topIP1(
    ap_uint<576> link_in[N_INPUT_LINKS],
    ap_uint<576> link_out[N_OUTPUT_LINKS]
    ) {
#pragma HLS PIPELINE II=8

    l1ct::PuppiObj pfselne[N_SECTORS][NNEUTRALS];
    l1ct::HadCaloObj H_in_regionized[N_SECTORS][NCALO];

#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
#pragma HLS INTERFACE ap_ctrl_hs port=return
#pragma HLS ARRAY_PARTITION variable=H_in_regionized complete dim=0
#pragma HLS ARRAY_PARTITION variable=pfselne complete dim=0

    ap_uint<12> start, end;

    // define the  sector boundaries and overlaps
    l1ct::PFRegion region[N_SECTORS];
regions_init:
    for(int i=0 ; i < N_SECTORS ; i++) {
        region[i].hwEtaCenter = l1ct::glbeta_t(12);
        region[i].hwEtaHalfWidth = l1ct::eta_t(12);
        region[i].hwEtaExtra = l1ct::eta_t(0);
        region[i].hwPhiExtra = l1ct::phi_t(2);
        region[i].hwPhiHalfWidth = l1ct::phi_t(6);
        region[i].hwPhiCenter = l1ct::glbphi_t(12*i+5);
    }



    // initialize the neutral hadron objects
hadcalo_init:
    for( loop  i= 0; i < N_SECTORS; i++) {
    #pragma HLS unroll
        for( loop j = 0; j < NCALO ; j++) {
        #pragma HLS unroll
            H_in_regionized[i][j].clear();
        }
    }

    // unpack and regionize the PF Clusters to the neutral hadron objects
    
    unpackLinksToHadCalo(link_in, H_in_regionized, region);

    // invoke puppi for the regionized clusters
    for( loop i=0 ;i < N_SECTORS ; i++)
    {
        for( loop j=0 ;j < NCALO ; j++)
         {
             std::cout<<"HCALOobj "<<i<<" , "<<j
              <<" : phi "<<H_in_regionized[i][j].hwPhi
              <<" : eta "<<H_in_regionized[i][j].hwEta
              <<" : et  "<<H_in_regionized[i][j].hwPt
              <<"\n";

         }
    }


    compute(pfselne, H_in_regionized, region);


    // pack the puppi outputs into links
    packLinks(pfselne, link_out);
}
