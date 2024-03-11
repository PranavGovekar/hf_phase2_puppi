#include "linpuppi.h"
#include "algo_topIP1.h"
#include <cstdint>


#include <bitset>

// #define DEBUGME 1
// packs the puppi objects to the output links
void packLinks(l1ct::PuppiObj puppiobjs[N_SECTORS][NNEUTRALS], ap_uint<576> link_out[N_OUTPUT_LINKS] ) {
    #pragma HLS INLINE
    const ap_uint<8> BW = 64;
    ap_uint<576> temp_link[N_OUTPUT_LINKS] = {0};
#pragma HLS ARRAY_PARTITION dim=1 type=complete variable=temp_link

pack:
    for(loop i=0; i<N_SECTORS; i++)
    {
     #pragma HLS unroll
        for(loop k=0; k<NNEUTRALS; k++)
        {
         #pragma HLS unroll
            ap_uint<6> i_link = (i*NNEUTRALS+k)/N_PUPPI_LINK;
            ap_uint<5> offset = (i*NNEUTRALS+k) % N_PUPPI_LINK;
            ap_uint<12> start = offset*BW;

            temp_link[i_link].range(start+BW,start) = puppiobjs[i][k].pack();
        }
    }
//    std::copy(std::begin(temp_link), std::end(temp_link), std::begin(link_out));
    memcopy: for(loop i=0; i<N_OUTPUT_LINKS; i++){
#pragma HLS UNROLL
    	link_out[i] = temp_link[i];
	}
}

void unpackLinks( const ap_uint<576> link_in[N_INPUT_LINKS] ,
                  ap_uint<12>  ET[N_INPUT_LINKS][N_PF_LINK] ,
                  ap_uint<8>  ETA[N_INPUT_LINKS][N_PF_LINK] ,
		          ap_int<9>   PHI[N_INPUT_LINKS][N_PF_LINK]
		)
{
    #pragma HLS INLINE
unpack_link_linkLoop:
    for(loop i=0; i<N_INPUT_LINKS; i++){
    #pragma HLS unroll
    ap_uint<576> word= link_in[i];
        unpack_link_objLoop : for(loop j=0; j<N_PF_LINK; j++)
        {
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
//     }

}

void regionizeToHCALObjects(
                  const ap_uint<12>  ET[N_INPUT_LINKS][N_PF_LINK] ,
                  const ap_uint<8>  ETA[N_INPUT_LINKS][N_PF_LINK] ,
		  const ap_int<9>   Phi[N_INPUT_LINKS][N_PF_LINK] ,
		  const l1ct::PFRegion region[N_SECTORS],
	 	  l1ct::HadCaloObj H_in_regionized[N_SECTORS][NCALO]
    )
{
    #pragma HLS INLINE
    // NEED TO DECIDE WHERE TO ACTUALLY PUT THIS , maybe we can move the defenition of "REGIONS" here ?
    const uint8_t NLINKS_TO_SCAN(3);
    const ap_uint<4> REGION_TO_LINK_MAP[N_SECTORS][NLINKS_TO_SCAN]={
					{0,1,5},
					{1,2,0},
					{2,3,1},
					{3,4,2},
					{4,5,3},
					{5,0,4}
			} ;
    #pragma HLS ARRAY_PARTITION variable=REGION_TO_LINK_MAP dim=0
    int nClusterInRegion[N_SECTORS] ; // can be changed to uint4
    #pragma HLS ARRAY_PARTITION variable=nClusterInRegion dim=0
    for( loop i=0 ;i < N_SECTORS ; i++)
    {
        #pragma HLS unroll
        nClusterInRegion[i]=0;
    }

 sector_loop:
   for( loop i=0 ;i < N_SECTORS ; i++)
    {
      #pragma hls unroll
		#ifdef DEBUGME
		  std::cout<<" > Setting region : "<<i
					<<" [ "<<region[i].hwPhiCenter<<"+/-"<< region[i].hwPhiHalfWidth <<"] "<<"\n";
		#endif
	  link_loop_per_sector:
		for( loop j=0 ;j < NLINKS_TO_SCAN ; j++)
		{

			int phi_offset=region[i].hwPhiCenter;
			const loop lnk = REGION_TO_LINK_MAP[i][j];
			// Assumes that the regions are defines such that the there is only division in phi
			// Assumes that the regions have monotonically increasing phi boundaries
			if(  i==0             and    lnk==(N_INPUT_LINKS-1)      ) phi_offset+=72;
			else if(  i==(N_SECTORS-1) and    lnk==0                 ) phi_offset-=72;
		  paticle_loop:
			for(loop k=0; k < N_PF_LINK ; k++)
			{
				#pragma HLS PIPELINE II=2
				nClusterMaxCheck  :	if( nClusterInRegion[i] < NCALO)
				{
					localPhi_def_sum  :	ap_int<9>  local_phi = Phi[lnk][k]-  phi_offset;
					localEta_def_subs :	ap_int<9>  local_eta = ETA[lnk][k]-region[i].hwEtaCenter;
					isInsideEval :      bool isInside = region[i].isInside(local_eta, local_phi )    ;
					ap_uint<4> hid= nClusterInRegion[i];
					isInsideCheck1 :    if(isInside)
					NClusIncrement :     {
												nClusterInRegion[i]++;
										 }
					#ifdef DEBUGME
					std::cout<<"  > Checking (r"<<i<<",l"<<lnk<<",p"<<k<<") : "
							<<"("<<Phi[lnk][k]<<"->"<<local_phi<<", "<<ETA[lnk][k]<<"->"<<local_eta
							<<") is inside region "<<i<<" : "<<region[i].isInside(local_eta, local_phi)
							<<"\n" ;
					#endif
					isInsideCheck2  :   if(isInside){
					HCALObjWriteEta :        	    H_in_regionized[i][hid].hwEta = local_eta;
					HCALObjWritePhi :       	    H_in_regionized[i][hid].hwPhi = local_phi;
					HCALObjWriteEt  :        	    H_in_regionized[i][hid].hwPt =  ET[lnk][k];
					#ifdef DEBUGME
									std::cout<<"    > Adds cluster with et "<<ET[lnk][k]
									<<" as "<<nClusterInRegion[i]<<" item \n";
					#endif
								}
				}
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

// Unpacks the
void unpackLinksToHadCalo(
             const ap_uint<576> link_in[N_INPUT_LINKS] ,
             l1ct::HadCaloObj H_in_regionized[N_SECTORS][NCALO],
             const l1ct::PFRegion region[N_SECTORS])
{
    #pragma HLS INLINE
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

void algo_topIP1(ap_uint<576> link_in[N_INPUT_LINKS],ap_uint<576> link_out[N_OUTPUT_LINKS]) {

    l1ct::PuppiObj pfselne[N_SECTORS][NNEUTRALS];
    l1ct::HadCaloObj H_in_regionized[N_SECTORS][NCALO];

    #pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
    #pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
    #pragma HLS INTERFACE ap_ctrl_hs port=return
    #pragma HLS ARRAY_PARTITION variable=H_in_regionized complete dim=0
    #pragma HLS ARRAY_PARTITION variable=pfselne complete dim=0

//	#pragma HLS PIPELINE

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



    // initialize the neutral hadron objects
//hadcalo_init:
//    for( loop  i= 0; i < N_SECTORS; i++) {
//     #pragma HLS unroll
//        for( loop j = 0; j < NCALO ; j++)
//        {
//         #pragma HLS unroll
//            H_in_regionized[i][j].clear();
//        }
//    }

    // unpack and regionize the PF Clusters to the neutral hadron objects

    unpackLinksToHadCalo(link_in, H_in_regionized, region);


    // invoke puppi for the regionized clusters
puppi:
    for(loop i=0; i<N_SECTORS; i++) {
     #pragma HLS unroll
         fwdlinpuppi(region[i], H_in_regionized[i], pfselne[i]);
//  proxy for puppi
//      proxyLoop :  for(loop j=0; j<NNEUTRALS ;j++)
//        {
//	   #pragma HLS unroll
//            pfselne[i][j].hwPt  = H_in_regionized[i][j].hwPt  ;
//            pfselne[i][j].hwEta = region[i].hwGlbEta( H_in_regionized[i][j].hwEta );
//            pfselne[i][j].hwPhi = region[i].hwGlbPhi(H_in_regionized[i][j].hwPhi );
//        }


    }

    // pack the puppi outputs into links
    packLinks(pfselne, link_out);
}

void ReadWrite( ap_uint<64> in[54], 
                ap_uint<64> out[54]
            ){

    #pragma HLS INTERFACE m_axi port=in bundle=aximm1
    #pragma HLS INTERFACE m_axi port=out bundle=aximm2

    const ap_uint<8> BW = 64;
    const ap_uint<8> N_WORDS = 9;

    ap_uint<576> link_in[N_INPUT_LINKS];
    ap_uint<576> link_out[N_OUTPUT_LINKS];
    ap_uint<576> tempLink;
    ap_uint<12> start = 0;

    for(loop i=0; i<N_INPUT_LINKS; i++) {
        for(loop j=0; j<N_WORDS; j++) {
            #pragma HLS PIPELINE
            #ifdef DEBUGME
            std::cout << "\nelement : " << i*N_WORDS+j;
            std::cout << "\n" << start << endl;
            std::cout << in[(i*N_WORDS)+j] <<endl;
            std::cout << std::bitset<64>(in[(i*N_WORDS)+j]) << endl;
            #endif
            tempLink.range(start+(BW-1),start) = in[(i*N_WORDS)+j];
            start = start + BW;
        }
        #ifdef DEBUGME
        std::cout << "\ntempLink " << i << " : " <<  tempLink;
        std::cout << "\ntempLink binary " << i << " : " <<  std::bitset<576>(tempLink);
        #endif

        link_in[i] = tempLink;
        start = 0;
    }
    
    #ifdef DEBUGME
    std::cout << "HERE+++++++++++++++++++++++++++" << endl;
    for(loop i=0; i<N_OUTPUT_LINKS; i++) {
        std::cout << std::bitset<576>(link_in[i]) << endl;
    }
    std::cout << "HERE+++++++++++++++++++++++++++" << endl;
    #endif

    algo_topIP1(link_in,link_out);

    #ifdef DEBUGME
    std::cout << "HERE+++++++++++++++++++++++++++" << endl;
    for(loop i=0; i<N_OUTPUT_LINKS; i++) {
        std::cout << std::bitset<576>(link_out[i]) << endl;
    }
    std::cout << "HERE+++++++++++++++++++++++++++" << endl;
    #endif

    for(loop i=0; i<N_OUTPUT_LINKS; i++) {
        tempLink = link_out[i];
        for(loop j=0; j<N_WORDS; j++) {
            #pragma HLS PIPELINE
            out[(i*N_WORDS)+j] = tempLink.range(start+(BW-1),start);
            start = start + BW;
        }
        start = 0;
    }
}
