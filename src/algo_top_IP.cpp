#ifndef ALGO_TOPIP1_CPP
#define ALGO_TOPIP1_CPP

//#define DEBUGME

#ifdef DEBUGME
#include <bitset>
#endif

#include "linpuppi.h"
#include "algo_topIP1.h"
#include <cstdint>

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
#pragma HLS INLINE	

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
        nClusterInRegion[i]=0;
    }

sector_loop: 
    for( loop i=0 ;i < N_SECTORS ; i++) {
    #pragma hls unroll
        #ifdef DEBUGME
        std::cout<<" > Setting region : "<<i
                <<" [ "<<region[i].hwPhiCenter<<"+/-"<< region[i].hwPhiHalfWidth <<"] "<<"\n";
		#endif
    link_loop_per_sector:
		for( loop j=0 ;j < NLINKS_TO_SCAN ; j++) {
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
    const l1ct::PFRegion region[N_SECTORS]
    ) {

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

void algo_topIP1(
    ap_uint<576> link_in[N_INPUT_LINKS],
    ap_uint<576> link_out[N_OUTPUT_LINKS]
    ) {
#pragma HLS PIPELINE

    l1ct::PuppiObj pfselne[N_SECTORS][NNEUTRALS];
    l1ct::HadCaloObj H_in_regionized[N_SECTORS][NCALO];

#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0
// #pragma HLS INTERFACE ap_ctrl_hs port=return
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

// streaming interface -->

void stream_in(
    hls::stream<axi_stream> &inputStream,
    ap_uint<576> *link_in
    ) { 
    
    axi_stream buffer;
    ap_uint<576> tempLink;
    ap_uint<BIT_WIDTH> tempArray[N_WORDS];
    ap_uint<12> start = 0;

    for(loop i=0; i<N_WORDS; i++) {
        buffer = inputStream.read();
        tempLink.range(start+(BIT_WIDTH-1),start) = buffer.data;

        #ifdef DEBUGME

        std::cout << "in_stream buff :" << buffer.data << endl;
        std::cout << "in_stream buff :" << tempLink.range(start+(BIT_WIDTH-1),start) << endl;
        std::cout << "tempLink inside stream_in : " << std::bitset<64>(tempLink.range(start+(BIT_WIDTH-1),start)) << endl;
        std::cout << "start+B , start :" << start+(BIT_WIDTH-1) << "," << start << endl;
        #endif
        
        start = start + BIT_WIDTH;

        if(i == N_WORDS-1){
            start = 0;
        }
    }

    #ifdef DEBUGME
    std::cout << "tempLink inside stream_in : " << std::bitset<576>(tempLink) << endl;
    #endif

    *link_in = tempLink;
}

void stream_out(
    hls::stream<axi_stream> &outputStream,
    ap_uint<576> *link_out
    ) {

    axi_stream buffer;
    ap_uint<576> tempLink;
    ap_uint<12> start = 0;

    tempLink = *link_out;

    #ifdef DEBUGME
    std::cout << "*link_out : " << std::bitset<576>(*link_out) << endl;
    std::cout << "tempLink : " << std::bitset<576>(tempLink) << endl;
    #endif

    for(loop i=0; i<N_WORDS; i++) {
        #ifdef DEBUGME
        std::cout << "start+B , start :" << start+(BIT_WIDTH-1) << "," << start << endl;
        std::cout << "tempLink ranged: " << std::bitset<64>(tempLink.range(start+(BIT_WIDTH-1),start)) << endl;
        #endif
        
        buffer.data = tempLink.range(start+(BIT_WIDTH-1),start);
        start = start + BIT_WIDTH;

        #ifdef DEBUGME
        std::cout << "out_stream :" << buffer.data << endl;
        std::cout << "i insdie the stream_out for loop :" << i << endl;
        #endif

        if (i == N_WORDS-1){
            #ifdef DEBUGME
            cout << "i insdie the stream_out for loop :" << i << endl;            
            #endif
            
            buffer.last = 1;
            start = 0;
        }
        else {
            buffer.last = 0;
        }

        outputStream.write(buffer);
    } 
}

void AXIStream_wrapper(
    hls::stream<axi_stream> &inputStream0,
    hls::stream<axi_stream> &inputStream1,
    hls::stream<axi_stream> &inputStream2,
    hls::stream<axi_stream> &inputStream3,
    hls::stream<axi_stream> &inputStream4,
    hls::stream<axi_stream> &inputStream5,
    // if it works, it works.
    hls::stream<axi_stream> &outputStream0,
    hls::stream<axi_stream> &outputStream1,
    hls::stream<axi_stream> &outputStream2,
    hls::stream<axi_stream> &outputStream3,
    hls::stream<axi_stream> &outputStream4,
    hls::stream<axi_stream> &outputStream5
    ) {

#pragma HLS INTERFACE mode=ap_ctrl_hs port=return
#pragma HLS DATAFLOW
#pragma HLS INTERFACE mode=axis register_mode=both port=inputStream0 register
#pragma HLS INTERFACE mode=axis register_mode=both port=inputStream1 register
#pragma HLS INTERFACE mode=axis register_mode=both port=inputStream2 register
#pragma HLS INTERFACE mode=axis register_mode=both port=inputStream3 register
#pragma HLS INTERFACE mode=axis register_mode=both port=inputStream4 register
#pragma HLS INTERFACE mode=axis register_mode=both port=inputStream5 register
#pragma HLS INTERFACE mode=axis register_mode=both port=outputStream0 register
#pragma HLS INTERFACE mode=axis register_mode=both port=outputStream1 register
#pragma HLS INTERFACE mode=axis register_mode=both port=outputStream2 register
#pragma HLS INTERFACE mode=axis register_mode=both port=outputStream3 register
#pragma HLS INTERFACE mode=axis register_mode=both port=outputStream4 register
#pragma HLS INTERFACE mode=axis register_mode=both port=outputStream5 register
    
    ap_uint<576> link_in[N_INPUT_LINKS];
    ap_uint<576> link_out[N_OUTPUT_LINKS];

    stream_in(inputStream0, &link_in[0]);
    stream_in(inputStream1, &link_in[1]);
    stream_in(inputStream2, &link_in[2]);
    stream_in(inputStream3, &link_in[3]);
    stream_in(inputStream4, &link_in[4]);
    stream_in(inputStream5, &link_in[5]);

     algo_topIP1(link_in, link_out);

//   for(loop i=0; i<6; i++) {
//   	link_out[i] = link_in[i];
//   }

    stream_out(outputStream0, &link_out[0]);
    stream_out(outputStream1, &link_out[1]);
    stream_out(outputStream2, &link_out[2]);
    stream_out(outputStream3, &link_out[3]);
    stream_out(outputStream4, &link_out[4]);
    stream_out(outputStream5, &link_out[5]);

}

#endif
