#include "caloObjects.h"

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
                    superTowers[i+1][j+1].eta = i+1;
                    superTowers[i+1][j+1].phi = j+1;
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

void findMaxEnergySuperTowerInPhi(const hftower PhiTowers[24],
							hftower& phiC)
{
    hftower tempArray[32];
    #pragma HLS ARRAY_PARTITION variable=tempArray complete dim=0

    for(loop i=0; i<24; i++){
        tempArray[i] = PhiTowers[i];
    }

    for(loop i=32; i>1; i=(i/2))
    {
        for(loop j=0; j < i/2; j++)
        {
        	tempArray[j] = bestOf2(tempArray[j*2], tempArray[(j*2) + 1]);
        }
    }

    phiC = tempArray[0];
}


void findMaxEnergySuperTowerInEta(const hftower EtaTowers[4], hftower& etaC)
{
    hftower tempArray[4];
    #pragma HLS ARRAY_PARTITION variable=tempArray complete dim=0

    for(loop i=0; i<4; i++) {
        tempArray[i] = EtaTowers[i];
    }


    for(loop i=4; i>1; i=(i/2))
    {
        for(loop j=0; j < i/2; j++)
        {
        	tempArray[j] = bestOf2(tempArray[j*2], tempArray[(j*2) + 1]);
        }
    }

    etaC = tempArray[0];

}

void findMaxEnergySuperTower(const hftower HFRegion[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2],
                        ap_uint<5>& etaC,
                        ap_uint<8>& phiC)
{
//#pragma HLS INLINE
    #pragma HLS ARRAY_PARTITION variable=HFRegion complete dim=0

	hftower towersPhi[24];
    #pragma HLS ARRAY_PARTITION variable=towersPhi complete dim=0
	hftower tempPhi;
    for(ap_uint<5> phi = 1; phi < 25; phi++)
    {
    	hftower towersEta[4];
        #pragma HLS ARRAY_PARTITION variable=towersEta complete dim=0
        for(loop eta = 1; eta < 5; eta++)
        {
            towersEta[eta-1] = HFRegion[eta][phi];
        }
        findMaxEnergySuperTowerInEta(towersEta, towersPhi[phi-1]);
    }

    findMaxEnergySuperTowerInPhi(towersPhi, tempPhi);

    etaC = tempPhi.eta;
    phiC = tempPhi.phi;
#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL >7 )
    {
        std::cout<<"  Center found at : "<<etaC<<" , "<<phiC<<"\n" ;
    }
#endif
}

void formJetsAndZeroOut(hftower superTowers[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2],
                        ap_uint<10> etaC,
                        ap_uint<10> phiC,
                        ap_uint<18>& etaSum) {

    #pragma HLS ARRAY_PARTITION variable=superTowers complete dim=0

    ap_uint<1> mask[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2] = {0};
    ap_uint<1> zero[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2] = {0};

    for (ap_uint<10> i = 0; i < (TOWERS_IN_ETA / 3) + 2; i++) {
        #pragma HLS UNROLL
        if (i + 1 >= etaC && i <= etaC + 1) {
            for (ap_uint<10> j = 0; j < (TOWERS_IN_PHI / 3) + 2; j++) {
                #pragma HLS UNROLL
                if (j + 1 >= phiC && j <= phiC + 1) {
                    mask[i][j] = 1;
                }
            }
        }
    }

    for (ap_uint<10> i = 0; i < (TOWERS_IN_ETA / 3) + 2; i++) {
        #pragma HLS UNROLL
        for (ap_uint<10> j = 0; j < (TOWERS_IN_PHI / 3) + 2; j++) {
            #pragma HLS UNROLL
            etaSum += superTowers[i][j].energy * mask[i][j];
        }
    }

    for (ap_uint<10> i = 0; i < (TOWERS_IN_ETA / 3) + 2; i++) {
        #pragma HLS UNROLL
        for (ap_uint<10> j = 0; j < (TOWERS_IN_PHI / 3) + 2; j++) {
            #pragma HLS UNROLL
            zero[i][j] = (ap_uint<1>)(1 - mask[i][j]);
        }
    }

    for (ap_uint<10> i = 0; i < (TOWERS_IN_ETA / 3) + 2; i++) {
        #pragma HLS UNROLL
        for (ap_uint<10> j = 0; j < (TOWERS_IN_PHI / 3) + 2; j++) {
            #pragma HLS UNROLL
            superTowers[i][j].energy *= zero[i][j];
        }
    }
}

void selectTaus(const jets Jet[9], jets Taus[9])
{

#pragma HLS ARRAY_PARTITION variable=Jet complete dim=0
#pragma HLS ARRAY_PARTITION variable=Taus complete dim=0

    ap_uint<4> count = 0;
    for(loop idx = 0; idx < 9; idx++)
    {
        // todo : fix to multiplication
        if((Jet[idx].seedET * 10) >= (Jet[idx].ET * 7))
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

#pragma HLS ARRAY_PARTITION variable=Jet complete dim=0


	hftower __superTowers[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2];
#pragma HLS ARRAY_PARTITION variable=__superTowers complete dim=0

    for (ap_uint<10> i = 0; i < (TOWERS_IN_ETA/3)+2; i++) {
        for (ap_uint<10> j = 0; j < (TOWERS_IN_PHI/3)+2; j++) {
        	__superTowers[i][j] = superTowers[i][j];
        }
    }

    for(loop idx = 0; idx < 9; idx++)
    {
        ap_uint<8> phiC = 1;
        ap_uint<5> etaC = 1;
        ap_uint<18> etSum = 0;

        findMaxEnergySuperTower(__superTowers, etaC, phiC);

        Jet[idx].seedET = __superTowers[etaC][phiC].energy;
        formJetsAndZeroOut(__superTowers, etaC, phiC, etSum);

//        if(etSum ==0 )
//        {
//            etaC=0;
//            phiC=0;
//        }
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

void Exy(const ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS],
		ap_fixed<32,16>& Ex,
		ap_fixed<32,16>& Ey) {

	ap_fixed<32,16> sin_LUT[36] = {0.0872, 0.2588, 0.4226, 0.5736, 0.7071, 0.8192, 0.9063, 0.9659, 0.9962, 0.9962, 0.9659, 0.9063,
	                               0.8192, 0.7071, 0.5736, 0.4226, 0.2588, 0.0872, -0.0872, -0.2588, -0.4226, -0.5736, -0.7071, -0.8192,
	                               -0.9063, -0.9659, -0.9962, -0.9962, -0.9659, -0.9063, -0.8192, -0.7071, -0.5736, -0.4226, -0.2588, -0.0872};

    ap_fixed<32,16> Ex_temp = 0;
    ap_fixed<32,16> Ey_temp = 0;

    for(loop j = 0; j < TOWERS_IN_PHI/2; j=j+2) {
        #pragma HLS UNROLL
        for(loop i = 0; i < TOWERS_IN_ETA - 2; i++) {
            #pragma HLS UNROLL

            ap_uint<9> A_Energy = link_in[j/2].range(i*10 + 7, i*10);
            ap_uint<9> B_Energy = link_in[j/2].range(i*10 + 117, i*10 + 110);

            A_Energy = A_Energy << 1;
            B_Energy = B_Energy << 1;

            Ey_temp += A_Energy * sin_LUT[j];
            Ex_temp += A_Energy * sin_LUT[(j + 9) % 36];

            Ey_temp += B_Energy * sin_LUT[j+1];
            Ex_temp += B_Energy * sin_LUT[(j+10) % 36];

        }

        ap_uint<8> A10 = link_in[j/2].range(107, 100);
        ap_uint<8> B10 = link_in[j/2].range(217, 210);

        Ey_temp += (A10 * sin_LUT[j])*2;
        Ex_temp += (A10 * sin_LUT[(j+9) % 36])*2;

        Ey_temp += (B10 * sin_LUT[j+1])*2;
        Ey_temp += (B10 * sin_LUT[(j+10) % 36])*2;
    }


    Ex = Ex_temp >> 1;
    Ey = Ey_temp >> 1;


    #ifndef __SYNTHESIS__
    std::cout << "\nEx: " << Ex << "\nEy: " << Ey << std::endl;
    #endif
}

void makeCaloObjects(const ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS],
                     jets Jets[9], jets Taus[9], ap_fixed<32,16>& Ex, ap_fixed<32,16>& Ey, ap_uint<12>& HT)
{

//    #pragma HLS PIPELINE II=9

	ap_uint<LINK_WIDTH> __link_in[N_INPUT_LINKS];
#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
	for(int link = 0 ; link < N_INPUT_LINKS ; link++){
		__link_in[link] = link_in[link];
	}

	hftower superTowerGrid[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2];
#pragma HLS ARRAY_PARTITION variable=superTowerGrid complete dim=0
    jets _Jets[9];
#pragma HLS ARRAY_PARTITION variable=_Jets complete dim=0
    jets _Taus[9];
#pragma HLS ARRAY_PARTITION variable=superTowerGrid complete dim=0
    ap_uint<12> temp_HT;

    unpackToSuperTowers(link_in, superTowerGrid);

    Exy(__link_in, Ex, Ey);

    makeJets(superTowerGrid, _Jets);

    selectTaus(_Jets, _Taus);

    for(loop idx = 0; idx < 9; idx++)
    {
        Jets[idx] = _Jets[idx];
        temp_HT += _Jets[idx].ET;
    }

    for(loop idx = 0; idx < 9; idx++)
    {
        Taus[idx] = _Taus[idx];
    }

    HT = temp_HT;
}





