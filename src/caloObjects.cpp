#include "caloObjects.h"

// Unpacker function
// This function  unpacks the input set of links into fixed-grid Super Towers

void unpackToSuperTowers(const ap_uint<LINK_WIDTH> link_in[11],
                         hftower superTowers[(TOWERS_IN_ETA/3)+2][((TOWERS_IN_PHI/3)/N_SECTORS_ST)+2])
{
	 #pragma HLS ARRAY_PARTITION variable=superTowers type=complete dim=0
     #pragma HLS ARRAY_PARTITION variable=link_in type=complete

	hftower HFTowers[TOWERS_IN_ETA][11*4];
     #pragma HLS ARRAY_PARTITION variable=HFTowers type=complete dim=0

    unpackToSuperTowers_1_1:
    for(loop link=0; link<11; link++)
    {

        hftower HFTowers_temp[TOWERS_IN_ETA][TOWERS_IN_PHI/N_INPUT_LINKS];
        #pragma HLS ARRAY_PARTITION variable=HFTowers_temp type=complete dim=0
        processInputLink(link_in[link], HFTowers_temp);

        unpackToSuperTowers_1_2:
        for(loop eta=0; eta<TOWERS_IN_ETA; eta++)
        {
        	unpackToSuperTowers_1_3:
            for(loop phi=0; phi<TOWERS_IN_PHI/N_INPUT_LINKS; phi++)
            {
                HFTowers[eta][(link*(TOWERS_IN_PHI/N_INPUT_LINKS))+phi] = HFTowers_temp[eta][phi];
            }
        }
    }

	unpackToSuperTowers_2_2:
	for (loop eta = 0; eta < 4; ++eta){
		unpackToSuperTowers_2_3:
		for (loop phi = 0; phi < 13; ++phi){

			unpackToSuperTowers_2_4:
			for (loop m = 0; m < 3; ++m){
				unpackToSuperTowers_2_5:
				for (loop n = 0; n < 3; ++n){
					superTowers[eta+1][phi].energy += HFTowers[eta * 3 + m][(phi * 3 + n)+1].energy;
				}
			}
			superTowers[eta+1][phi].eta = eta;
			superTowers[eta+1][phi].phi = phi;
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
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    if(DEBUG_LEVEL > 2)
    {
        std::cout<<"HF Super Towers : "<<"\n";
        for(loop eta=0; eta<(TOWERS_IN_ETA/3) + 2; eta++)
        {
            for(loop phi=0; phi<(TOWERS_IN_PHI/6) + 2; phi++)
            {
                std::cout<< "["<<eta<<","<<phi<<"] :" << superTowers[eta][phi].energy;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    if(DEBUG_LEVEL > 0)
    {
        for(loop eta=0; eta<(TOWERS_IN_ETA/3) + 2; eta++)
        {
            std::cout<<"@@HFSuperTowers | "<<eta<<" | ";
            for(loop phi=0; phi<(TOWERS_IN_PHI/6) + 2; phi++)
            {
                std::cout<<std::setw(3)<<superTowers[eta][phi].energy<<" |";
            }
            std::cout << std::endl;
        }
    }

#endif
}

void findMaxEnergySuperTowerInPhi(const hftower PhiTowers[24],
							hftower& phiC) {
#pragma HLS INLINE off
    hftower tempArray[32];
     #pragma HLS ARRAY_PARTITION variable=tempArray type=complete

    findMaxEnergySuperTowerInPhi_1_1:
    for(loop i=0; i<24; i++){
        tempArray[i] = PhiTowers[i];
    }

    findMaxEnergySuperTowerInPhi_2_1:
    for(loop i=32; i>1; i=(i/2))
    {
    	findMaxEnergySuperTowerInPhi_2_2:
        for(loop j=0; j < i/2; j++)
        {
        	tempArray[j] = bestOf2(tempArray[j*2], tempArray[(j*2) + 1]);
        }
    }

    phiC = tempArray[0];
}


void findMaxEnergySuperTowerInEta(const hftower EtaTowers[4], hftower& etaC) {
#pragma HLS INLINE off
    hftower tempArray[4];
     #pragma HLS ARRAY_PARTITION variable=tempArray type=complete

    findMaxEnergySuperTowerInEta_1_1:
    for(loop i=0; i<4; i++) {
        tempArray[i] = EtaTowers[i];
    }

    findMaxEnergySuperTowerInEta_2_1:
    for(loop i=4; i>1; i=(i/2))
    {
    	findMaxEnergySuperTowerInEta_2_2:
        for(loop j=0; j < i/2; j++)
        {
        	tempArray[j] = bestOf2(tempArray[j*2], tempArray[(j*2) + 1]);
        }
    }

    etaC = tempArray[0];

}

void findMaxEnergySuperTower(const hftower HFRegion[(TOWERS_IN_ETA/3)+2][((TOWERS_IN_PHI/3)/N_SECTORS_ST)+2],
                        ap_uint<5>& etaC,
                        ap_uint<8>& phiC)
{
//#pragma HLS INLINE
     #pragma HLS ARRAY_PARTITION variable=HFRegion type=complete dim=0

	hftower towersPhi[24];
     #pragma HLS ARRAY_PARTITION variable=towersPhi type=complete
	hftower tempPhi;
	findMaxEnergySuperTower_1_1:
    for(ap_uint<5> phi = 1; phi < 13; phi++)
    {
    	hftower towersEta[4];
         #pragma HLS ARRAY_PARTITION variable=towersEta type=complete
    	findMaxEnergySuperTower_1_2:
        for(loop eta = 1; eta < 5; eta++)
        {
            towersEta[eta-1] = HFRegion[eta][phi];
        }
        findMaxEnergySuperTowerInEta(towersEta, towersPhi[phi-1]);
    }

    findMaxEnergySuperTowerInPhi(towersPhi, tempPhi);

    etaC = tempPhi.eta+1;
    phiC = tempPhi.phi;
#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL >7 )
    {
        std::cout<<"  Center found at : "<<etaC<<" , "<<phiC<<"\n" ;
    }
#endif
}


void formJetsAndZeroOut(hftower superTowers[(TOWERS_IN_ETA/3)+2][((TOWERS_IN_PHI/3)/N_SECTORS_ST)+2],
                        ap_uint<10> etaC,
                        ap_uint<10> phiC,
                        ap_uint<18>& etaSum) {
#pragma HLS ARRAY_PARTITION variable=superTowers type=complete dim=0

    ap_uint<1> mask[(TOWERS_IN_ETA/3)+2][((TOWERS_IN_PHI/3)/N_SECTORS_ST)+2];
#pragma HLS ARRAY_PARTITION variable=mask type=complete dim=0
    ap_uint<1> zero[(TOWERS_IN_ETA/3)+2][((TOWERS_IN_PHI/3)/N_SECTORS_ST)+2];
#pragma HLS ARRAY_PARTITION variable=zero type=complete dim=0


    for (ap_uint<10> eta = 0; eta < (TOWERS_IN_ETA/3)+2; eta++) {
    	for (ap_uint<10> phi = 0; phi < ((TOWERS_IN_PHI/3)/N_SECTORS_ST)+2; phi++) {
    		mask[eta][phi] = 0;
    		zero[eta][phi] = 0;
    	}
    }

    formJetsAndZeroOut_1_1:
    for (ap_uint<10> i = 0; i < (TOWERS_IN_ETA / 3) + 2; i++) {
        if (i + 1 >= etaC && i <= etaC + 1) {
        	formJetsAndZeroOut_1_2:
            for (ap_uint<10> j = 0; j < (TOWERS_IN_PHI / 6) + 2; j++) {
                if (j + 1 >= phiC && j <= phiC + 1) {
                    mask[i][j] = 1;
                }
            }
        }
    }

    formJetsAndZeroOut_2_1:
    for (ap_uint<10> i = 0; i < (TOWERS_IN_ETA / 3) + 2; i++) {
    	formJetsAndZeroOut_2_2:
        for (ap_uint<10> j = 0; j < (TOWERS_IN_PHI / 6) + 2; j++) {
            etaSum += superTowers[i][j].energy * mask[i][j];
        }
    }

    formJetsAndZeroOut_3_1:
    for (ap_uint<10> i = 0; i < (TOWERS_IN_ETA / 3) + 2; i++) {
    	formJetsAndZeroOut_3_2:
        for (ap_uint<10> j = 0; j < (TOWERS_IN_PHI / 6) + 2; j++) {
            zero[i][j] = (ap_uint<1>)(1 - mask[i][j]);
        }
    }

    formJetsAndZeroOut_4_1:
    for (ap_uint<10> i = 0; i < (TOWERS_IN_ETA / 3) + 2; i++) {
    	formJetsAndZeroOut_4_2:
        for (ap_uint<10> j = 0; j < (TOWERS_IN_PHI / 6) + 2; j++) {
            superTowers[i][j].energy *= zero[i][j];
        }
    }
}

void selectTaus(const jets Jet[9], jets Taus[9])
{
#pragma HLS INLINE off
 #pragma HLS ARRAY_PARTITION variable=Jet type=complete
 #pragma HLS ARRAY_PARTITION variable=Taus type=complete

    selectTaus_1_1:
    for(loop idx = 0; idx < 9; idx++)
    {
        if((Jet[idx].seedET * 10) >= (Jet[idx].ET * 7))
        {
            Taus[idx] = Jet[idx];
        }
    }

#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 2)
    {
        std::cout<< "Taus : " << std::endl;
        for(loop cluster=0; cluster<9; cluster++)
        {
            std::cout<<"Tau " << cluster << " E = "<< Taus[cluster].ET << " Seed = " << Taus[cluster].seedET
                     <<", center = ("<< Taus[cluster].Eta <<","<< Taus[cluster].Phi << ")\n";
        }
        std::cout<<std::endl;
    }
#endif

}

void makeJets(const ap_uint<LINK_WIDTH> link_in[11],
		const ap_uint<1> sector,
		jets Jet[5])
{
#pragma HLS INLINE off
 #pragma HLS ARRAY_PARTITION variable=link_in type=complete
 #pragma HLS ARRAY_PARTITION variable=Jet type=complete

	hftower superTowers[(TOWERS_IN_ETA/3)+2][((TOWERS_IN_PHI/3)/N_SECTORS_ST)+2];
 #pragma HLS ARRAY_PARTITION variable=superTowers type=complete dim=0

	unpackToSuperTowers(link_in, superTowers);

    for(loop idx = 0; idx < 5; idx++)
    {
        ap_uint<8> phiC = 1;
        ap_uint<5> etaC = 1;
        ap_uint<18> etSum = 0;

        findMaxEnergySuperTower(superTowers, etaC, phiC);

        Jet[idx].seedET = superTowers[etaC][phiC].energy;
        formJetsAndZeroOut(superTowers, etaC, phiC, etSum);

//        if(etSum ==0 )
//        {
//            etaC=0;
//            phiC=0;
//        }
        Jet[idx].ET  = etSum;
        Jet[idx].Eta = etaC-1;
        Jet[idx].Phi = phiC-1 + (12*sector);

    }

#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 2)
    {
        std::cout<< "Jets : " << std::endl;
        for(loop cluster=0; cluster<5; cluster++)
        {
            std::cout<<"Jet " << cluster << " E = "<< Jet[cluster].ET << " Seed = " << Jet[cluster].seedET
                     <<", center = ("<< Jet[cluster].Eta <<","<< Jet[cluster].Phi << ")\n";
        }
        std::cout<<std::endl;
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

    Exy_1_1:
    for(loop j = 0; j < TOWERS_IN_PHI/2; j=j+2) {
    	Exy_1_2:
        for(loop i = 0; i < TOWERS_IN_ETA - 2; i++) {

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


void swap_1(const jets Jets_in[10], jets Jets_out[10]) {
#pragma HLS INLINE off
#pragma HLS PIPELINE
	jetGreaterSmaller res;

    res = AscendDescend(Jets_in[0], Jets_in[1]);
    Jets_out[0] = res.greater;
    Jets_out[1] = res.smaller;

    res = AscendDescend(Jets_in[2], Jets_in[3]);
    Jets_out[2] = res.greater;
    Jets_out[3] = res.smaller;

    res = AscendDescend(Jets_in[4], Jets_in[5]);
    Jets_out[4] = res.greater;
    Jets_out[5] = res.smaller;

    res = AscendDescend(Jets_in[6], Jets_in[7]);
    Jets_out[6] = res.greater;
    Jets_out[7] = res.smaller;

    res = AscendDescend(Jets_in[8], Jets_in[9]);
    Jets_out[8] = res.greater;
    Jets_out[9] = res.smaller;
}

void swap_2(const jets Jets_in[10], jets Jets_out[10]) {
#pragma HLS INLINE off
	jetGreaterSmaller res;

    Jets_out[0] = Jets_in[0];
    res = AscendDescend(Jets_in[1], Jets_in[2]);
    Jets_out[1] = res.greater;
    Jets_out[2] = res.smaller;

    res = AscendDescend(Jets_in[3], Jets_in[4]);
    Jets_out[3] = res.greater;
    Jets_out[4] = res.smaller;

    res = AscendDescend(Jets_in[5], Jets_in[6]);
    Jets_out[5] = res.greater;
    Jets_out[6] = res.smaller;

    res = AscendDescend(Jets_in[7], Jets_in[8]);
    Jets_out[7] = res.greater;
    Jets_out[8] = res.smaller;

    Jets_out[9] = Jets_in[9];
}

void sortJets(const jets Jets_in[N_SECTORS_ST][5], jets Jets_out[9]) {
#pragma HLS INLINE off

#pragma HLS ARRAY_PARTITION variable=Jets_in type=complete dim=0
#pragma HLS ARRAY_PARTITION variable=Jets_out type=complete

	jets __Jets_in[10];
#pragma HLS ARRAY_PARTITION variable=__Jets_in type=complete

	for(int sector = 0 ; sector < N_SECTORS_ST; sector++){
		for(int jet = 0 ; jet < 5; jet++){
			__Jets_in[5*sector + jet] = Jets_in[sector][jet];
		}
	}


    jets temp_1[10], temp_2[10], temp_3[10], temp_4[10], temp_5[10], temp_6[10];
    jets temp_7[10], temp_8[10], temp_9[10], temp_10[10];

    swap_1(__Jets_in, temp_1);
    swap_2(temp_1, temp_2);
    swap_1(temp_2, temp_3);
    swap_2(temp_3, temp_4);
    swap_1(temp_4, temp_5);
    swap_2(temp_5, temp_6);
    swap_1(temp_6, temp_7);
    swap_2(temp_7, temp_8);
    swap_1(temp_8, temp_9);
    swap_2(temp_9, temp_10);

    for (int i = 0; i < 9; i++) {
        Jets_out[i] = temp_10[i];
    }
}


void makeCaloObjects(const ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS],
                     jets Jets[9], jets Taus[9], ap_fixed<32,16>& Ex, ap_fixed<32,16>& Ey, ap_uint<12>& HT)
{
 #pragma HLS ARRAY_PARTITION variable=link_in type=complete

	ap_uint<LINK_WIDTH> __link_in[N_INPUT_LINKS];
 #pragma HLS ARRAY_PARTITION variable=__link_in type=complete
	for(int link = 0 ; link < N_INPUT_LINKS ; link++){
		__link_in[link] = link_in[link];
	}

	ap_uint<LINK_WIDTH> linksInSector[N_SECTORS_ST][11];
 #pragma HLS ARRAY_PARTITION variable=linksInSector type=complete dim=0
	for(int sector = 0 ; sector < N_SECTORS_ST; sector++){
		for(int link = 0 ; link < N_INPUT_LINKS/N_SECTORS_ST; link++){
			linksInSector[sector][link+1] = link_in[9*sector + link];
		}
	}

	linksInSector[0][0] = link_in[17];
	linksInSector[0][10] = link_in[9];

	linksInSector[1][0] = link_in[8];
	linksInSector[1][10] = link_in[0];

	hftower superTowerGrid[2][(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/6)+2];
 #pragma HLS ARRAY_PARTITION variable=superTowerGrid type=complete dim=0

    jets _Jets[N_SECTORS_ST][5];
 #pragma HLS ARRAY_PARTITION variable=_Jets type=complete dim=0

    jets sorted_Jets[9];
 #pragma HLS ARRAY_PARTITION variable=sorted_Jets type=complete

    jets _Taus[9];
 #pragma HLS ARRAY_PARTITION variable=_Taus type=complete

    ap_uint<12> temp_HT;


    Exy(__link_in, Ex, Ey);

    for(int sector = 0 ; sector < N_SECTORS_ST; sector++){
    	makeJets(linksInSector[sector], sector, _Jets[sector]);
    }

    sortJets(_Jets, sorted_Jets);

#ifndef __SYNTHESIS__
    if(DEBUG_LEVEL > 2)
    {
        std::cout<< "Sorted Jets : " << std::endl;
        for(loop cluster=0; cluster<9; cluster++)
        {
            std::cout<<"Jet " << cluster << " E = "<< sorted_Jets[cluster].ET << " Seed = " << sorted_Jets[cluster].seedET
                     <<", center = ("<< sorted_Jets[cluster].Eta <<","<< sorted_Jets[cluster].Phi << ")\n";
        }
        std::cout<<std::endl;
    }
#endif

    selectTaus(sorted_Jets, _Taus);

    makeCaloObjects_2_1:
    for(loop idx = 0; idx < 9; idx++)
    {
        Jets[idx] = sorted_Jets[idx];
        temp_HT += sorted_Jets[idx].ET;
    }

    makeCaloObjects_3_1:
    for(loop idx = 0; idx < 9; idx++)
    {
        Taus[idx] = _Taus[idx];
    }

    HT = temp_HT;
}





