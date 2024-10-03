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

//
//void findMaxEnergySuperTower(const hftower HFRegion[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2],
//                             ap_uint<10>& etaC,
//                             ap_uint<10>& phiC)
//{
//
//    #pragma HLS PIPELINE
//    #pragma HLS ARRAY_PARTITION variable=HFRegion complete dim=0
//
//    hftower towersPhi[4];
//    ap_uint<10> PhiCenters[4];
//    ap_uint<10> EtaCenter;
//
//    #pragma HLS ARRAY_PARTITION variable=towersPhi complete dim=0
//    #pragma HLS ARRAY_PARTITION variable=PhiCenters complete dim=0
//
//    for(ap_uint<5> eta = 1; eta < 5; eta++)
//    {
//        //std::cout<<"   ETA SCAN !! "<<eta<<" \n";
//        findMaxEnergyTowerInArray<24,32>(&HFRegion[eta][1], PhiCenters[eta-1]);
//        PhiCenters[eta-1] += 1;
//        towersPhi[eta-1] = HFRegion[eta][PhiCenters[eta-1]];
//        //std::cout<<" eta = "<<eta<<" | ";
//        //for(int i =0 ; i < (TOWERS_IN_PHI/3)+2  ;  i++)
//        //{
//        //    std::cout<<HFRegion[eta][i].energy<<",";
//        //}
//        //std::cout<<"   | [ "<<towersPhi[PhiCenters[eta-1]].energy<<" / "<< PhiCenters[eta-1] <<" ] \n";
//    }
//
//    findMaxEnergyTowerInArray<4,4>(towersPhi, EtaCenter);
//    //std::cout<<"   Phi SCAN !! "<< EtaCenter <<"\n";
//
//    phiC = PhiCenters[EtaCenter];
//    etaC = EtaCenter+1;
//}

void findMaxEnergySuperTowerInPhi(const hftower HFRegion[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2],
                             const ap_uint<5> etaCenters[24], ap_uint<8>& phiC)
{
    #pragma HLS PIPELINE II=1
    hftower tempArray[32];
    ap_uint<5> index[32];

    #pragma HLS ARRAY_PARTITION variable=tempArray complete dim=0
    #pragma HLS ARRAY_PARTITION variable=index complete dim=0

    for(loop i=0; i<24; i++)
    {
    	index[i] = i;
        tempArray[i] = HFRegion[etaCenters[i]][i+1];
    }

    for(loop i=32; i>1; i=(i/2))
    {
        for(loop j=0; j < i/2; j++)
        {
            if(tempArray[index[j*2]].energy < tempArray[index[(j*2) + 1]].energy)
            {
                index[j] = index[(j*2)+1];
            }
            else
            {
                index[j] = index[j*2];
            }
        }
    }

    phiC = index[0];
}


void findMaxEnergySuperTowerInEta(const hftower EtaTowers[4], ap_uint<5>& etaC)
{
    #pragma HLS PIPELINE II=1
    hftower tempArray[4];
    ap_uint<5> index[4];

    #pragma HLS ARRAY_PARTITION variable=tempArray complete dim=0
    #pragma HLS ARRAY_PARTITION variable=index complete dim=0

    for(loop i=0; i<4; i++)
    {
        tempArray[i] = EtaTowers[i];
        index[i] = i;
    }


    for(loop i=4; i>1; i=(i/2))
    {
        for(loop j=0; j < i/2; j++)
        {
            if(tempArray[index[j*2]].energy < tempArray[index[(j*2) + 1]].energy)
            {
                index[j] = index[(j*2)+1];
            }
            else
            {
                index[j] = index[j*2];
            }
        }
    }

    etaC = index[0];
}


void findMaxEnergySuperTower(const hftower HFRegion[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2],
                        ap_uint<5>& etaC,
                        ap_uint<8>& phiC)
{
    #pragma HLS ARRAY_PARTITION variable=HFRegion complete dim=0

    ap_uint<5> towersPhi[24];
    #pragma HLS ARRAY_PARTITION variable=towersPhi complete dim=0
    ap_uint<8> tempPhi;
    for(ap_uint<5> phi = 1; phi < 25; phi++)
    {
        hftower towersEta[12];
        #pragma HLS ARRAY_PARTITION variable=towersEta complete dim=0
        for(loop eta = 1; eta < 5; eta++)
        {
            towersEta[eta] = HFRegion[eta][phi];
        }
        findMaxEnergySuperTowerInEta(towersEta, towersPhi[phi-1]);
    }

    findMaxEnergySuperTowerInPhi(HFRegion, towersPhi, tempPhi);

    etaC = towersPhi[tempPhi];
    phiC = tempPhi + 1;
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
                        ap_uint<12>& etaSum)
{

    #pragma HLS ARRAY_PARTITION variable=superTowers complete dim=0

    etaSum = superTowers[etaC][phiC].energy +
             superTowers[etaC+1][phiC].energy +
             superTowers[etaC-1][phiC].energy +
             superTowers[etaC][phiC+1].energy +
             superTowers[etaC][phiC-1].energy +
             superTowers[etaC+1][phiC+1].energy +
             superTowers[etaC-1][phiC+1].energy +
             superTowers[etaC+1][phiC-1].energy +
             superTowers[etaC-1][phiC-1].energy;

    superTowers[etaC][phiC].energy = 0;
    superTowers[etaC+1][phiC].energy = 0;
    superTowers[etaC-1][phiC].energy = 0;
    superTowers[etaC][phiC+1].energy = 0;
    superTowers[etaC][phiC-1].energy = 0;
    superTowers[etaC+1][phiC+1].energy = 0;
    superTowers[etaC-1][phiC+1].energy = 0;
    superTowers[etaC+1][phiC-1].energy = 0;
    superTowers[etaC-1][phiC-1].energy = 0;
}

void selectTaus(const jets Jet[9], jets Taus[9])
{

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

    for(loop idx = 0; idx < 9; idx++)
    {
        ap_uint<8> phiC = 1;
        ap_uint<5> etaC = 1;
        ap_uint<12> etSum = 0;

        findMaxEnergySuperTower(superTowers, etaC, phiC);

        Jet[idx].seedET = superTowers[etaC][phiC].energy;
        formJetsAndZeroOut(superTowers, etaC, phiC, etSum);

        if(etSum ==0 )
        {
            etaC=0 ;
            phiC=0;
        }
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

void Exy (hftower superTowers[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2],
		ap_uint<32>& HT, ap_fixed<32,16>& MET){

	ap_fixed<32,16> sin_LUT[24] = {0.1305,0.3827,0.6088,0.7934,0.9239,0.9914,0.9914,0.9239,0.7934,0.6088,0.3827,0.1305,
			-0.1305,-0.3827,-0.6088,-0.7934,-0.9239,-0.9914,-0.9914,-0.9239,-0.7934,-0.6088,-0.3827,-0.1305};

//	[(phi + 6) % 24];

	ap_fixed<32,16> Ex;
	ap_fixed<32,16> Ey;
	ap_uint<32> temp_HT;
	for(loop eta = 1; eta < 5; eta++) {
		for(loop phi = 1; phi < 25; phi++) {
			Ex += (superTowers[eta][phi].energy) * sin_LUT[((phi + 6) % 24)-1];
			Ey += (superTowers[eta][phi].energy) * sin_LUT[phi-1];
			if(superTowers[eta][phi].energy > 20){
				temp_HT += superTowers[eta][phi].energy;
			}
			std::cout  << "eta : " << eta << " phi : " << phi << " [Ex : " << Ex << " Ey : " << Ey << "]\n";
		}
	}

	ap_fixed<32,16> MET2 = (Ex*Ex) + (Ey*Ey);
	std::cout << "MET^2 : " << MET2 << endl;
	MET = hls::sqrt((Ex*Ex) + (Ey*Ey));
	HT = temp_HT;

#ifndef __SYNTHESIS__
	for(loop eta = 1; eta < 5; eta++) {
		for(loop phi = 1; phi < 25; phi++) {
			std::cout << "\t" << superTowers[eta][phi].energy;
		}
		std::cout<<endl;
	}
	std::cout<<endl;
	std::cout << "Ex : " << Ex << "\nEy : " << Ey << endl;
	std::cout << "MET : " << MET << "\nHT : " << HT << endl;
	std::cout<<endl;
#endif
}

void makeCaloObjects(const ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS],
                     jets Jets[9], jets Taus[9] )
{

    #pragma HLS PIPELINE

    hftower superTowerGrid[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2];
    jets _Jets[9];
    jets _Taus[9];

	ap_uint<32> HT;
	ap_fixed<32,16> MET;

    unpackToSuperTowers(link_in, superTowerGrid);

    Exy(superTowerGrid, HT, MET);

    makeJets(superTowerGrid, _Jets);

    selectTaus(_Jets, _Taus);

    for(loop idx = 0; idx < 9; idx++)
    {
        Jets[idx] = _Jets[idx];
    }

    for(loop idx = 0; idx < 9; idx++)
    {
        Taus[idx] = _Taus[idx];
    }

}





