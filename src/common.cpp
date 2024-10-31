#include "common.h"

// Unpacker function for a wedge
// This function takes one link and makes 12x4 grid of hftower Array
void processInputLink( ap_uint<LINK_WIDTH> link_in,
                       hftower towerGrid[TOWERS_IN_ETA][TOWERS_IN_PHI/N_INPUT_LINKS])
{

    #pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
    #pragma HLS ARRAY_PARTITION variable=towerGrid complete dim=0

    for(loop j=0; j<(TOWERS_IN_PHI/N_INPUT_LINKS)/2; j=j+2)
    {
        for(loop i=0; i<TOWERS_IN_ETA-1; i++)
        {
            hftower halfEnergy = hftower(link_in.range((i*10)+9, i*10));
            halfEnergy.energy = halfEnergy.energy << 1;

            towerGrid[i][j] = halfEnergy;
            towerGrid[i][j+1] = halfEnergy;

            hftower halfEnergy_0 = hftower(link_in.range((i*10+110)+9, i*10+110));
            halfEnergy_0.energy = halfEnergy_0.energy << 1;

            towerGrid[i][j+2] = halfEnergy_0;
            towerGrid[i][j+3] = halfEnergy_0;
        }
    }

    hftower A10 = towerGrid[TOWERS_IN_ETA-2][0] ;
    hftower B10 = towerGrid[TOWERS_IN_ETA-2][2] ;

    A10.energy = A10.energy >> 1 ;
    B10.energy = B10.energy >> 1 ;

    towerGrid[TOWERS_IN_ETA-2][0] = A10;
    towerGrid[TOWERS_IN_ETA-2][1] = A10;
    towerGrid[TOWERS_IN_ETA-2][2] = A10;
    towerGrid[TOWERS_IN_ETA-2][3] = A10;

    towerGrid[TOWERS_IN_ETA-1][0] = B10;
    towerGrid[TOWERS_IN_ETA-1][1] = B10;
    towerGrid[TOWERS_IN_ETA-1][2] = B10;
    towerGrid[TOWERS_IN_ETA-1][3] = B10;


#ifndef __SYNTHESIS__

    if (DEBUG_LEVEL > 10)
    {
        for(loop eta=0; eta<TOWERS_IN_ETA; eta++)
        {
            for(loop phi=0; phi<TOWERS_IN_PHI/N_INPUT_LINKS; phi++)
            {
                std::cout<< "\t" << towerGrid[eta][phi].energy;
            }
            std::cout << endl;
        }
        std::cout << endl;
    }
#endif
}


template <int SIZE, int N_POWER_2>
void findMaxEnergyTowerInArray(const hftower Towers[SIZE], ap_uint<10>& center)
{
    #pragma HLS PIPELINE
//
//    // Calculate the next power of two
//    int N_POWER_2 = SIZE;
//    N_POWER_2 -= 1;
//    N_POWER_2 |= N_POWER_2 >> 1;
//    N_POWER_2 |= N_POWER_2 >> 2;
//    N_POWER_2 |= N_POWER_2 >> 4;
//    N_POWER_2 += 1;

    hftower tempArray[N_POWER_2];
    ap_int<10> index[N_POWER_2];

    #pragma HLS ARRAY_PARTITION variable=tempArray complete dim=0
    #pragma HLS ARRAY_PARTITION variable=index complete dim=0

    for(ap_uint<10> i=0; i<SIZE; i++)
    {
        tempArray[i] = Towers[i];
        index[i] = i;
    }
    for(ap_uint<10> i=SIZE; i<N_POWER_2; i++)
    {
        index[i] = SIZE;
    }


    for(ap_uint<10> i=N_POWER_2; i>1; i=(i/2))
    {
        //std::cout<<" -- > for i = "<<i<<"\n";
        //for(ap_uint<10> j=0; j < i; j++)
        //    std::cout<<setw(3)<<tempArray[index[j]].energy<<" , ";
        //std::cout<<std::endl;
        //for(ap_uint<10> j=0; j < i; j++)
        //    std::cout<<setw(3)<<index[j]<<" , ";
        //std::cout<<std::endl;
        for(ap_uint<10> j=0; j < i/2; j++)
        {
            //if(tempArray[index[j*2]].energy < tempArray[index[(j*2) + 1]].energy)
            if(tempArray[index[j]].energy < tempArray[index[j + i/2]].energy)
            {
                index[j] = index[j+i/2];
            }
            else
            {
                index[j] = index[j];
            }
        }
    }

    center = index[0];
}

