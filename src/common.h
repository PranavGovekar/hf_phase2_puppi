#ifndef COMMON_H
#define COMMON_H

#include "ap_int.h"
#include "bitonicSort16.h"
#include "pfCluster.h"
#include "hfTowers.h"

#define DEBUG_LEVEL 10

#define TOWERS_IN_ETA 12
#define TOWERS_IN_PHI 72
#define EXTRA_IN_PHI 4
#define EXTRA_IN_ETA 1
#define MIN_CLUSTER_SEED_ENERGY 1

#define LINK_WIDTH 220

#define N_INPUT_LINKS  18
#define N_OUTPUT_LINKS  6


#define N_PF_LINK 8
#define N_PUPPI_LINK 8
#define N_SECTORS 18
#define N_PF 48
#define N_EXTRA (NCALO - NNEUTRALS)


#define N_PF_CLUSTERS 4
#define NTOWER_IN_ETA_PER_SECTOR (TOWERS_IN_ETA + EXTRA_IN_ETA*2)
#define NTOWER_IN_PHI_PER_SECTOR ((TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2)
#define LINKS_PER_REGION ((N_INPUT_LINKS/N_SECTORS) + 2)
#define N_SORT_ELEMENTS 16


typedef ap_uint<10> loop;

void processInputLink(ap_uint<LINK_WIDTH> link_in,
                      hftower towerGrid[TOWERS_IN_ETA][TOWERS_IN_PHI/N_INPUT_LINKS]);

template <int SIZE, int N_POWER_2>
void findMaxEnergyTowerInArray(const hftower Towers[SIZE], ap_uint<10>& center);


#endif
