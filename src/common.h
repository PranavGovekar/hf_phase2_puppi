#ifndef COMMON_H
#define COMMON_H

#include "ap_int.h"
#include "pfCluster.h"
#include "hfTowers.h"
#include "jets.h"


#define TOWERS_IN_ETA 12
#define TOWERS_IN_PHI 72
#define EXTRA_IN_PHI 4
#define EXTRA_IN_ETA 1
#define MIN_CLUSTER_SEED_ENERGY 1

#define LINK_WIDTH 220

#define N_INPUT_LINKS  18
#define N_OUTPUT_LINKS  6

#define DEBUG_LEVEL 1000

#define N_PF_LINK 8
#define N_PUPPI_LINK 8
#define N_SECTORS 18
#define N_PF 48
#define N_EXTRA (NCALO - NNEUTRALS)
#define N_SECTORS_PF 6

#define N_PF_CLUSTERS 4
#define NTOWER_IN_ETA_PER_SECTOR (TOWERS_IN_ETA + EXTRA_IN_ETA*2)
#define NTOWER_IN_PHI_PER_SECTOR ((TOWERS_IN_PHI/N_SECTORS) + EXTRA_IN_PHI*2)
#define LINKS_PER_REGION ((N_INPUT_LINKS/N_SECTORS) + 2)
#define N_SORT_ELEMENTS 16

class GreaterSmaller
{
public:
	PFcluster greater, smaller;
};

class jetGreaterSmaller
{
public:
	jets greater, smaller;
};

typedef ap_uint<10> loop;

void processInputLink(ap_uint<LINK_WIDTH> link_in,
                      hftower towerGrid[TOWERS_IN_ETA][TOWERS_IN_PHI/N_INPUT_LINKS]);

template <int SIZE, int N_POWER_2>
void findMaxEnergyTowerInArray(const hftower Towers[SIZE], ap_uint<10>& center);

hftower bestOf2(const hftower& ecaltp0, const hftower& ecaltp1);
PFcluster  bestOf2(const PFcluster& ecaltp0, const PFcluster& ecaltp1);

GreaterSmaller AscendDescend(const PFcluster &x, const PFcluster &y);
jetGreaterSmaller AscendDescend(const jets &x, const jets &y);

#endif
