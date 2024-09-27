#ifndef ALGO_TOPIP1_H
#define ALGO_TOPIP1_H

#include "firmware/linpuppi.h"
#include <iostream>
#include "ap_int.h"
#include <algorithm>
#include <utility>
#include <stdint.h>

#include "layer1_objs.h"
#include "pf.h"
#include "puppi.h"
#include "layer1_multiplicities.h"
#include "firmware/linpuppi_bits.h"

#include "bitonicSort16.h"
#include "pfCluster.h"

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
#define N_SORT_ELEMENTS N

using namespace std;
typedef ap_uint<10> loop;



void algo_topIP1(ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], ap_uint<576> link_out[6]);

//void processInputLink( ap_uint<LINK_WIDTH> link_in,
//						hftower HFTowers[TOWERS_IN_ETA][TOWERS_IN_PHI/N_INPUT_LINKS]);
//
//void makeRegion(const ap_uint<LINK_WIDTH> link_in[LINKS_PER_REGION],
//						hftower HFRegions[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR]);
//
//void findMaxEnergyTowerInEta(const hftower EtaTowers[TOWERS_IN_ETA], ap_uint<5>& etaC);
//
//void findMaxEnergyTowerInPhi(const hftower HFRegion[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR],
//								const ap_uint<5> etaCenters[6],
//								ap_uint<8>& phiC);
//void findMaxEnergyTower(const hftower HFRegion[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR],
//				ap_uint<5>& etaC,
//				ap_uint<8>& phiC);
//
//void formClusterAndZeroOut(hftower HFRegion[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR],
//		ap_uint<5> etaC,
//		ap_uint<8> phiC,
//		ap_uint<12>& etaSum);
//
//void clustrizer (const ap_uint<LINK_WIDTH> regionLinks[LINKS_PER_REGION],
//		PFcluster Clusters[N_PF_CLUSTERS]);
//
//void packer(PFcluster Clusters[N_PF_CLUSTERS], const ap_uint<576>& link_out_sector, const ap_uint<7> sector);


#endif
