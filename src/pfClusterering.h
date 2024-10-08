#ifndef PUPPI_H
#define PUPPI_H

#define REG_HF

#include "layer1_objs.h"
#include "puppi.h"
#include "layer1_multiplicities.h"
#include "firmware/linpuppi.h"
#include "firmware/linpuppi_bits.h"
#include "pf.h"

#include "common.h"

void makeRegion(const ap_uint<LINK_WIDTH> link_in[LINKS_PER_REGION],
                hftower HFRegions[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR]);

void findMaxEnergyTowerInPhi(const hftower HFRegion[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR],
                             const ap_uint<5> etaCenters[6], ap_uint<8>& phiC);

void findMaxEnergyTowerInEta(const hftower EtaTowers[TOWERS_IN_ETA], ap_uint<5>& etaC);

void findMaxEnergyTower(const hftower HFRegion[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR],
                        ap_uint<5>& etaC, ap_uint<8>& phiC);

void formClusterAndZeroOut(hftower HFRegion[NTOWER_IN_ETA_PER_SECTOR][NTOWER_IN_PHI_PER_SECTOR],
                           ap_uint<5> etaC, ap_uint<8> phiC, ap_uint<12>& etaSum);

void makeCaloClusters(const ap_uint<LINK_WIDTH> regionLinks[LINKS_PER_REGION], PFcluster Clusters[N_PF_CLUSTERS]);

void packer(PFcluster Clusters[N_PF_CLUSTERS], const ap_uint<576>& link_out_sector, const ap_uint<7> sector);

void makePFClusters(const ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS], PFcluster ClustersOut[N_OUTPUT_LINKS][N_SORT_ELEMENTS]);


#endif
