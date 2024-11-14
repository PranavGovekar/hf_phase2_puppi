#ifndef CALO_OBJS_H
#define CALO_OBJS_H

#include "common.h"
#include "jets.h"
#include <hls_math.h>

void unpackToSuperTowers(const ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS],
		hftower superTowers[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2]);

void findMaxEnergySuperTower(const hftower HFRegion[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2],
                        ap_uint<5>& etaC,
                        ap_uint<8>& phiC);

void formJetsAndZeroOut(hftower superTowers[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2],
                        ap_uint<10> etaC, ap_uint<10> phiC, ap_uint<18>& etaSum);

void selectTaus(const jets Jet[9], jets Taus[9]);

void makeJets(hftower superTowers[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2], jets Jet[9]);

//void Exy(hftower superTowers[(TOWERS_IN_ETA/3)+2][(TOWERS_IN_PHI/3)+2],
//         ap_uint<32>& HT, ap_fixed<32,16>& MET);

void makeCaloObjects(const ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS],
                     jets Jets[9], jets Taus[9], ap_fixed<32,16>& Ex, ap_fixed<32,16>& Ey, ap_uint<12>& HT);


#endif
