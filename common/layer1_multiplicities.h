#ifndef FIRMWARE_dataformats_layer1_multiplicities_h
#define FIRMWARE_dataformats_layer1_multiplicities_h

// DEFINE MULTIPLICITIES
#if defined(REG_HGCal)
#define NTRACK 30
#define NCALO 20
#define NCALOFWD 12 // FOR TMUX18 running
#define NMU 4
#define NPV 1
#define NSELCALO 20
#define NALLNEUTRALS NSELCALO
#define NPUPPIFINALSORTED 18
#define NEMCALO 10
#define NPHOTON NEMCALO 
// not used but must be there because used in header files
#define NNEUTRALS 1
// Configuration of EG algo follows
#define NEMCALO_EGIN 10
#define NTRACK_EGIN 10
// Composite EG params
#define NTRACK_PER_EMCALO_EGCOMP 3

#define NEM_EGOUT 5

//--------------------------------
#elif defined(REG_HGCalNoTK)
#define NCALO 12
#define NNEUTRALS 8
#define NALLNEUTRALS NCALO
#define NPUPPITOSORT NALLNEUTRALS
#define NPUPPIFINALSORTED NALLNEUTRALS
// dummy
#define NMU 1
#define NPV 1
#define NTRACK 1
#define NEMCALO 1
#define NPHOTON NEMCALO
#define NSELCALO 1
#define NTRACK_PER_EMCALO_EGCOMP 1

//--------------------------------
#elif defined(REG_HF)
#define NCALO 16
#define NNEUTRALS 8
#define NALLNEUTRALS NCALO
#define NPUPPIFINALSORTED NNEUTRALS 
// dummy
#define NMU 1
#define NPV 1
#define NTRACK 1
#define NEMCALO 1
#define NPHOTON NEMCALO
#define NSELCALO 1
#define NTRACK_PER_EMCALO_EGCOMP 1

//--------------------------------
#else // BARREL
#ifndef REG_Barrel
#ifndef CMSSW_GIT_HASH
#warning                                                                       \
    "No region defined, assuming it's barrel (#define REG_Barrel to suppress this)"
#endif
#endif

#define NTRACK 22
#define NCALO 15
#define NEMCALO 12
#define NMU 2
#define NPV 1
#define NPHOTON NEMCALO
#define NSELCALO NCALO
#define NALLNEUTRALS (NPHOTON + NSELCALO)
#define NNEUTRALS NALLNEUTRALS
#define NPUPPIFINALSORTED 18
// Configuration of EG algo follows
#define NEMCALO_EGIN 10
#define NTRACK_EGIN 13
#define NEM_EGOUT 10
//dummy
#define NTRACK_PER_EMCALO_EGCOMP 1

#define NREGIONS 18 //added for barrel sorter

#endif // region

// unique number for now
#define NL1_EGOUT 16

#if defined(BOARD_KU15P)
#define PACKING_DATA_SIZE 64
#define PACKING_NCHANN 42
#elif defined(BOARD_Serenity)
#define PACKING_DATA_SIZE 64
#define PACKING_NCHANN 120
#elif defined(BOARD_APD1)
#define PACKING_DATA_SIZE 64
#define PACKING_NCHANN 96
#else
#define PACKING_DATA_SIZE 64
#endif


namespace l1ct {

    template <int N> struct ct_log2_ceil {
      enum { value = ct_log2_ceil<(N / 2) + (N % 2)>::value + 1 };
    };
    template <> struct ct_log2_ceil<2> {
      enum { value = 1 };
    };
    template <> struct ct_log2_ceil<1> {
      enum { value = 0 };
    };

} // namespace


#endif
