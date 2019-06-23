
#include "algo_inputs_layer2_real_v3.h"

template<int NB>
ap_uint<NB> dr2_int_cap(etaphi_t eta1, etaphi_t phi1, etaphi_t eta2, etaphi_t phi2, ap_uint<NB> max) {
    #pragma HLS inline
    #pragma HLS PIPELINE
    //hardcode for etaphi size
    int tmpe = eta1-eta2;
    ap_uint<NB> deta = (tmpe > 0 ? tmpe : -tmpe);
    int tmpp = phi1-phi2;
    ap_uint<NB> dphi = (tmpp > 0 ? tmpp : -tmpp);
    int dr2 = max;
    if ((deta >> (NB/2))==0 && (dphi >> (NB/2))==0) {
        ap_uint<NB> deta2 = deta*deta;
        ap_uint<NB> dphi2 = dphi*dphi;
        dr2 = deta2 + dphi2;
    }
    return (dr2 < int(max) ? ap_uint<NB>(dr2) : max);
}
template<typename T, int NIn, int NOut>
void ptsort_hwopt_ind(T in[NIn], T out[NOut]) { 
    #pragma HLS PIPELINE
    T tmp[NOut];
    #pragma HLS ARRAY_PARTITION variable=tmp complete
    for (int iout = 0; iout < NOut; ++iout) {
        #pragma HLS unroll
        tmp[iout].hwPt = 0;
    }

    for (int it = 0; it < NIn; ++it) {
        for (int iout = NOut-1; iout >= 0; --iout) {
            if (tmp[iout].hwPt <= in[it].hwPt) {
                if (iout == 0 || tmp[iout-1].hwPt > in[it].hwPt) {
                    tmp[iout] = in[it];
                } else {
                    tmp[iout] = tmp[iout-1];
                }
            }
        }

    }
    for (int iout = 0; iout < NOut; ++iout) {
        out[iout] = tmp[iout];
    }

}
template<unsigned int N,unsigned int Offset> 
inline void mp7_unpack(axi_t data[], PFChargedObj emcalo[]) {
    #pragma HLS inline
    #pragma HLS PIPELINE
    for (unsigned int i = 0; i < N; ++i) {
      emcalo[i].hwPt       = data[Offset+i](15, 0);
      emcalo[i].hwEta      = data[Offset+i](24,16);
      emcalo[i].hwPhi      = data[Offset+i](34,25);
      emcalo[i].hwId       = data[Offset+i](38,35);
      emcalo[i].hwZ0       = data[Offset+i](48,39);
    }
}
template<unsigned int N,unsigned int Offset,unsigned int NR> 
inline void mp7_unpack2(axi_t data[], unsigned int iRegion,PFChargedObj emcalo[]) {
    #pragma HLS inline
    #pragma HLS PIPELINE
    for (unsigned int i = 0; i < N; ++i) {
      emcalo[iRegion*N+i].hwPt       = data[Offset+i](15, 0);
      emcalo[iRegion*N+i].hwEta      = data[Offset+i](24,16);
      emcalo[iRegion*N+i].hwPhi      = data[Offset+i](34,25);
      emcalo[iRegion*N+i].hwId       = data[Offset+i](38,35);
      emcalo[iRegion*N+i].hwZ0       = data[Offset+i](48,39);
    }
}
inline void mp7_unpack_seed(axi_t data, PFChargedObj &seed) {
    #pragma HLS inline
    #pragma HLS PIPELINE
    seed.hwPt          = data(15, 0);
    seed.hwEta         = data(24,16);
    seed.hwPhi         = data(34,25);
    seed.hwId          = data(38,35);
    seed.hwZ0          = data(48,39);
}
template<int N_TABLE>
void init_map_table(ap_int<8>  arr[N_TABLE], int iId){
  //4 blocks across*2 8 blocks                                                                                                                                                                                                                                                 
  for(int iregion = 0; iregion < N_TABLE; iregion++) {
    bool pUp    = true;
    bool pLeft  = true;
    bool pEdgeLeft  = false;
    bool pEdgeRight = false;
    bool pEdgeTop   = false;
    bool pEdgeBot   = false;
    int  pRow    = iregion/8;
    int  pX      = (iregion % 8)/2;
    int  pY      = iregion/16;
    int  pRegion = pY*4+pX;
    if(iregion % 2 == 1) pLeft        = false;
    if(pRow    % 2 == 1) pUp          = false;
    if(iregion % 8 == 0) pEdgeLeft    = true;
    if(iregion % 8 == 7) pEdgeRight   = true;
    if(iregion < 8     ) pEdgeTop     = true;
    if(iregion > N_TABLE-9) pEdgeBot  = true;
    if(pUp && pLeft) {
      if(iId == 0) arr[iregion] = pRegion;
      if(iId == 1) arr[iregion] = pRegion-4;
      if(iId == 2) arr[iregion] = pRegion-5;
      if(iId == 3) arr[iregion] = pRegion-1;
    }
    if(pUp && !pLeft) {
      if(iId == 0) arr[iregion] = pRegion;
      if(iId == 1) arr[iregion] = pRegion-4;
      if(iId == 2) arr[iregion] = pRegion-3;
      if(iId == 3) arr[iregion] = pRegion+1;
    }
    if(!pUp && pLeft) {
      if(iId == 0) arr[iregion] = pRegion;
      if(iId == 1) arr[iregion] = pRegion+4;
      if(iId == 2) arr[iregion] = pRegion+3;
      if(iId == 3) arr[iregion] = pRegion-1;
    }
    if(!pUp && !pLeft) {
      if(iId == 0) arr[iregion] = pRegion;
      if(iId == 1) arr[iregion] = pRegion+4;
      if(iId == 2) arr[iregion] = pRegion+5;
      if(iId == 3) arr[iregion] = pRegion+1;
    }
    if(pEdgeTop && iId == 1) arr[iregion] += N_TABLE/4;
    if(pEdgeTop && iId == 2) arr[iregion] += N_TABLE/4;
    if(pEdgeBot && iId == 1) arr[iregion] -= N_TABLE/4;
    if(pEdgeBot && iId == 2) arr[iregion] -= N_TABLE/4;
    if(pEdgeLeft || pEdgeRight) {
      if(iId == 2) arr[iregion] = 0;
      if(iId == 3) arr[iregion] = 0;
    }
  }
}

template<int N_TABLE>
void init_map_table0(ap_int<8>  arr[N_TABLE]){
  //4 blocks across*2 8 blocks                                                                                                                                                                                                                                                 
  for(int iregion = 0; iregion < N_TABLE; iregion++) {
    bool pUp    = true;
    bool pLeft  = true;
    bool pEdgeLeft  = false;
    bool pEdgeRight = false;
    bool pEdgeTop   = false;
    bool pEdgeBot   = false;
    int  pRow    = iregion/8;
    int  pX      = (iregion % 8)/2;
    int  pY      = iregion/16;
    int  pRegion = pY*4+pX;
    if(iregion % 2 == 1) pLeft        = false;
    if(pRow    % 2 == 1) pUp          = false;
    if(iregion % 8 == 0) pEdgeLeft    = true;
    if(iregion % 8 == 7) pEdgeRight   = true;
    if(iregion < 8     ) pEdgeTop     = true;
    if(iregion > N_TABLE-9) pEdgeBot  = true;
    if(pUp && pLeft) {
      arr[iregion] = pRegion;
    }
    if(pUp && !pLeft) {
      arr[iregion] = pRegion;
    }
    if(!pUp && pLeft) {
      arr[iregion] = pRegion;      
    }
    if(!pUp && !pLeft) {
      arr[iregion] = pRegion;
    }
  }
}
template<int N_TABLE>
void init_map_table1(ap_int<8>  arr[N_TABLE]){
  //4 blocks across*2 8 blocks                                                                                                                                                                                                                                                 
  for(int iregion = 0; iregion < N_TABLE; iregion++) {
    bool pUp    = true;
    bool pLeft  = true;
    bool pEdgeLeft  = false;
    bool pEdgeRight = false;
    bool pEdgeTop   = false;
    bool pEdgeBot   = false;
    int  pRow    = iregion/8;
    int  pX      = (iregion % 8)/2;
    int  pY      = iregion/16;
    int  pRegion = pY*4+pX;
    if(iregion % 2 == 1) pLeft        = false;
    if(pRow    % 2 == 1) pUp          = false;
    if(iregion % 8 == 0) pEdgeLeft    = true;
    if(iregion % 8 == 7) pEdgeRight   = true;
    if(iregion < 8     ) pEdgeTop     = true;
    if(iregion > N_TABLE-9) pEdgeBot  = true;
    if(pUp && pLeft) {
      arr[iregion] = pRegion-4;
    }
    if(pUp && !pLeft) {
      arr[iregion] = pRegion-4;
    }
    if(!pUp && pLeft) {
      arr[iregion] = pRegion+4;
    }
    if(!pUp && !pLeft) {
      arr[iregion] = pRegion+4;
    }
    if(pEdgeTop) arr[iregion] += N_TABLE/4;
    if(pEdgeBot) arr[iregion] -= N_TABLE/4;
  }
}

template<int N_TABLE>
void init_map_table2(ap_int<8>  arr[N_TABLE]){
  //4 blocks across*2 8 blocks                                                                                                                                                                                                                                                 
  for(int iregion = 0; iregion < N_TABLE; iregion++) {
    bool pUp    = true;
    bool pLeft  = true;
    bool pEdgeLeft  = false;
    bool pEdgeRight = false;
    bool pEdgeTop   = false;
    bool pEdgeBot   = false;
    int  pRow    = iregion/8;
    int  pX      = (iregion % 8)/2;
    int  pY      = iregion/16;
    int  pRegion = pY*4+pX;
    if(iregion % 2 == 1) pLeft        = false;
    if(pRow    % 2 == 1) pUp          = false;
    if(iregion % 8 == 0) pEdgeLeft    = true;
    if(iregion % 8 == 7) pEdgeRight   = true;
    if(iregion < 8     ) pEdgeTop     = true;
    if(iregion > N_TABLE-9) pEdgeBot  = true;
    if(pUp && pLeft) {
      arr[iregion] = pRegion-5;
    }
    if(pUp && !pLeft) {
      arr[iregion] = pRegion-3;
    }
    if(!pUp && pLeft) {
      arr[iregion] = pRegion+3;
    }
    if(!pUp && !pLeft) {
      arr[iregion] = pRegion+5;
    }
    if(pEdgeTop) arr[iregion] += N_TABLE/4;
    if(pEdgeBot) arr[iregion] -= N_TABLE/4;
    if(pEdgeLeft || pEdgeRight) {
      arr[iregion] = 0;
    }
  }
}

template<int N_TABLE>
void init_map_table3(ap_int<8>  arr[N_TABLE]){
  //4 blocks across*2 8 blocks                                                                                                                                                                                                                                                 
  for(int iregion = 0; iregion < N_TABLE; iregion++) {
    bool pUp    = true;
    bool pLeft  = true;
    bool pEdgeLeft  = false;
    bool pEdgeRight = false;
    bool pEdgeTop   = false;
    bool pEdgeBot   = false;
    int  pRow    = iregion/8;
    int  pX      = (iregion % 8)/2;
    int  pY      = iregion/16;
    int  pRegion = pY*4+pX;
    if(iregion % 2 == 1) pLeft        = false;
    if(pRow    % 2 == 1) pUp          = false;
    if(iregion % 8 == 0) pEdgeLeft    = true;
    if(iregion % 8 == 7) pEdgeRight   = true;
    if(iregion < 8     ) pEdgeTop     = true;
    if(iregion > N_TABLE-9) pEdgeBot  = true;
    if(pUp && pLeft) {
      arr[iregion] = pRegion-1;
    }
    if(pUp && !pLeft) {
      arr[iregion] = pRegion+1;
    }
    if(!pUp && pLeft) {
      arr[iregion] = pRegion-1;
    }
    if(!pUp && !pLeft) {
      arr[iregion] = pRegion+1;
    }
    if(pEdgeLeft || pEdgeRight) {
      arr[iregion] = 0;
    }
  }
}
void arraymap0(etaphi_t iEta,etaphi_t iPhi,ap_int<8> *arr) {
 #pragma HLS PIPELINE
#ifdef __HLS_SYN__
  bool initialized = false;
  ap_int<8> arr0[NREGIONS*4];
  #pragma HLS ARRAY_PARTITION variable=arr0 complete
#else
  static bool initialized = false;
  static ap_int<8> arr0[NREGIONS*4];
  #pragma HLS ARRAY_PARTITION variable=arr0 complete
#endif
  if (!initialized) {
    init_map_table0<NREGIONS*4>(arr0);
    initialized = true;
  }
  ap_int<8> pEta = ap_int<16>(iEta*0.0125) % 8;
  ap_int<8> pPhi = ap_int<16>(iPhi*0.0125) % (NREGIONS/4);
  *arr = arr0[pPhi*8+pEta];
}
void arraymap1(etaphi_t iEta,etaphi_t iPhi,ap_int<8> *arr) {
 #pragma HLS PIPELINE
#ifdef __HLS_SYN__
  bool initialized = false;
  ap_int<8> arr1[NREGIONS*4];
  #pragma HLS ARRAY_PARTITION variable=arr1 complete
#else
  static bool initialized = false;
  static ap_int<8> arr1[NREGIONS*4];
  #pragma HLS ARRAY_PARTITION variable=arr1 complete
#endif
  if (!initialized) {
    init_map_table1<NREGIONS*4>(arr1);
    initialized = true;
  }
  ap_int<8> pEta = ap_int<16>(iEta*0.0125) % 8;
  ap_int<8> pPhi = ap_int<16>(iPhi*0.0125) % (NREGIONS/4);
  *arr = arr1[pPhi*8+pEta];
}
void arraymap2(etaphi_t iEta,etaphi_t iPhi,ap_int<8> *arr) {
 #pragma HLS PIPELINE
#ifdef __HLS_SYN__
  bool initialized = false;
  ap_int<8> arr2[NREGIONS*4];
  #pragma HLS ARRAY_PARTITION variable=arr2 complete
#else
  static bool initialized = false;
  static ap_int<8> arr2[NREGIONS*4];
  #pragma HLS ARRAY_PARTITION variable=arr2 complete
#endif
  if (!initialized) {
    init_map_table2<NREGIONS*4>(arr2);
    initialized = true;
  }
  ap_int<8> pEta = ap_int<16>(iEta*0.0125) % 8;
  ap_int<8> pPhi = ap_int<16>(iPhi*0.0125) % (NREGIONS/4);
  *arr = arr2[pPhi*8+pEta];
}
void arraymap3(etaphi_t iEta,etaphi_t iPhi,ap_int<8> *arr) {
 #pragma HLS PIPELINE
#ifdef __HLS_SYN__
  bool initialized = false;
  ap_int<8> arr3[NREGIONS*4];
  #pragma HLS ARRAY_PARTITION variable=arr3 complete
#else
  static bool initialized = false;
  static ap_int<8> arr3[NREGIONS*4];
  #pragma HLS ARRAY_PARTITION variable=arr3 complete
#endif
  if (!initialized) {
    init_map_table3<NREGIONS*4>(arr3);
    initialized = true;
  }
  ap_int<8> pEta = ap_int<16>(iEta*0.0125) % 8;
  ap_int<8> pPhi = ap_int<16>(iPhi*0.0125) % (NREGIONS/4);
  *arr = arr3[pPhi*8+pEta];
}
template<unsigned int N0,unsigned int N1, unsigned int N2, unsigned int N3,unsigned int NOUT>
void tausort(PFChargedObj allparts[NOUT], PFChargedObj pfch[N0], PFChargedObj pfpho[N1], PFChargedObj pfne[N2], PFChargedObj pfmu[N3]) {
  #pragma HLS PIPELINE
  #pragma HLS ARRAY_PARTITION variable=allparts complete
  for(int i = 0; i < N0; i++) { 
    #pragma HLS UNROLL
    allparts[i] = pfch[i];
  }
  for(int i = 0; i < N1; i++) { 
    #pragma HLS UNROLL
    allparts[i+N0] = pfpho[i];
  }
  for(int i = 0; i < N2; i++) { 
    #pragma HLS UNROLL
    allparts[i+N0+N1] = pfne[i];
  }
  for(int i = 0; i < N3; i++) { 
    #pragma HLS UNROLL
    allparts[i+N0+N1+N2] = pfmu[i];
  }
}
template<unsigned int N0,unsigned int N1, unsigned int N2, unsigned int N3,unsigned int NOUT>
void tausort_in(PFChargedObj pfch[N0], PFChargedObj pfpho[N1], PFChargedObj pfne[N2], PFChargedObj pfmu[N3],
		hls::stream<PFChargedObj > axis_in[DATA_SIZE]) {
  #pragma HLS PIPELINE II=1
  for(int i = 0; i < N0; i++) {
    #pragma HLS UNROLL
    axis_in[i].write(pfch[i]);
  }
  for(int i = 0; i < N1; i++) {
    #pragma HLS UNROLL
    axis_in[i+N0].write(pfpho[i]);
  }
  for(int i = 0; i < N2; i++) {
    #pragma HLS UNROLL
    axis_in[i+N0+N1].write(pfne[i]);
  }
  for(int i = 0; i < N3; i++) {
    #pragma HLS UNROLL
    axis_in[i+N0+N1+N2].write(pfmu[i]);
  }
}
template<unsigned int N,unsigned int NR,unsigned int NMAX>
void deltaR_in(int iOffSet, etaphi_t seedeta,etaphi_t seedphi,PFChargedObj pfch[],hls::stream<axi_t > axis_in[]) {
  #pragma HLS inline
  #pragma HLS PIPELINE
  const ap_int<16> eDR2MAX = DR2MAX;
  PFChargedObj dummyc; dummyc.hwPt = 0; dummyc.hwEta = 0; dummyc.hwPhi = 0; dummyc.hwId = 0; dummyc.hwZ0 = 0;
  for (int i = 0; i < NMAX; i++) {
    #pragma HLS UNROLL
    int drcheck = dr2_int_cap<12>(seedeta,seedphi,pfch[i].hwEta,pfch[i].hwPhi,eDR2MAX);
    if(drcheck < DRCONE) {
      dummyc = pfch[i];
    } else {
      dummyc.hwPt = 0; dummyc.hwEta = 0; dummyc.hwPhi = 0; dummyc.hwId = 0; dummyc.hwZ0 = 0;
    }
    axi_t dummy = (dummyc.hwId,dummyc.hwPhi,dummyc.hwEta,dummyc.hwPt); 
    //std::cout << "writting ==> " << iOffSet+i << std::endl;
    axis_in[iOffSet+i].write(dummy);
  }
}
void algo_inputs_layer2_real_v3(hls::stream<axi_t> link_in[3*NTRACK+3*NCALO],hls::stream<axi_t> allparts_in [NTRACK*4+NCALO*4]) { 
  #pragma HLS PIPELINE 
  #pragma HLS stream variable=link_in     depth=36
  #pragma HLS stream variable=allparts_in depth=6

  PFChargedObj pfch1[DEPTH][NTRACK];
  PFChargedObj pfne1[DEPTH][NCALO];
  PFChargedObj pfch2[DEPTH][NTRACK];
  PFChargedObj pfne2[DEPTH][NCALO];
  PFChargedObj pfch3[DEPTH][NTRACK];
  PFChargedObj pfne3[DEPTH][NCALO];
  //#pragma HLS ARRAY_RESHAPE variable=pfch  complete
  //#pragma HLS ARRAY_RESHAPE variable=pfne  complete
  //Below has to be the same as in DATA note that #define doesn't seem to work
  #pragma HLS ARRAY_RESHAPE variable=pfch1   block factor=10 dim=2
  #pragma HLS ARRAY_RESHAPE variable=pfne1   block factor=6  dim=2
  //#pragma HLS RESOURCE      variable=pfch1   core=RAM_2P_BRAM
  //#pragma HLS RESOURCE      variable=pfne1   core=RAM_2P_BRAM
  #pragma HLS DEPENDENCE    variable=pfch1   intra false
  #pragma HLS DEPENDENCE    variable=pfne1   intra false

  #pragma HLS ARRAY_RESHAPE variable=pfch2   block factor=10 dim=2
  #pragma HLS ARRAY_RESHAPE variable=pfne2   block factor=6  dim=2
  //#pragma HLS RESOURCE      variable=pfch2   core=RAM_2P_BRAM
  //#pragma HLS RESOURCE      variable=pfne2   core=RAM_2P_BRAM
  #pragma HLS DEPENDENCE    variable=pfch2   intra false
  #pragma HLS DEPENDENCE    variable=pfne2   intra false

  #pragma HLS ARRAY_RESHAPE variable=pfch3   block factor=10 dim=2
  #pragma HLS ARRAY_RESHAPE variable=pfne3   block factor=6  dim=2
  //#pragma HLS RESOURCE      variable=pfch3   core=RAM_2P_BRAM
  //#pragma HLS RESOURCE      variable=pfne3   core=RAM_2P_BRAM
  #pragma HLS DEPENDENCE    variable=pfch3   intra false
  #pragma HLS DEPENDENCE    variable=pfne3   intra false

  //First take the seed from each region
  PFChargedObj chseed[NREGIONS];
  #pragma HLS ARRAY_PARTITION variable=chseed complete
  //#pragma HLS DEPENDENCE variable=chseed intra false
  LoopA:
  for(unsigned int idepth =0; idepth < DEPTH; idepth++) { 
    #pragma HLS PIPELINE II=1
    axi_t chin1[NTRACK];
    axi_t nein1[NCALO];
    axi_t chin2[NTRACK];
    axi_t nein2[NCALO];
    axi_t chin3[NTRACK];
    axi_t nein3[NCALO];
    for(int i0 = 0; i0 < NTRACK;  i0++) chin1[i0] = link_in[i0].read();
    for(int i0 = 0; i0 < NTRACK;  i0++) chin2[i0] = link_in[i0+NTRACK].read();
    for(int i0 = 0; i0 < NTRACK;  i0++) chin3[i0] = link_in[i0+2*NTRACK].read();
    for(int i0 = 0; i0 < NCALO;   i0++) nein1[i0] = link_in[3*NTRACK+i0].read();
    for(int i0 = 0; i0 < NCALO;   i0++) nein2[i0] = link_in[3*NTRACK+i0+NCALO].read();
    for(int i0 = 0; i0 < NCALO;   i0++) nein3[i0] = link_in[3*NTRACK+i0+2*NCALO].read();
    mp7_unpack_seed    (chin1[0],chseed[idepth]);
    mp7_unpack_seed    (chin2[0],chseed[DEPTH+idepth]);
    mp7_unpack_seed    (chin3[0],chseed[2*DEPTH+idepth]);
    mp7_unpack<NTRACK,0>       (chin1,pfch1[idepth]);
    mp7_unpack<NTRACK,0>       (chin2,pfch2[idepth]);
    mp7_unpack<NTRACK,0>       (chin3,pfch3[idepth]);
    mp7_unpack<NCALO,0>        (nein1,pfne1[idepth]);
    mp7_unpack<NCALO,0>        (nein2,pfne2[idepth]);
    mp7_unpack<NCALO,0>        (nein3,pfne3[idepth]);
  }
  //Next iterate through seeds and compute DR
  PFChargedObj seeds[NTAU];
  #pragma HLS ARRAY_PARTITION variable=seeds complete
  ptsort_hwopt_ind<PFChargedObj,NREGIONS,NTAU>(chseed, seeds);
  const ap_int<16> eDR2MAX = DR2MAX;
  ap_uint<6> seedsort[NTAU];
  #pragma HLS ARRAY_PARTITION variable=seedsort complete
  LoopB:
  for(int iseed0=0; iseed0 < NTAU; iseed0++) {
   #pragma HLS UNROLL
    seedsort[iseed0]=iseed0;
  }
  LoopC:
  for(int iseed0=0; iseed0 < NTAU; iseed0++) {
    #pragma HLS UNROLL
    for(int iseed1 = 1; iseed1 < NTAU; iseed1++) {
      if(iseed1 < iseed0+1) continue;                                                                                                                                                                                                                                       
      int drcheck = dr2_int_cap<16>(seeds[iseed0].hwEta,seeds[iseed0].hwPhi,seeds[iseed1].hwEta,seeds[iseed0].hwPhi,eDR2MAX);
      if(drcheck < DRCONE && iseed1 > iseed0) {
        seedsort[iseed1] = 0;
      }
    }
  }
  //#pragma HLS DEPENDENCE variable=allparts_in   intra false
  PFChargedObj dummyc; dummyc.hwPt = 0; dummyc.hwEta = 0; dummyc.hwPhi = 0; dummyc.hwId = 0; dummyc.hwZ0 = 0;
  LoopD:
  for(int itau = 0; itau < NTAU; itau++) {
    #pragma HLS UNROLL
    etaphi_t seedeta = seeds[itau].hwEta;
    etaphi_t seedphi = seeds[itau].hwPhi;
    ap_int<8> arr[4];
    #pragma HLS ARRAY_PARTITION variable=arr  complete 
    arraymap0(seedeta,seedphi,&(arr[0]));
    arraymap1(seedeta,seedphi,&(arr[1]));
    arraymap2(seedeta,seedphi,&(arr[2]));
    arraymap3(seedeta,seedphi,&(arr[3]));
    for(int iregion = 0; iregion < 4; iregion++) {
      #pragma HLS UNROLL
      int pRegion = arr[iregion];
      int NTRACKOFFSET   = iregion*NTRACK;//NTRACK;
      int NSELCALOOFFSET = iregion*NCALO+4*NTRACK;//NSELCALO;
      if(pRegion < 0) continue;
      if(pRegion < DEPTH) { 
       deltaR_in<NTRACK,NREGIONS,NTRACK>(NTRACKOFFSET,  seedeta,seedphi,pfch1[pRegion], allparts_in);
       deltaR_in<NCALO,NREGIONS,NCALO>  (NSELCALOOFFSET,seedeta,seedphi,pfch1[pRegion], allparts_in);
      } else if (pRegion < 2*DEPTH)  { 
       pRegion = pRegion-DEPTH;
       deltaR_in<NTRACK,NREGIONS,NTRACK>(NTRACKOFFSET,  seedeta,seedphi,pfch2[pRegion], allparts_in);
       deltaR_in<NCALO,NREGIONS,NCALO>  (NSELCALOOFFSET,seedeta,seedphi,pfne2[pRegion], allparts_in);
      } else { 
       pRegion = pRegion-2*DEPTH;
       deltaR_in<NTRACK,NREGIONS,NTRACK>(NTRACKOFFSET,  seedeta,seedphi,pfch3[pRegion], allparts_in);
       deltaR_in<NCALO,NREGIONS,NCALO>  (NSELCALOOFFSET,seedeta,seedphi,pfne3[pRegion], allparts_in);
      }
    }
  }
}
