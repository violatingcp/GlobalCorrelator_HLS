#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "ap_int.h"
#include "firmware/algo_inputs_layer2_real_v2.h"
//#include "utils/pattern_serializer.h"
//#include "utils/test_utils.h"

using namespace std;

typedef ap_uint<64>           axi_t;
typedef hls::stream<axi_t> stream_t;

int main(int argc, char ** argv) {
  hls::stream<axi_t>     link_in[NTRACK+NCALO];
  hls::stream<axi_t>     link_out[4*NTRACK+4*NCALO];

  axi_t     arr_link_in[DEPTH][NTRACK+NCALO];
  axi_t     arr_link_out[NTAU][4*NTRACK+4*NCALO];

  axi_t chparts[DEPTH][NTRACK];  
  axi_t neparts[DEPTH][NCALO];  
  PFChargedObj chparts_out[DEPTH][NTRACK];  
  PFChargedObj neparts_out[DEPTH][NCALO];  
  PFChargedObj parts_out  [NTAU][4*NTRACK+4*NCALO];

  float dphi  = 0.02; 
  float deta  = 0.02; 
  float pt    = 1.;
  float pterr = 0.1;
  for(int idepth = 0; idepth < DEPTH; idepth++) { 
     for(int i0 = 0; i0 < NTRACK; i0++) { 
      int pShift = int(i0/NTRACK);
      int   pEta  = idepth + pShift % 4;
      int   pPhi  = (idepth+ pShift)/3;
      float pTEta = 0. + pEta*0.7+dphi*int(i0/2);
      float pTPhi = 0. + pPhi*0.7+deta*(1+int(i0/2));
      axi_t pTmp;
      pTmp.range(15, 0) = pt*PT_SCALE;
      pTmp.range(24,16) = pTEta*ETAPHI_SCALE;
      pTmp.range(34,25) = pTPhi*ETAPHI_SCALE;
      pTmp.range(38,35) = 0;
      pTmp.range(48,39) = 0;
      chparts[idepth+pShift][i0-NTRACK*pShift] = pTmp;

      PFChargedObj pPartTmp;
      pPartTmp.hwPt  = pt*PT_SCALE;
      pPartTmp.hwEta = pTEta*ETAPHI_SCALE;
      pPartTmp.hwPhi = pTPhi*ETAPHI_SCALE;
      pPartTmp.hwId  = 0;
      chparts_out[idepth+pShift][i0-NTRACK*pShift] = pPartTmp;
    }
    for(int i0 = 0; i0 < NCALO; i0++) { 
      int pShift = int(i0/NCALO);
      int   pEta  = idepth + pShift % 4;
      int   pPhi  = (idepth+ pShift)/3;
      float pTEta = 0. + pEta*0.7+dphi*(2+int(i0/2));
      float pTPhi = 0. + pPhi*0.7+deta*(6+int(i0/2));
      axi_t pTmp;
      pTmp.range(15, 0) = pt*PT_SCALE;
      pTmp.range(24,16) = pTEta*ETAPHI_SCALE;
      pTmp.range(34,25) = pTPhi*ETAPHI_SCALE;
      pTmp.range(38,35) = 2;
      pTmp.range(55,39) = pt*PT_SCALE;
      neparts[idepth+pShift][i0-NCALO*pShift] = pTmp;
    
      PFChargedObj pPartTmp;
      pPartTmp.hwPt  = pt*PT_SCALE;
      pPartTmp.hwEta = pTEta*ETAPHI_SCALE;
      pPartTmp.hwPhi = pTPhi*ETAPHI_SCALE;
      pPartTmp.hwId  = 1;
      neparts_out[idepth+pShift][i0-NCALO*pShift] = pPartTmp;
    }
  } 
  //MP7PatternSerializer serInPatterns  ("mp7_input_patterns.txt",1);
  //MP7PatternSerializer serOutPatterns ("mp7_output_patterns.txt",1);
  //HumanReadablePatternSerializer serHR("human_readable_patterns.txt");
  for(int idepth = 0; idepth < DEPTH; idepth++) { 
    for(int i0 = 0; i0 < NTRACK;  i0++) link_in[i0         ].write(chparts[idepth][i0]);  
    //for(int i0 = 0; i0 < NTRACK;  i0++) link_in[i0+NTRACK  ].write(chparts[idepth+1][i0]);  
    //for(int i0 = 0; i0 < NTRACK;  i0++) link_in[i0+2*NTRACK].write(chparts[idepth+2][i0]);  
    for(int i0 = 0; i0 < NCALO;   i0++) link_in[i0+NTRACK].write(neparts[idepth][i0]);  
    //for(int i0 = 0; i0 < NCALO;   i0++) link_in[i0+2*NTRACK+NCALO].write(neparts[idepth+1][i0]);  
    //for(int i0 = 0; i0 < NCALO;   i0++) link_in[i0+3*NTRACK+2*NCALO].write(neparts[idepth+2][i0]);  

    for(int i0 = 0; i0 < NTRACK;  i0++) arr_link_in[idepth][i0                       ] = chparts[idepth][i0];
    //for(int i0 = 0; i0 < NTRACK;  i0++) arr_link_in[idepth][i0              +  NTRACK] = chparts[idepth+1][i0];
    //for(int i0 = 0; i0 < NTRACK;  i0++) arr_link_in[idepth][i0              +2*NTRACK] = chparts[idepth+2][i0];
    for(int i0 = 0; i0 < NCALO;   i0++) arr_link_in[idepth][i0              +NTRACK] = neparts[idepth][i0];  
    //for(int i0 = 0; i0 < NCALO;   i0++) arr_link_in[idepth][i0+NCALO        +2*NTRACK] = neparts[idepth+1][i0];  
    //for(int i0 = 0; i0 < NCALO;   i0++) arr_link_in[idepth][i0+2*NCALO      +3*NTRACK] = neparts[idepth+2][i0];  
  }
  algo_inputs_layer2_real_v2(link_in,  link_out);
  for(int idepth = 0; idepth < NTAU; idepth++) {
    for(int ipart = 0; ipart < 4*(NTRACK+NCALO); ipart++) { 
      axi_t pTmpOut;
      PFChargedObj pTmp;
      link_out[ipart].read(pTmpOut);
      arr_link_out[idepth][ipart] = pTmpOut;
      float pPt  = pTmpOut.range(15,0);  pTmp.hwPt  = pPt;  pPt/=PT_SCALE;
      float pEta = pTmpOut.range(31,16); pTmp.hwEta = pEta; pEta/=ETAPHI_SCALE;
      float pPhi = pTmpOut.range(47,32); pTmp.hwPhi = pPhi; pPhi/=ETAPHI_SCALE;
      int   pId  = pTmpOut.range(63,48); pTmp.hwId  = pId;
      parts_out[idepth][ipart] = pTmp;
      std::cout << "===> tau part " << idepth << " -- part " << ipart << " vector " << pPt << "-- " << pEta << " -- " << pPhi << " - Id - " << pId << std::endl;
    }
  }
  /*
  for(int idepth = 0; idepth < DEPTH; idepth++) { 
    serInPatterns(arr_link_in[idepth]); serOutPatterns(arr_link_out[idepth]);
    serHR(chparts_out[idepth], emparts_out[idepth], neparts_out[idepth], muparts_out[idepth],parts_out[idepth]);
  } 
  */
}

