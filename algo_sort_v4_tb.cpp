#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "ap_int.h"
#include "firmware/algo_sort_v4.h"
#include "utils/pattern_serializer.h"
#include "utils/test_utils.h"

using namespace std;

int main(int argc, char ** argv) {
  float dphi  = 0.02; 
  float deta  = 0.02; 
  float pt    = 1.;
  float pterr = 0.1;
  PFChargedObj parts[DEPTH][DATA_SIZE];  
  for(int idepth = 0; idepth < DEPTH; idepth++) { 
    int   pEta  = idepth % 4;
    int   pPhi  = idepth/3;
    for(int i0 = 0; i0 < DATA_SIZE; i0++) { 
      int shift   = (i0/20);
      float pTEta = 0. + pEta*0.7+dphi*(shift+int(i0/2));
      float pTPhi = 0. + pPhi*0.7+deta*(2*shift+1+int(i0/2));
      PFChargedObj pTmp;
      pTmp.hwPt  = pt*PT_SCALE;
      pTmp.hwEta = pTEta*ETAPHI_SCALE;
      pTmp.hwPhi = pTPhi*ETAPHI_SCALE;
      pTmp.hwId  = shift;
      pTmp.hwZ0  = 0;
      parts[idepth][i0] = pTmp;
    }
  } 
  MP7PatternSerializer serInPatterns  ("mp7_input_patterns.txt",1);
  MP7PatternSerializer serOutPatterns ("mp7_output_patterns.txt",1);
  HumanReadablePatternSerializer serHR("human_readable_patterns.txt");
  MP7DataWord input[DEPTH][MP7_NCHANN];
  MP7DataWord output[DEPTH][MP7_NCHANN];
  for(int idepth = 0; idepth < DEPTH; idepth++) { 
    for(int i0 = 0; i0 < MP7_NCHANN; i0++) { 
      if(i0 > DATA_SIZE-1) break; //safety check
      input[idepth][i0] = ( parts[idepth][i0].hwId,  parts[idepth][i0].hwPhi , parts[idepth][i0].hwEta, parts[idepth][i0].hwPt  );
    } 
    algo_sort_v4(input[idepth],output[idepth]);
  } 
  //tmpaxi_t tau[DEPTH][NPART];  
  PFChargedObj parts_out[DEPTH][NTAU];  
  for(int idepth = 0; idepth < 1; idepth++) {
    for(int i0 = 0; i0 < MP7_NCHANN; i0++) { 
      PFChargedObj pTmp;
      pTmp.hwPt       = output[idepth][2*i0+0](15, 0);
      pTmp.hwEta      = output[idepth][2*i0+0](31,16);
      pTmp.hwPhi      = output[idepth][2*i0+1](15, 0);
      pTmp.hwId       = output[idepth][2*i0+1](31,16);
      float pPt  = pTmp.hwPt;pPt/=PT_SCALE;
      float pEta = pTmp.hwEta;pEta/=ETAPHI_SCALE;
      float pPhi = pTmp.hwPhi;pPhi/=ETAPHI_SCALE;
      if(i0 < NTAU) parts_out[idepth][i0] = pTmp;
      std::cout << "===> depth " << idepth << " -- part " << i0 << " vector " << pPt << "-- " << pEta << " -- " << pPhi << std::endl;
    }
  }
  for(int idepth = 0; idepth < DEPTH; idepth++) { 
    serInPatterns(input[idepth]); serOutPatterns(output[idepth]);
    serHR(parts[idepth],parts_out[idepth]);
  }
}

