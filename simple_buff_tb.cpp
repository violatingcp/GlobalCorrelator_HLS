#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "ap_int.h"
#include "firmware/simple_buff.h"

using namespace std;

typedef ap_int<16> pt_t;
typedef ap_int<10>  etaphi_t;
typedef ap_int<5>  vtx_t;
typedef ap_uint<3>  particleid_t;
typedef ap_int<10> z0_t;  // 40cm / 0.1
struct PFChargedObj {
  pt_t hwPt;
  etaphi_t hwEta, hwPhi; // relative to the region center, at calo                                                                                                                                                                     
  particleid_t hwId;
  z0_t hwZ0;
};

#define DATA_SIZE 36
#define DEPTH     128

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
      pt = rand() % 100;
      pTmp.hwPt  = pt;
      pTmp.hwEta = pTEta;
      pTmp.hwPhi = pTPhi;
      pTmp.hwId  = shift;
      pTmp.hwZ0  = 0;
      parts[idepth][i0] = pTmp;
    }
  } 
  for(int i0 = 0; i0 < 128; i0++) { 
    MP7DataWord  input    [DEPTH][MP7_NCHANN];
    MP7DataWord  output   [DEPTH][MP7_NCHANN];
    for(int idepth = 0; idepth < DEPTH; idepth++) { 
      for(int i0 = 0; i0 < MP7_NCHANN; i0++) { 
	if(i0 > DATA_SIZE-1) break; //safety check
	input[idepth][i0] = ( parts[idepth][i0].hwId,  parts[idepth][i0].hwPhi , parts[idepth][i0].hwEta, parts[idepth][i0].hwPt  );
      } 
      simple_buff(input[idepth],output[idepth]);
      //tmpaxi_t tau[DEPTH][NPART];  
      for(int idepth = 0; idepth < 1; idepth++) {
	for(int i0 = 0; i0 < MP7_NCHANN; i0++) { 
	  PFChargedObj pTmp;
	  pTmp.hwPt       = output[idepth][2*i0+0](15, 0);
	  pTmp.hwEta      = output[idepth][2*i0+0](31,16);
	  pTmp.hwPhi      = output[idepth][2*i0+1](15, 0);
	  pTmp.hwId       = output[idepth][2*i0+1](31,16);
	  float pPt  = pTmp.hwPt;//pPt/=PT_SCALE;
	  float pEta = pTmp.hwEta;//pEta/=ETAPHI_SCALE;
	  float pPhi = pTmp.hwPhi;//pPhi/=ETAPHI_SCALE;
	  //if(i0 < NTAU) parts_out[idepth][i0] = pTmp;
	  //std::cout << "===> depth " << idepth << " -- part " << i0 << " vector " << pPt << "-- " << pEta << " -- " << pPhi << std::endl;
	}
      }
    }
  }
}

