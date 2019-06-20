#ifndef ALGO_INPUTS_LAYER2_REAL_V2_H
#define ALGO_INPUTS_LAYER2_REAL_V2_H

#define AP_INT_MAX_W 1600
#include "data.h"

#define DATA_SIZE 64
#define NTAU  6
#define NREGIONS 36
#define NPART 25
#define DEPTH NREGIONS
#define NTAUPARTS  10
#define DRCONE 8410
#define DR2MAX 10000
#define MP7_NCHANN 72

static float PT_SCALE = 4.0;     // quantize in units of 0.25 GeV (can be changed)
static float ETAPHI_FACTOR = 4;  // size of an ecal crystal in phi in integer units (our choice)
static float ETAPHI_SCALE = ETAPHI_FACTOR*(180./M_PI);  // M_PI/180 is the size of an ECal crystal; we make a grid that is 4 times that size
static int16_t PHI_WRAP = 360*ETAPHI_FACTOR;            // what is 3.14 in integer
//typedef ap_axis <64,1,1,1> axi_t;
void algo_inputs_layer2_real_v2(hls::stream<axi_t> link_in[NTRACK+NCALO],hls::stream<axi_t> allparts_in [4*NCALO+4*NTRACK]);

#endif
