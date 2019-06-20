#ifndef __SORTNING_NETWORK_H__
#define __SORTNING_NETWORK_H__

//#include "GlobalCorrelator_HLS/firmware/simple_fullpfalgo.h"
#include "hls_stream.h"
#include "data.h"

#define NTAU 6
#define DATA_SIZE 64
#define NTAUPARTS 10

void swap1(axi_t &data1,axi_t &data2);
void swap2(axi_t &data1,axi_t &data2);
void sorting_network_64_in(axi_t datas[DATA_SIZE]);
void sorting_network_64(hls::stream<axi_t > data_in[], hls::stream<axi_t > data_out[]);
void sorting_network_64_nn(axi_t data_in[DATA_SIZE],axi_t data_out[NTAU]);
void sorting_network_128_in(axi_t datas[DATA_SIZE]);
void sorting_network_128(hls::stream<axi_t > data_in[], hls::stream<axi_t > data_out[]);
void sorting_network_128_nn(hls::stream<axi_t> data_in[DATA_SIZE],hls::stream<PFChargedObj> &data_out);
#endif
