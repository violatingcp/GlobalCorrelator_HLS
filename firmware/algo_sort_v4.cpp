
#include "algo_sort_v4.h"
//#include "sorting_network_axi.h"
//#include "bitonic_sort.h"
#include "tau_nn.h"

//16+10+10+3+10 for pt,eta,phi,particleid,z0
//typedef ap_axis <64*NPART,1,1,1> axi_t;
//typedef hls::stream<axi_t> hls::stream<axi_t>;
//stream depth is TM6 and 240 MHz gives 24 clocks

template<typename T, int NIn, int NOut>
void ptsort_hwopt_ind(T in[NIn], T out[NOut]) { 
    #pragma HLS PIPELINE
    T tmp[NOut];
    #pragma HLS ARRAY_PARTITION variable=tmp complete
    for (int iout = 0; iout < NOut; ++iout) {
        #pragma HLS unroll
        tmp[iout] = 0;
    }

    for (int it = 0; it < NIn; ++it) {
        for (int iout = NOut-1; iout >= 0; --iout) {
            if (tmp[iout] <= in[it]) {
                if (iout == 0 || tmp[iout-1] > in[it]) {
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
template<unsigned int N,unsigned int OFFSET, unsigned int OFFSET2> 
inline void mp7_unpack(MP7DataWord data[MP7_NCHANN], PFChargedObj data_in[DATA_SIZE]) {
    #pragma HLS inline
    #pragma HLS PIPELINE
    for (unsigned int i = 0; i < N; ++i) {
      PFChargedObj pTmp;
      //32bit below
      //pTmp.hwPt       = data[OFFSET+2*i+0](15, 0);
      //pTmp.hwEta      = data[OFFSET+2*i+0](31,16);
      //pTmp.hwPhi      = data[OFFSET+2*i+1](15, 0);
      //pTmp.hwId       = data[OFFSET+2*i+1](31,16);
      //64bit below
      pTmp.hwPt       = data[OFFSET+i](15, 0);
      pTmp.hwEta      = data[OFFSET+i](31,16);
      pTmp.hwPhi      = data[OFFSET+i](47,32);
      pTmp.hwId       = data[OFFSET+i](63,48);
      data_in[i+OFFSET2] = pTmp;
    }
}
template<unsigned int N>
void make_inputs(input_t nn_data[N*8], axi_t pf[DATA_SIZE]) {
  //#pragma HLS inline
  #pragma HLS PIPELINE
  for (int i = 0; i < N; i++) {
    #pragma HLS PIPELINE II=1
    axi_t tmpobj = pf[i];
    nn_data[i*8+0] = input_t(tmpobj(15, 0));
    nn_data[i*8+1] = input_t(tmpobj(31,16));
    nn_data[i*8+2] = input_t(tmpobj(47,32));
    ap_uint<16>  id = tmpobj(63, 48);
    nn_data[i*8+3] = input_t(id == 2 ? 1 : 0);
    nn_data[i*8+4] = input_t(id == 3 ? 1 : 0);
    nn_data[i*8+5] = input_t(id == 4 ? 1 : 0);
    nn_data[i*8+6] = input_t(id == 1 ? 1 : 0);
    nn_data[i*8+7] = input_t(id == 0 ? 1 : 0);
  }
} 
void algo_sort_v4(MP7DataWord input[MP7_NCHANN],MP7DataWord output[MP7_NCHANN]) { 
//void algo_sort_v4(hls::stream<axi_t> input[MP7_NCHANN],hls::stream<axi_t> output[MP7_NCHANN]) { 
		  //MP7DataWord input[MP7_NCHANN],MP7DataWord output[MP7_NCHANN]) { 
  #pragma HLS ARRAY_PARTITION variable=input complete
  #pragma HLS ARRAY_PARTITION variable=output complete
  //#pragma HLS stream variable=input     depth=12
  #pragma HLS INTERFACE ap_none port=output
  //#pragma HLS DATAFLOW 
  #pragma HLS PIPELINE II=12
  for(int i1 = 0; i1 < NTAU; i1++) { 
   axi_t data[64];
   #pragma HLS ARRAY_PARTITION variable=data complete
   axi_t data_out[64];
   #pragma HLS ARRAY_PARTITION variable=data_out complete
   for(int i0 = 0; i0 < 64; i0++) data[i0] = input[i0];//.read(input[i0]);
   //sorting_network_64_in(data);
   //BitonicSortOptimized(data);
   ptsort_hwopt_ind<axi_t,64,10>(data,data_out);
   input_t  nn_data[NTAUPARTS*8];
   make_inputs<NTAUPARTS>(nn_data,data_out);
   result_t taunn[N_OUTPUTS];
   axi_t dummyc = data_out[0];
   tau_nn(nn_data,taunn);
   dummyc(15,0) = taunn[0]*100;
   output[i1] = dummyc;
   //for(int i0 = 1; i0 < 10; i0++) output[i0] = data_out[i0];//.read(data[i0]);
  }
}

