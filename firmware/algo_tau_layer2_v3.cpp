#include "algo_tau_layer2_v3.h"
#include "tau_nn.h"
#include "sorting_network.h"

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
template<unsigned int N,unsigned int OFFSET, unsigned int OFFSET2> 
inline void mp7_unpack(MP7DataWord data[MP7_NCHANN], hls::stream<PFChargedObj> data_in[DATA_SIZE]) {
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
      data_in[i+OFFSET2].write(pTmp);
    }
}
template<unsigned int N>
void make_inputs(input_t nn_data[N*8], hls::stream<PFChargedObj> pf[N]) {
  //#pragma HLS inline
  #pragma HLS PIPELINE
  for (int i = 0; i < N; i++) {
    #pragma HLS PIPELINE II=1
    PFChargedObj tmpobj = pf[i].read();
    nn_data[i*8+0] = input_t(tmpobj.hwPt);
    nn_data[i*8+1] = input_t(tmpobj.hwEta);
    nn_data[i*8+2] = input_t(tmpobj.hwPhi);
    nn_data[i*8+3] = input_t(tmpobj.hwId == 2 ? 1 : 0);
    nn_data[i*8+4] = input_t(tmpobj.hwId == 3 ? 1 : 0);
    nn_data[i*8+5] = input_t(tmpobj.hwId == 4 ? 1 : 0);
    nn_data[i*8+6] = input_t(tmpobj.hwId == 1 ? 1 : 0);
    nn_data[i*8+7] = input_t(tmpobj.hwId == 0 ? 1 : 0);
  }
} 
void algo_tau_layer2_v3(MP7DataWord input[MP7_NCHANN],MP7DataWord output[MP7_NCHANN]) { 
  #pragma HLS ARRAY_PARTITION variable=input complete
  #pragma HLS ARRAY_PARTITION variable=output complete
  #pragma HLS INTERFACE ap_none port=output
  #pragma HLS PIPELINE II=12
  //#pragma HLS INTERFACE axis port=allparts_in
  //#pragma HLS INTERFACE axis port=allparts_out
  //#pragma HLS INTERFACE s_axilite port=return bundle=ctrl

  hls::stream<PFChargedObj > allparts_in[DATA_SIZE];
  hls::stream<PFChargedObj > allparts_out;
  #pragma HLS stream variable=allparts_in   depth=12
  #pragma HLS stream variable=allparts_out  depth=6
  //#pragma HLS stream variable=allparts_out depth=3  
  //#pragma HLS stream variable=link_out     depth=5

  //Doing 64 for now
  mp7_unpack<64,0,0>  (input,allparts_in);
  //mp7_unpack<64,0,64>  (input,allparts_in);
  //mp7_unpack<28,0,36> (input,allparts_in);
  //  mp7_unpack<36,0,72> (input,allparts_in);
  //  mp7_unpack<20,0,108>(input,allparts_in);
  sorting_network_64_nn(allparts_in,allparts_out);
  for(int i0 = 0; i0 < NTAU; i0++) { 
    PFChargedObj tmp_out; 
    allparts_out.read(tmp_out);
    output[i0] = ( tmp_out.hwId,  tmp_out.hwPhi , tmp_out.hwEta, tmp_out.hwPt );
  }
}

