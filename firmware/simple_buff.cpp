#include "simple_buff.h"

void buffer_input(MP7DataWord seed[MP7_NCHANN], MP7DataWord seed_out[MP7_NCHANN]) {
  hls::stream<MP7DataWord> seed_sub_tmp[MP7_NCHANN];
  //MP7DataWord seed_sub_tmp[MP7_NCHANN];
  #pragma HLS DATA_PACK variable=seed_sub_tmp
  #pragma HLS ARRAY_PARTITION variable=seed_sub_tmp complete
  #pragma HLS STREAM variable=seed_sub_tmp depth=10

  for (int icalo = 0; icalo < MP7_NCHANN; ++icalo) {
    //#pragma HLS latency min=1
    #pragma HLS LOOP UNROLL
    //seed_sub_tmp[icalo] = seed[icalo];
    seed_sub_tmp[icalo].write(seed[icalo]);
  }
  for (int icalo = 0; icalo < MP7_NCHANN; ++icalo) {
    //#pragma HLS latency min=1
    #pragma HLS LOOP UNROLL
    //seed_out[icalo] = seed_sub_tmp[icalo];
    seed_out[icalo] = seed_sub_tmp[icalo].read();
  }
}

void simple_buff(MP7DataWord input[MP7_NCHANN],MP7DataWord output[MP7_NCHANN]) { 
  #pragma HLS ARRAY_PARTITION variable=input complete
  #pragma HLS ARRAY_PARTITION variable=output complete
  #pragma HLS INTERFACE ap_none port=output
  #pragma HLS PIPELINE 

  //MP7DataWord input_buf[MP7_NCHANN];
  //#pragma HLS ARRAY_PARTITION variable=input_buf complete
  //buffer_input(input,input_buf);
  
  buffer_input(input,output);

}
