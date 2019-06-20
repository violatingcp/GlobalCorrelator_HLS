//
//    rfnoc-hls-neuralnet: Vivado HLS code for neural-net building blocks
//
//    Copyright (C) 2017 EJ Kreinar
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef NNET_LAYER_H_
#define NNET_LAYER_H_

#include "nnet_common.h"
#include "hls_stream.h"
#include <math.h>

namespace nnet {

struct layer_config
{
    // Internal data type definitions
    typedef float bias_t;
    typedef float weight_t;
    typedef float accum_t;

    // Layer Sizes
    static const unsigned n_in = 10;
    static const unsigned n_out = 10;

    // Resource reuse info
    static const unsigned io_type = io_parallel;
    static const unsigned reuse_factor = 1;
    static const bool store_weights_in_bram = false;
    static const unsigned n_zeros = 0;
    // partitioning arrays cyclically to go with roll factors?
};


#define DIV_ROUNDUP(n,d) ((n + d - 1) / d)
#define ADD_LAT 5

 template<class data_T, class res_T, typename CONFIG_T>
void compute_layer(
    data_T    data[CONFIG_T::n_in],
    res_T     res[CONFIG_T::n_out],
    typename CONFIG_T::weight_t  weights[CONFIG_T::n_in*CONFIG_T::n_out],
    typename CONFIG_T::bias_t    biases[CONFIG_T::n_out])
{
    data_T cache;
    typename CONFIG_T::accum_t mult[CONFIG_T::n_in*CONFIG_T::n_out];
    typename CONFIG_T::accum_t acc[CONFIG_T::n_out];

    // Use a function_instantiate in case it helps to explicitly optimize unchanging weights/biases
    //#pragma HLS function_instantiate variable=weights,biases

    if (CONFIG_T::io_type == io_parallel){
        // For parallel inputs:
        //   - completely partition arrays -- target fabric
        //   - if we have an unroll factor, limit number of multipliers
        #pragma HLS PIPELINE II=CONFIG_T::reuse_factor

        static const int block_factor       = (CONFIG_T::n_in*CONFIG_T::n_out/CONFIG_T::reuse_factor); //DIV_ROUNDUP(CONFIG_T::n_in*CONFIG_T::n_out, CONFIG_T::reuse_factor);
	std::cout << " ---> " << CONFIG_T::n_in*CONFIG_T::n_out << " -- " << CONFIG_T::reuse_factor << " --- " << (CONFIG_T::n_in*CONFIG_T::n_out/CONFIG_T::reuse_factor) << std::endl;
        //#pragma HLS ARRAY_PARTITION variable=weights block factor=block_factor // remove this line for now, it breaks compression sometimes
        //#pragma HLS ARRAY_PARTITION variable=weights complete // remove this line for now, it breaks compression sometimes
        #pragma HLS ARRAY_PARTITION variable=biases complete
        //#pragma HLS ARRAY_RESHAPE   variable=weights block factor=block_factor
        #pragma HLS ARRAY_RESHAPE   variable=mult block factor=block_factor
        #pragma HLS ARRAY_PARTITION   variable=mult complete
        #pragma HLS ARRAY_PARTITION variable=acc complete

	int multiplier_limit  = ceil(float(CONFIG_T::n_in*CONFIG_T::n_out) / float(CONFIG_T::reuse_factor)) - floor(float(CONFIG_T::n_zeros) / float(CONFIG_T::reuse_factor));
	#pragma HLS ALLOCATION instances=mul limit=multiplier_limit operation

    } else if (CONFIG_T::io_type == io_serial){
        #pragma HLS ARRAY_RESHAPE variable=weights complete dim=2
        #pragma HLS ARRAY_PARTITION variable=mult complete dim=2
        #pragma HLS ARRAY_PARTITION variable=acc complete dim=1
        #pragma HLS DATAFLOW
        #pragma HLS STREAM variable=mult depth=1
        #pragma HLS STREAM variable=acc depth=1
    }

    // Do the matrix-multiply
    Product1: for(int ii = 0; ii < CONFIG_T::n_in; ii++) {
        if (CONFIG_T::io_type == io_serial){
            #pragma HLS PIPELINE
        }
        cache = data[ii];
        Product2: for(int jj = 0; jj < CONFIG_T::n_out; jj++) {
            if (CONFIG_T::io_type == io_serial) {
                int multiplier_limit  = ceil(float(CONFIG_T::n_out) / float(CONFIG_T::reuse_factor));
                #pragma HLS ALLOCATION instances=mul limit=multiplier_limit operation
            }
            mult[ii*CONFIG_T::n_out+jj] = cache * weights[ii*CONFIG_T::n_out+jj];
        }
    }

    // Initialize accumulator with input biases
    ResetAccum: for(int iacc = 0; iacc < CONFIG_T::n_out; iacc++) {
        if (CONFIG_T::io_type == io_serial){
            #pragma HLS UNROLL
        }
        acc[iacc] = (typename CONFIG_T::accum_t) biases[iacc];
    }

    // special loop for accumulation
    typename CONFIG_T::accum_t acc_lat[CONFIG_T::n_out][ADD_LAT];
    #pragma HLS ARRAY_PARTITION variable=acc_lat complete dim=0
    AddLatencyInit: 
    for (int ii = 0; ii < CONFIG_T::n_out; ii++){
      for (int ij= 0; ij < ADD_LAT; ij++){
        #pragma UNROLL
	acc_lat[ii][ij] = 0;
      }
    }
    for(int ii = 0; ii < CONFIG_T::n_in; ii++) {
	for (int io = 0; io < CONFIG_T::n_out; io++){
         #pragma HLS UNROLL
         for (int ia = 0; ia < ADD_LAT; ia++){
          #pragma HLS UNROLL
	  int index = ii*CONFIG_T::n_out+io;
	  acc_lat[io][ia] += mult[index];
	}
      }
    }

 FullAccum: 
    for (int ii = 0; ii < CONFIG_T::n_out; ii++){
            #pragma HLS UNROLL
      for (int ij= 0; ij < ADD_LAT; ij++){
	acc[ii] += acc_lat[ii][ij];
      }
     }
    // Cast to "res_t" type
    Result: for(int ires = 0; ires < CONFIG_T::n_out; ires++){
        if (CONFIG_T::io_type == io_serial){
            #pragma HLS UNROLL
        }
        res[ires] = (res_T) (acc[ires]);
    }    
}

}

#endif
