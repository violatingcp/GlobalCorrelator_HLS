#ifndef SIMPLE_BUFF_H
#define SIMPLE_BUFF_H

#include <ap_int.h>
#include <hls_stream.h>

typedef ap_uint<64> MP7DataWord;
#define MP7_NCHANN 72

void simple_buff(MP7DataWord input[MP7_NCHANN],MP7DataWord output[MP7_NCHANN]);

#endif
