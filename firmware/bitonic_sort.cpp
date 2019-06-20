#include "bitonic_sort.h"

template<class array_t, class index_t>
void exchange(array_t a[], index_t i, index_t j) {
    #pragma HLS INLINE
  array_t t = a[i];
  a[i] = a[j];
  a[j] = t;
}

void BitonicSortOptimized(axi_t work_array[64]) {
  // HLS csynth says it can't meeting timing requirements if II<4 for N=16
  // HLS csynth says it can't meeting timing requirements if II<5 for N=32
  // HLS csynth says it can't meeting timing requirements if II<7 for N=64
  // HLS csynth says it can't meeting timing requirements if II<10 for N=128
  #pragma HLS ARRAY_PARTITION variable=work_array complete
  int i,j;

 LOOP2: for (j=1; j>0; j=j>>1) {
  ILOOP2: for (i=0; i<64; i++) {
      int ij=i^j;
      if ((ij)>i) {
	if ((i&2)==0 && work_array[i] > work_array[ij]) 
	  exchange(work_array,i,ij);
	if ((i&2)!=0 && work_array[i] < work_array[ij])
	  exchange(work_array,i,ij);
      }
    }
  }

 LOOP4: for (j=2; j>0; j=j>>1) {
  ILOOP4: for (i=0; i<64; i++) {
      int ij=i^j;
      if ((ij)>i) {
	if ((i&4)==0 && work_array[i] > work_array[ij]) 
	  exchange(work_array,i,ij);
	if ((i&4)!=0 && work_array[i] < work_array[ij])
	  exchange(work_array,i,ij);
      }
    }
  }

 LOOP8: for (j=4; j>0; j=j>>1) {
  ILOOP8: for (i=0; i<64; i++) {
      int ij=i^j;
      if ((ij)>i) {
	if ((i&8)==0 && work_array[i] > work_array[ij]) 
	  exchange(work_array,i,ij);
	if ((i&8)!=0 && work_array[i] < work_array[ij])
	  exchange(work_array,i,ij);
      }
    }
  }

 LOOP16: for (j=8; j>0; j=j>>1) {
  ILOOP16: for (i=0; i<64; i++) {
      int ij=i^j;
      if ((ij)>i) {
	if ((i&16)==0 && work_array[i] > work_array[ij]) 
	  exchange(work_array,i,ij);
	if ((i&16)!=0 && work_array[i] < work_array[ij])
	  exchange(work_array,i,ij);
      }
    }
  }
 LOOP32: for (j=16; j>0; j=j>>1) {
  ILOOP32: for (i=0; i<64; i++) {
      int ij=i^j;
      if ((ij)>i) {
	if ((i&32)==0 && work_array[i] > work_array[ij]) 
	  exchange(work_array,i,ij);
	if ((i&32)!=0 && work_array[i] < work_array[ij])
	  exchange(work_array,i,ij);
      }
    }
  }
 LOOP64: for (j=32; j>0; j=j>>1) {
  ILOOP64: for (i=0; i<64; i++) {
      int ij=i^j;
      if ((ij)>i) {
	if ((i&64)==0 && work_array[i] > work_array[ij]) 
	  exchange(work_array,i,ij);
	if ((i&64)!=0 && work_array[i] < work_array[ij])
	  exchange(work_array,i,ij);
      }
    }
  }
}
