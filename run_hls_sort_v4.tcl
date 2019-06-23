open_project -reset proj_sort_v4
set_top algo_sort_v4
add_files firmware/tau_nn.cpp
add_files firmware/algo_sort_v4.cpp 
#add_files firmware/bitonic_sort.cpp
add_files -tb utils/pattern_serializer.cpp -cflags "-DTESTMP7"
add_files -tb algo_sort_v4_tb.cpp 
open_solution -reset "solution1"

#Specify FPGA and clock constraints
set_part {xcku115-flvd1517-2-i}
create_clock -period 3.0 -name default
#create_clock -period 4.166667 -name default
#create_clock -period 6.25 -name default
set_clock_uncertainty 1.5
config_interface -trim_dangling_port
#catch {config_array_partition -maximum_size 4096}

csim_design 
csynth_design
#cosim_design  -trace_level all
export_design -flow syn -format ip_catalog -vendor "cern-cms" -version 1.0 -description algo_sort_v4

