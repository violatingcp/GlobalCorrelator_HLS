open_project -reset proj_inputs_v3
set_top algo_inputs_layer2_v3
add_files firmware/algo_inputs_layer2_v3.cpp -cflags "-DTESTMP7"
add_files -tb algo_inputs_layer2_v3_tb.cpp 
add_files -tb utils/pattern_serializer.cpp -cflags "-DTESTMP7"
open_solution -reset "solution1"

#Specify FPGA and clock constraints
set_part {xcku115-flvd1517-2-i}
create_clock -period 4.166667 -name default
set_clock_uncertainty 1.5
config_interface -trim_dangling_port

csim_design 
csynth_design
#cosim_design  -trace_level all
export_design -format ip_catalog -vendor "cern-cms" -version 1.0 -description algo_inputs_layer2_v3
