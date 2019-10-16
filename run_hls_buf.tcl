open_project -reset simple_buff
set_top simple_buff
add_files firmware/simple_buff.cpp 
add_files -tb simple_buff_tb.cpp 
open_solution -reset "solution1"

#Specify FPGA and clock constraints
#set_part {xcku115-flvd1517-2-i}
set_part {xcvu9p-flga2104-2L-e}
create_clock -period 3.0 -name default
#create_clock -period 4.166667 -name default
#create_clock -period 6.25 -name default
#set_clock_uncertainty 1.5
config_interface -trim_dangling_port

csim_design 
csynth_design
#cosim_design  -trace_level all
#export_design -format ip_catalog -vendor "cern-cms" -version 1.0 -description algo_inputs_layer2_v3
export_design  -flow syn -format ip_catalog -vendor "cern-cms" -version 1.0 -description algo_inputs_layer2_v3

