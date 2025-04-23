#/*
#Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
#SPDX-License-Identifier: X11
#*/

//#include "cmdlineparser.h"
#include <iostream>
#include <cstring>

#include "ap_int.h"
#include <bitset>

// XRT includes
#include "xrt/xrt_bo.h"
#include <experimental/xrt_xclbin.h>
#include "xrt/xrt_device.h"
#include "xrt/xrt_kernel.h"

#define DATA_SIZE 54

typedef unsigned long long int ulli;

int main(int argc, char** argv)
{

    std::cout << "argc = " << argc << std::endl;
    for(int i=0; i < argc; i++)
    {
        std::cout << "argv[" << i << "] = " << argv[i] << std::endl;
    }

    // Read settings
    std::string binaryFile = "./puppi.xclbin";
    int device_index = 0;

    std::cout << "Open the device" << device_index << std::endl;
    auto device = xrt::device(device_index);
    std::cout << "Load the xclbin " << binaryFile << std::endl;
    auto uuid = device.load_xclbin("./puppi.xclbin");

    size_t vector_size_bytes = sizeof(ulli) * DATA_SIZE;

    //auto krnl = xrt::kernel(device, uuid, "vadd");
    auto krnl = xrt::kernel(device, uuid, "ReadWrite", xrt::kernel::cu_access_mode::exclusive);

    std::cout << "Allocate Buffer in Global Memory\n";
    auto boIn1 = xrt::bo(device, vector_size_bytes, krnl.group_id(0)); //Match kernel arguments to RTL kernel
    auto boIn2 = xrt::bo(device, vector_size_bytes, krnl.group_id(1));

    // Map the contents of the buffer object into host memory
    auto bo0_map = boIn1.map<ulli*>();
    auto bo1_map = boIn2.map<ulli*>();
    std::fill(bo0_map, bo0_map + DATA_SIZE, 0);
    std::fill(bo1_map, bo1_map + DATA_SIZE, 0);

    const ulli input[] = {0x14a01f,0x0d001d,0x14401c,0x124006,0x16e006,0x172006,0x66005,0x086005,0x0,0x1ce022,0x19201e,0x2e6009,0x18e006,0x1ae006,0x2ee006,0x270006,0x290006,0x0,0x30602a,0x3ce025,0x46600a,0x3ae007,0x472007,0x3e8006,0x408006,0x36a006,0x0,0x548035,0x486030,0x590020,0x49201d,0x528009,0x52a008,0x54a008,0x4a8007,0x0,0x748031,0x650025,0x728009,0x6e6008,0x706008,0x766008,0x628008,0x648008,0x0,0x0886029,0x084e028,0x0812022,0x786008,0x082e008,0x0866007,0x08e6006,0x7ec006,0x0};

    // Create the test data
    for (int i = 0; i < DATA_SIZE; ++i)
    {
        bo0_map[i] = input[i];
    }

    for (int i = 0; i < DATA_SIZE; ++i)
    {
        std::cout << bo0_map[i] << std::endl;
    }
    std::cout << "\noutput: \n" ;
    // Synchronize buffer content with device side
    std::cout << "synchronize input buffer data to device global memory\n";
    boIn1.sync(XCL_BO_SYNC_BO_TO_DEVICE);

    std::cout << "Execution of the kernel\n";
    auto run = krnl(boIn1,boIn2); //DATA_SIZE=size
    run.wait();

    // Get the output;
    std::cout << "Get the output data from the device" << std::endl;
    boIn2.sync(XCL_BO_SYNC_BO_FROM_DEVICE);

    // Validate results
    for (int i = 0; i < DATA_SIZE; ++i)
    {
        std::cout << bo1_map[i] << std::endl;
        std::cout << std::bitset<64>(bo1_map[i]) << std::endl;
    }

    return 0;
}


