#include "algo_topIP1.h"
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include<string>
#include <dirent.h>
#include <unistd.h>

#define N_TEST_FILES 1

int main() {
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != nullptr) {
        std::cout << "Current working directory: " << cwd << std::endl;
    } else {
        std::cerr << "Error: Unable to get current working directory" << std::endl;
        return 1;
    }

    loop start, end ;

    ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS] ;
    ap_uint<LINK_WIDTH> link_out[N_OUTPUT_LINKS] ;

    std::string patternFile_dir = "../../../../../data/inputPatterFiles/";

    DIR* dir = opendir(patternFile_dir.c_str());
    if (dir == nullptr) {
        std::cerr << "Error: Unable to open directory" << std::endl;
        return 1;
    }

    dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0 || strstr(entry->d_name,".links") ==0)
            continue;

        std::string file_path = patternFile_dir + entry->d_name;
        std::cout << "File path: " << file_path << std::endl;	
        
	std::ifstream read(file_path);
        if (!read.is_open()) {
            std::cout << "Error: Failed to open file " << file_path << std::endl;
            continue;
        }

        std::string line;
        std::getline(read, line);
        //std::cout << line << std::endl;
        std::getline(read, line);
        //std::cout << line << std::endl;

        std::istringstream word(line);
        std::string link;
        int index = 0;
        while (word >> link && index < N_INPUT_LINKS) {
        	//std::cout << index << std::endl;
        	//std::cout << "B: " << link << std::endl;

        	link_in[index] = ap_uint<LINK_WIDTH>(link.c_str(),2);

        	//std::cout << "\nA: ";
        	//for(int i = 0; i < LINK_WIDTH; i++){
        	//	std::cout << link_in[index][LINK_WIDTH-(i+1)];
        	//}
        	//std::cout << "\n";
        	index++;
        }

        algo_topIP1(link_in, link_out);

        hftower HFTowers;
	/*
	std::cout<<"Output Towers ! \n";
        for(loop i = 0; i < N_OUTPUT_LINKS; i++){
        	for(loop j = 0; j < TOWERS_IN_ETA-1; j++){
    			ap_uint<10> start =j*10;
    			ap_uint<10> end = start+9;
    			HFTowers.fillhftower(((ap_uint<10>) link_out[i].range(end, start))) ;
    			cout << "  > tower  [ "
    			 << " phi slice : " << i<<" | "<< " eta slice : "<<j<<" ] "
			 <<" energy  :" <<  HFTowers.energy<< endl ;
		}
          }
        */
    }

    closedir(dir);

    return 0;
}
