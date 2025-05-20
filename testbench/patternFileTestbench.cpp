#include "algo_topIP1.h"
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <dirent.h>
#include <unistd.h>
#include <bitset>

#define N_TEST_FILES 1

int main() {
    char cwd[1024];
    if (getcwd(cwd, sizeof(cwd)) != nullptr) {
        std::cout << "Current working directory: " << cwd << std::endl;
    } else {
        std::cerr << "Error: Unable to get current working directory" << std::endl;
        return 1;
    }

    loop start, end;

    ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS];
    ap_uint<576> link_out[6];

    std::string patternFile_dir = "../../../../../data/inputPatternFiles/";

    DIR* dir = opendir(patternFile_dir.c_str());
    if (dir == nullptr) {
        std::cerr << "Error: Unable to open directory" << std::endl;
        return 1;
    }

    dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0 || strstr(entry->d_name, ".links") == 0)
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
        std::getline(read, line);

        std::istringstream word(line);
        std::string link;
        int index = 0;
        while (word >> link && index < N_INPUT_LINKS) {
            link_in[index] = ap_uint<LINK_WIDTH>(link.c_str(), 2);
            index++;
        }

        algo_topIP1(link_in, link_out);

        std::cout << std::endl;
        for (loop i = 0; i < 10; i++) {
            std::cout << std::hex << (link_out[i]) << std::dec << std::endl;
        }
        std::cout << std::endl;

        for (loop i = 0; i < 10; i++) {
            for (loop j = 0; j < 576; j++) {
                std::cout << (link_out[i][j]);
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;

        for (loop i = 0; i < 6; i++) {
            ap_uint<12> start = 0;
            for (loop k = 0; k < N_PUPPI_LINK; k++) {
                l1ct::PuppiObj puppiOutObj = l1ct::PuppiObj::unpack(link_out[i].range(start + 63, start));
                start = start + 64;

                std::cout << "OUT : " << puppiOutObj.hwPt << ",";
                std::cout << puppiOutObj.hwEta << "," << puppiOutObj.hwPhi << "\n" << start << std::endl;
            }
        }

        std::string base_name = entry->d_name;
        size_t pos = base_name.find(".links");
        if (pos != std::string::npos) {
            base_name = base_name.substr(0, pos);
            std::string des_path = patternFile_dir + base_name + ".des";

            std::ifstream des_file(des_path);
            if (des_file.is_open()) {
                std::cout << "\nTrue Clusters of " << base_name << ":\n";
                while (std::getline(des_file, line)) {
                    std::cout << line << std::endl;
                }
                std::cout << std::endl;
                des_file.close();
            } else {
                std::cout << "Warning: No corresponding .des file found for " << entry->d_name << std::endl;
            }
        }
	std::cout << "\n@@@NEW@@@\n";
    }

    closedir(dir);

    
    return 0;
}

