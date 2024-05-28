#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "ap_int.h"
#include "algo_topIP1.h"

#define LINK_WIDTH 576

const int ROWS_PER_EVENT = 9;
const ap_uint<8> BW = 64;

void readFile(const std::string& filename, std::vector<std::string>& lines) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open the file " << filename << std::endl;
        exit(1);
    }
    std::string line;
    while (std::getline(infile, line)) {
        if (!line.empty() && line[0] != '#') {
            lines.push_back(line);
        }
    }
    infile.close();
}

void processLine(const std::string& line, int rowIndex, ap_uint<LINK_WIDTH>* link) {
    std::istringstream iss(line);
    std::string token;
    std::vector<std::string> hexValues;
    iss >> token;

    while (iss >> token) {
        if (token.substr(0, 2) == "0x") {
            std::string hexValue;
            iss >> hexValue;
            if (hexValue.substr(0, 2) == "0x") {
                hexValues.push_back(hexValue.substr(2));
            }
        }
    }

    for (int i = 0; i < N_INPUT_LINKS; ++i) {
        if (rowIndex < ROWS_PER_EVENT) {
            ap_uint<64> value = std::stoull(hexValues[i], nullptr, 16);
            link[i].range(BW * (rowIndex + 1) - 1, BW * rowIndex) = value;
        }
    }
}

void resetLinks(ap_uint<LINK_WIDTH>* link) {
    for (int i = 0; i < N_INPUT_LINKS; ++i) {
        link[i] = 0;
    }
}

void compareLinks(ap_uint<LINK_WIDTH>* link_out, ap_uint<LINK_WIDTH>* expected_link_out) {
    for (int i = 0; i < N_INPUT_LINKS; ++i) {
        if (link_out[i] != expected_link_out[i]) {
            std::cout << "Mismatch at link " << i << ":\n";
            std::cout << "Computed: " << link_out[i] << "\n";
            std::cout << "Expected: " << expected_link_out[i] << "\n";
        }
    }
}

void printLinks(ap_uint<LINK_WIDTH>* link_out) {
	l1ct::PuppiObj puppiOutObj;
    PFcluster InPFcluster;
    for (int i = 0; i < N_OUTPUT_LINKS; ++i) {
        int start = 0;
        int end = start + BW;
        std::cout << "outLinks " << i << ": " << link_out[i] << std::endl;
        for (int k = 0; k < N_PUPPI_LINK; ++k) {
            puppiOutObj = l1ct::PuppiObj::unpack(link_out[i].range(end, start));
            InPFcluster.getPFcluster(link_out[i].range(end, start));
            start = start + BW;
            end = start + BW;

            std::cout << "Link : " << i << " cluster : " << k << " IN : " << InPFcluster.ET << ","
                      << InPFcluster.Eta << "," << InPFcluster.Phi << ",  "
                      << " OUT : " << puppiOutObj.hwPt << ","
                      << puppiOutObj.hwEta << "," << puppiOutObj.hwPhi << "\n";
        }
    }
}

int main() {
    std::string input_filename = "../../../../../data/l1HFPos-inputs_0.txt";
    std::string output_filename = "../../../../../data/l1HFPos-outputs_0.txt";

    std::vector<std::string> input_lines;
    std::vector<std::string> output_lines;

    readFile(input_filename, input_lines);
    readFile(output_filename, output_lines);

    int rowIndex = 0;
    ap_uint<LINK_WIDTH> link_in[N_INPUT_LINKS];
    ap_uint<LINK_WIDTH> link_out[N_INPUT_LINKS];
    ap_uint<LINK_WIDTH> expected_link_out[N_INPUT_LINKS];

    resetLinks(link_in);

    for (size_t lineIndex = 0; lineIndex < input_lines.size(); ++lineIndex) {
        processLine(input_lines[lineIndex], rowIndex, link_in);
        processLine(output_lines[lineIndex], rowIndex, expected_link_out);
        rowIndex++;

        if (rowIndex == ROWS_PER_EVENT) {
            std::cout << "Starting the Event processing!" << std::endl;
            algo_topIP1(link_in, link_out);

            printLinks(link_out);

            compareLinks(link_out, expected_link_out);

            rowIndex = 0;
            resetLinks(link_in);
        }
    }

    return 0;
}
