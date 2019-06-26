// uovie iostream
#ifndef IOSTREAM_H_
#define IOSTREAM_H_

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "mol_geom.h"

namespace uovie {
namespace Global {

    class process {
    public:
        std::string filename;
        std::ifstream in;
        std::ofstream out;
        system sys;

        void open(const std::string& file_name) {
            // check the extension
            filename = file_name;
            if (filename.rfind(".xyz") == std::string::npos)
                throw std::invalid_argument("Only .xyz file is a valid input file.");

            in.open(filename);

            if (in.fail()) {
                std::cout << "Can not open the file " << filename << '.' << std::endl;
                exit(EXIT_FAILURE);
            }

            std::string filename_no_extension;
            for (auto it = filename.begin(); it != filename.end() - 4; it++)
                filename_no_extension += *it;

            // open output file
            out.open(filename_no_extension + ".out");
        }

        void read();
        void print();

        void close() {
            // close files
            in.close();
            out.close();
        }

    };
    
}   // Global
}   // uovie

#endif