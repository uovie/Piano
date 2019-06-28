// uovie iostream
#ifndef IOSTREAM_H_
#define IOSTREAM_H_

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "simu_para.h"
#include "mol_geom.h"

namespace uovie {
namespace Global {

    class process {
    public:
        std::string filename;
        std::ifstream in;
        std::ofstream out;
        std::string job;
        basic_simu_para bsp;
        system sys;

        void open(const std::string& file_name);
        void read();
        void print();
        void close();
    };
    
}   // Global
}   // uovie
#endif