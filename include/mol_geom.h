// Molecular geometry
#ifndef MOL_GEOM_H_
#define MOL_GEOM_H_

#include <string>
#include <vector>

namespace uovie {
namespace Global {

    class atom {
    public:
        std::string symbol;                 // atomic symbol
        int atomic_number;                  // atomic number
        double m;                           // atomic mass
        std::vector<double> q;              // atomic position
        std::vector<double> p;              // atomic momentum
        std::vector<double> F;              // force acting on atom
    };

    class molecule {
    public:
        std::vector<atom> atoms;
    };

    class system {
    public:
        int dimension;
        int num_part;
        double volume;
        double temperature;
        double pressure;
        std::string model_type;
        std::vector<double> model_para;
        std::vector<molecule> molecules;
    };
}   // Global
}   // uovie

#endif