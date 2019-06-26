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
        double q[3];                        // atomic position
        double p[3] = { 1,1,1 };            // atomic momentum
        double F[3];                        // force acting on atom
    };

    class molecule {
    public:
        std::vector<atom> atoms;
    };

    class system {
    public:
        std::vector<molecule> molecules;
    };
}   // Global
}   // uovie

#endif