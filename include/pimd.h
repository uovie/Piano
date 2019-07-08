/* Path Integral Molecular Dynamics */
#ifndef PIMD_H_
#define PIMD_H_

// standard C++ headers
#include <fstream>

// uovie headers
#include "simu_para.h"
#include "mol_geom.h"

namespace uovie {
namespace pimd {

    // shorthand for physical constants
    constexpr double h_bar = uovie::phy_const::red_Planck_const;
    constexpr double k = uovie::phy_const::Boltzmann_const;
    
    class pimd_via_nhc {
    public:
        pimd_via_nhc() = default;
        pimd_via_nhc(std::ofstream& _out, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const int _nbead, const int _nchain) :
            out(_out), bsp(_bsp), sys(_sys), nbead(_nbead), nchain(_nchain) { }

        void implement();

    private:
        std::ofstream& out;
        const Global::basic_simu_para& bsp;
        const Global::system& sys;
        const int nbead;
        const int nchain;
    };

} // !pimd
} // !uovie
#endif // !PIMD_H_