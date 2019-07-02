/* Path Integral Molecular Dynamics */
#ifndef PIMD_H_
#define PIMD_H_

#include <vector>

#include "mol_geom.h"
#include "thermostat/nhc.h"

namespace uovie {
namespace pimd {

    class basic_vari {
    public:
        int num_time_inter;                 // number of time interval

    };


    class pimd_procedure {
    public:
        pimd_procedure() = default;
        pimd_procedure(Global::system& s) : sys(s) { }

        Global::system sys;


        // void implement(std::ofstream& out);

    private:
        int& d = sys.dimension;
        int& N = sys.num_part;
        double k = phy_const::Boltzmann_const;
        double& T = sys.temperature;

        // inverse staging transformayion



        void do_nhc(const int& bi);



    };


}   // pimd
}   // uovie


#endif