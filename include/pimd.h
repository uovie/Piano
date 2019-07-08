/* Path Integral Molecular Dynamics */
#ifndef PIMD_H_
#define PIMD_H_

// standard C++ headers
#include <vector>

// Piano headers
#include "phy_const.h"
#include "mol_geom.h"
#include "thermostat/nhc.h"



namespace uovie {
namespace pimd {
        
    class prim_sys_coll {
    public:
        prim_sys_coll() = default;
        prim_sys_coll(std::vector<Global::system>& ss): syss(ss) { }

        std::vector<Global::system> syss;
        std::vector<Global::system> stag_trans();       // staging transformation
    };

    class stra_sys_coll {
    public:
        stra_sys_coll() = default;
        stra_sys_coll(std::vector<Global::system>& ss) : syss(ss) { }

        std::vector<Global::system> syss;
        std::vector<Global::system> inve_stag_trans();  // inverse staging transformation
    };

    class pimd_procedure {
    public:
        pimd_procedure() = default;
        pimd_procedure(const Global::basic_simu_para& b, Global::system& s,
            const int nb) : bsp(b), sys(s), nbead(nb) {
            for (int bi = 0; bi < nbead; bi++)
                prsc.syss.push_back(s);
            stsc.syss = prsc.stag_trans();
        }

        void implement();
        void implement(std::ofstream& out, std::string& fn_no_ex);

    private:
        const Global::basic_simu_para& bsp;
        Global::system& sys;
        int nbead;
        prim_sys_coll prsc;  // primitive systerm collection
        stra_sys_coll stsc;  // staging transformed systerm collection

        int& d = sys.dimension;
        int& N = sys.num_part;
        const double k = uovie::phy_const::Boltzmann_const;
        double& T = sys.temperature;

        double tot_sys_energy = 0;
        double cano_prob_dens = 0;  // canonical probability density

        void print_pimd_cpd_title(std::ofstream& out);
        void print_pimd_cpd_data(std::ofstream& out, double& t);
    };

} // !pimd
} // !uovie
#endif