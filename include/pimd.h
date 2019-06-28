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

    class ficti_vari {
    public:
        std::vector<double> m_tilde;
        std::vector<double> r;
        std::vector<double> s;
        std::vector<double> m_bar;
        
    };

    class pimd_procedure {
    public:
        pimd_procedure() = default;
        pimd_procedure(Global::system& s) : sys(s) { }

        Global::system sys;
        std::vector<ficti_vari> fvs;


        // void implement(std::ofstream& out);

    private:
        // inverse staging transformayion


        // nhc
        

        /*** ================================================== ***/
        /*** Thermostat variables Declaration                   ***/
        /*** ================================================== ***/

        std::vector<thermostat::nhc::thermo_vari> tmvs;
        tmvs.push_back({ static_cast<double>(d * N), T, 0, 1 });
        tmvs.push_back({ 1, T, 0, -1 });
        tmvs.push_back({ 1, T, 0, 1 });
        tmvs.push_back({ 1, T, 0, -1 });

        /*** ================================================== ***/
        /*** Thermostat Factorization Scheme Declaration        ***/
        /*** ================================================== ***/

        thermostat::nhc::thermo_factor_scheme tfs(7, 1);
        tfs.init();

        /*** ================================================== ***/
        /*** NHC Simulation Execution                           ***/
        /*** ================================================== ***/

        thermostat::nhc::nhc_procedure nhc_proce(simulation.bsp,
            simulation.sys, tmvs, tfs);
        nhc_proce.implement(simulation.out);




    };


}   // pimd
}   // uovie


#endif