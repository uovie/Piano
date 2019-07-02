#include "..//include/pimd.h"


namespace uovie {
namespace pimd {

    void pimd_procedure::do_nhc(const int& bi)
    {
        // Thermostat variables Declaration
        std::vector<thermostat::nhc::thermo_vari> tmvs;
        tmvs.push_back({ static_cast<double>(d * N), T, 0, 1 });
        tmvs.push_back({ 1, T, 0, -1 });
        tmvs.push_back({ 1, T, 0, 1 });
        tmvs.push_back({ 1, T, 0, -1 });

        // Thermostat Factorization Scheme Declaration
        thermostat::nhc::thermo_factor_scheme tfs(7, 1);

        // NHC Simulation Execution
        thermostat::nhc::nhc_procedure_base nhc_proce(simulation.bsp,
            simulation.sys, tmvs, tfs, bi);
        nhc_proce.implement(simulation.out);

    }

}   // pimd
}   // uovie