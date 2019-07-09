// standard C++ headers
#include <string>
#include <vector>
#include <cmath>

// uovie headers
#include "process.h"
#include "simu_para.h"
#include "thermostat/ld.h"
#include "thermostat/at.h"
#include "thermostat/nhc.h"
#include "pimd.h"

using namespace uovie;

uovie::Global::process simulation;

int main(int argc, char* argv[])
{
    simulation.open(argv[1]);
    simulation.read();

    if (simulation.job == "ld") {
        if (simulation.des[0] == "side") {
            thermostat::ld::ld_procedure_side ld_proce(simulation.bsp, simulation.sys, 1);
            ld_proce.implement(simulation.out);
        }
        else if (simulation.des[0] == "middle") {
            thermostat::ld::ld_procedure_middle ld_proce(simulation.bsp, simulation.sys, 1);
            ld_proce.implement(simulation.out);
        }
    }
    else if (simulation.job == "at") {

        if (simulation.des[0] == "side") {
            thermostat::at::at_procedure_side at_proce(simulation.bsp, simulation.sys, 1.4);
            at_proce.implement(simulation.out);
        }
        else if (simulation.des[0] == "middle") {
            thermostat::at::at_procedure_middle at_proce(simulation.bsp, simulation.sys, 1.4);
            at_proce.implement(simulation.out);
        }

    }
    else if (simulation.job == "nhc") {

        thermostat::nhc::thermo_factor_scheme tfs(7, 1);
        if (simulation.des[0] == "global") {
            if (simulation.des[1] == "side") {
                thermostat::nhc::nhc_procedure_global_side nhc_proce(simulation.bsp, simulation.sys, tfs, 4);
                nhc_proce.implement(simulation.out);
            }
        }
        else if (simulation.des[0] == "local") {
            if (simulation.des[1] == "side") {
                thermostat::nhc::nhc_procedure_local_side nhc_proce(simulation.bsp, simulation.sys, tfs, 4);
                nhc_proce.implement(simulation.out);
            }
        }

    }
    else if (simulation.job == "pimd") {

        if (simulation.des[0] == "ld") {
            pimd::pimd_via_ld pimd_proce(simulation.out, simulation.bsp, simulation.sys, 8);
            pimd_proce.implement();
        }
        else if (simulation.des[0] == "nhc") {
            pimd::pimd_via_nhc pimd_proce(simulation.out, simulation.bsp, simulation.sys, 8, 4);
            pimd_proce.implement();
        }

    }

    simulation.print();
    simulation.close();

    return 0;
}