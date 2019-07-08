// standard C++ headers
#include <string>
#include <vector>
#include <cmath>

// uovie headers
#include "process.h"
#include "simu_para.h"
#include "phy_const.h"
#include "thermostat/nhc.h"
#include "pimd.h"

using namespace uovie;

uovie::Global::process simulation;

int main(int argc, char* argv[])
{
    simulation.open(argv[1]);
    simulation.read();

    if (simulation.job == "nhc") {

        thermostat::nhc::thermo_factor_scheme tfs(7, 1);
        if (simulation.des[0] == "global") {
            thermostat::nhc::nhc_procedure_global nhc_proce(simulation.bsp, simulation.sys, tfs, 4);
            nhc_proce.implement(simulation.out);
        }
        else if (simulation.des[0] == "local") {
            thermostat::nhc::nhc_procedure_local nhc_proce(simulation.bsp, simulation.sys, tfs, 4);
            nhc_proce.implement(simulation.out);
        }

    }
    else if (simulation.job == "pimd") {

        if (simulation.des[0] == "nhc") {
            pimd::pimd_via_nhc pimd_proce(simulation.out, simulation.bsp, simulation.sys, 8, 4);
            pimd_proce.implement();
        }

    }

    simulation.print();
    simulation.close();

    return 0;
}