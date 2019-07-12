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
            thermostat::ld::ld_side ld_proce(simulation.fn_no_ex, simulation.bsp, simulation.sys, 1);
            ld_proce.implement();
        }
        else if (simulation.des[0] == "middle") {
            thermostat::ld::ld_middle ld_proce(simulation.fn_no_ex, simulation.bsp, simulation.sys, 1);
            ld_proce.implement();
        }
    }
    else if (simulation.job == "at") {
        if (simulation.des[0] == "side") {
            thermostat::at::at_side at_proce(simulation.fn_no_ex, simulation.bsp, simulation.sys, 1.4);
            at_proce.implement();
        }
        else if (simulation.des[0] == "middle") {
            thermostat::at::at_middle at_proce(simulation.fn_no_ex, simulation.bsp, simulation.sys, 1.4);
            at_proce.implement();
        }

    }
    else if (simulation.job == "nhc") {

        thermostat::nhc::thermo_factor_scheme tfs(7, 1);
        if (simulation.des[0] == "global") {
            if (simulation.des[1] == "side") {
                thermostat::nhc::nhc_global_side nhc_proce(simulation.fn_no_ex, simulation.bsp, simulation.sys, tfs, 4);
                nhc_proce.implement();
            }
            else if (simulation.des[1] == "middle") {
                thermostat::nhc::nhc_global_middle nhc_proce(simulation.fn_no_ex, simulation.bsp, simulation.sys, tfs, 4);
                nhc_proce.implement();
            }
        }
        else if (simulation.des[0] == "local") {
            if (simulation.des[1] == "side") {
                thermostat::nhc::nhc_local_side nhc_proce(simulation.fn_no_ex, simulation.bsp, simulation.sys, tfs, 4);
                nhc_proce.implement();
            }
            else if (simulation.des[1] == "middle") {
                thermostat::nhc::nhc_local_middle nhc_proce(simulation.fn_no_ex, simulation.bsp, simulation.sys, tfs, 4);
                nhc_proce.implement();
            }
        }

    }
    else if (simulation.job == "pimd") {

        if (simulation.des[0] == "ld") {
            if (simulation.des[1] == "side") {
                pimd::pimd_via_ld_side pimd_proce(simulation.fn_no_ex, simulation.bsp, simulation.sys, 8, 1);
                pimd_proce.implement();
            }
            if (simulation.des[1] == "middle") {
                pimd::pimd_via_ld_middle pimd_proce(simulation.fn_no_ex, simulation.bsp, simulation.sys, 8, 2.41888433e-5);
                pimd_proce.implement();
            }
        }
        else if (simulation.des[0] == "at") {
            if (simulation.des[1] == "side") {
                pimd::pimd_via_at_side pimd_proce(simulation.fn_no_ex, simulation.bsp, simulation.sys, 8, 1);
                pimd_proce.implement();
            }
            if (simulation.des[1] == "middle") {
                pimd::pimd_via_at_middle pimd_proce(simulation.fn_no_ex, simulation.bsp, simulation.sys, 8, 1);
                pimd_proce.implement();
            }
        }
        else if (simulation.des[0] == "nhc") {
            thermostat::nhc::thermo_factor_scheme tfs(7, 1);
            if (simulation.des[1] == "side") {
                pimd::pimd_via_nhc_side pimd_proce(simulation.fn_no_ex, simulation.bsp, simulation.sys, 8, tfs, 4);
                pimd_proce.implement();
            }
            else if (simulation.des[1] == "middle") {
                pimd::pimd_via_nhc_middle pimd_proce(simulation.fn_no_ex, simulation.bsp, simulation.sys, 8, tfs, 4);
                pimd_proce.implement();
            }
        }
    }

    return 0;
}