// standard C++ headers
#include <string>
#include <vector>
#include <cmath>

// uovie headers
#include "process.h"
#include "simu_para.h"
#include "thermostat/nhc.h"

using namespace uovie;

uovie::Global::process simulation;

int main(int argc, char* argv[])
{
    /*** ================================================== ***/
    /*** Simulation Preparation                             ***/
    /*** ================================================== ***/

    simulation.open(argv[1]);
    simulation.read();

    /*** ================================================== ***/
    /*** Thermostat Factorization Scheme Declaration        ***/
    /*** ================================================== ***/

    thermostat::nhc::thermo_factor_scheme tfs(7, 1);

    /*** ================================================== ***/
    /*** NHC Simulation Execution                           ***/
    /*** ================================================== ***/

    thermostat::nhc::nhc_procedure_local_side nhc_proce(simulation.bsp,
        simulation.sys, tfs, 4);
    nhc_proce.implement(simulation.out);

    /*** ================================================== ***/
    /*** Simulation termination                             ***/
    /*** ================================================== ***/
    
    simulation.print();
    simulation.close();

    return 0;
}