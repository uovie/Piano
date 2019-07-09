// standard C++ headers
#include <string>
#include <vector>
#include <cmath>

// uovie headers
#include "process.h"
#include "simu_para.h"
#include "thermostat/ld.h"

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
    /*** LD Simulation Execution                            ***/
    /*** ================================================== ***/

    thermostat::ld::ld_procedure_side ld_proce(simulation.bsp, simulation.sys, 1);
    ld_proce.implement(simulation.out);

    /*** ================================================== ***/
    /*** Simulation termination                             ***/
    /*** ================================================== ***/

    simulation.print();
    simulation.close();

    return 0;
}