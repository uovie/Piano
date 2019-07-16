// standard C++ headers
#include <string>
#include <vector>
#include <cmath>
#include <cassert>

// uovie headers
#include "process.h"
#include "simu_para.h"
#include "pimd.h"

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
    /*** PIMD Simulation Execution                          ***/
    /*** ================================================== ***/

    pimd::pimd_via_nhc pimd_proce(simulation.out, simulation.bsp, simulation.sys, 8, 4);
    pimd_proce.implement();

    /*** ================================================== ***/
    /*** Simulation termination                             ***/
    /*** ================================================== ***/

    simulation.print();
    simulation.close();

    return 0;

}