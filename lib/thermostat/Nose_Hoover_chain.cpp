// Nose-Hoover Chain
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <chrono>

#include "../../include/simu_para.h"
#include "../../include/thermostat/nose_hoover_chain.h"

namespace uovie {
namespace thermostat {
namespace nhc {

    void print_nhc_procedure_title(std::ofstream& out, Global::system& sys, double& t);
    void print_nhc_procedure_data(std::ofstream& out, Global::system& sys, double& t);
    
    void nhc_procedure::implement(std::ofstream& out) {
        init();
        double t = 0;
        int ctr = 0;
        
        // print NHC procedure title
        print_nhc_procedure_title(out, sys, t);
        
        // print NHC procedure data
        print_nhc_procedure_data(out, sys, t);
        const auto tstart = std::chrono::high_resolution_clock::now();

        do {
            // NHC numerical evolution
            thermo_propagate();
            physic_propagate();
            thermo_propagate();

            if (ctr == bsp.data_coll_peri) {
                print_nhc_procedure_data(out, sys, t);
                ctr = 0;
            }
            ctr++;
            t += bsp.time_step_size;
        } while (t <= bsp.run_time);

        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;

        std::cout << "\nElapsed time of NHC procedure (s): " << time_elapsed.count() << std::endl;

        std::cout << "\nNormal termination. Congratulations!" << std::endl;

    }

    // print NHC procedure title
    void print_nhc_procedure_title(std::ofstream& out, Global::system& sys, double& t) {
        std::cout << "\n\nNHC Procedure:\n   Time";
        out << "\n\nNHC Procedure:\nTime";

        for (int mi = 0; mi < sys.molecules.size(); mi++) {
            for (int ai = 0; ai < sys.molecules[mi].atoms.size(); ai++) {
                for (int di = 0; di < sys.dimension; di++) {
                    std::cout << "    m[" << mi + 1 << "].a[" << ai + 1 << "].q[" << di + 1 << "]";
                    out << "\tm[" << mi + 1 << "].a[" << ai + 1 << "].q[" << di + 1 << "]";
                }
                for (int di = 0; di < sys.dimension; di++) {
                    std::cout << "\tm[" << mi + 1 << "].a[" << ai + 1 << "].p[" << di + 1 << "]";
                    out << "\tm[" << mi + 1 << "].a[" << ai + 1 << "].p[" << di + 1 << "]";
                }
            }
        }
    }

    // print NHC procedure data
    void print_nhc_procedure_data(std::ofstream& out, Global::system& sys, double& t) {
        std::cout << "\n" << std::fixed << std::setprecision(5) << std::setw(10) << t;
        out << "\n" << std::fixed << std::setprecision(5) << t;

        std::cout << std::setprecision(8);
        out << std::setprecision(8);
        for (int mi = 0; mi < sys.molecules.size(); mi++) {
            for (int ai = 0; ai < sys.molecules[mi].atoms.size(); ai++) {
                for (int di = 0; di < sys.dimension; di++) {
                    std::cout << std::setw(15) << sys.molecules[mi].atoms[ai].q[di];
                    out << "\t" << sys.molecules[mi].atoms[ai].q[di];
                }
                for (int di = 0; di < sys.dimension; di++) {
                    std::cout << std::setw(15) << sys.molecules[mi].atoms[ai].p[di];
                    out << "\t" << sys.molecules[mi].atoms[ai].p[di];
                }
            }
        }
    }

}   // uovie
}   // thermostat
}   // uovie