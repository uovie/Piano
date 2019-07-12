/* Thermostat */

// standard C++ headers
#include <iostream>
#include <cmath>
#include <iomanip>
#include <chrono>

#include "thermostat/thermostat.h"
#include "model.h"

namespace uovie {
namespace thermostat {

    // calculate physical forces
    void thermostat_base::calc_physic_force()
    {
        if (sys.model_type == "HO") { // simple harmonic forces
            model::harmonic_oscilator HO(sys.model_para[0]);
            F = -1 * m * pow(HO.ome(), 2) * q;
        }
        else if (sys.model_type == "LJ") { // Lennard_Jones forces
            F.setZero();
            model::Lennard_Jones LJ(sys.model_para[0], sys.model_para[1]);
            for (auto ni = 0; ni < N; ni++) {
                for (auto nj = 0; nj < N; nj++) {
                    if (ni == nj) continue;
                    F.block(d * ni, 0, d, 1) += LJ.F(q.block(d * ni, 0, d, 1), q.block(d * nj, 0, d, 1));
                }
            }
        }
    }

    // calculate certain energies
    void thermostat_base::calc_physic_energy()
    {
        // calculate kinetic energy
        kin_ene = (p.pow(2) / (2 * m)).sum();

        // calculate potential energy
        if (sys.model_type == "HO") { // simple harmonic forces
            model::harmonic_oscilator HO(sys.model_para[0]);
            pot_ene = 0.5 * pow(HO.ome(), 2) * (m * q.pow(2)).sum();
        }
        else if (sys.model_type == "LJ") { // Lennard_Jones forces
            model::Lennard_Jones LJ(sys.model_para[0], sys.model_para[1]);
            pot_ene = 0;
            for (auto ni = 0; ni < N - 1; ni++)
                for (auto nj = ni + 1; nj < N; nj++)
                    pot_ene += LJ.V(q.block(d * ni, 0, d, 1), q.block(d * nj, 0, d, 1));
        }
    }

    void thermostat_base::implement()
    {
        initialize();

        double t = 0;
        int ctr = 0;

        std::ofstream chk, out;
        chk.open(fn_no_ex + ".chk");
        out.open(fn_no_ex + ".out");

        print_ther_proce_title(chk, out);
        calc_physic_energy();
        print_ther_proce_data(chk, out, t);

        const auto tstart = std::chrono::high_resolution_clock::now();
        do {
            implement_one_step();
            if (ctr == bsp.data_coll_peri) {
                calc_physic_energy();
                print_ther_proce_data(chk, out, t);
                ctr = 0;
            }
            ctr++;
            t += Dt;
        } while (t <= bsp.run_time);
        const auto tstop = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> time_elapsed = tstop - tstart;

        print_conclusion_info(chk, out, time_elapsed);

        chk.close();
        out.close();
    }

    void thermostat_base::print_ther_proce_title(std::ofstream& chk, std::ofstream& out)
    {
        std::cout << "\n" + labels[0] + " procedure (" + labels[1] + ") is running." << std::endl;
        chk << labels[0] + " procedure (" + labels[1] + "):\n" << "           Time" << "              position"
            << "            momentum" << std::endl;
        out << labels[0] + " procedure (" + labels[1] + "):\n" << "           Time" << "              kin_ene"
            << "             pot_ene" << std::endl;
    }

    void thermostat_base::print_ther_proce_data(std::ofstream& chk, std::ofstream& out, double& t)
    {
        chk << std::scientific << std::setprecision(8) << std::setw(20) << t * uovie::phy_const::a_u_time * 1e15
            << std::setw(20) << q(0) << std::setw(20) << p(0) << std::endl;
        out << std::scientific << std::setprecision(8) << std::setw(20) << t * uovie::phy_const::a_u_time * 1e15
            << std::setw(20) << kin_ene << std::setw(20) << pot_ene << std::endl;
    }

    void thermostat_base::print_conclusion_info(std::ofstream& chk, std::ofstream& out,
        const std::chrono::duration<double>& time_elap) {
        chk << "\nElapsed time of " + labels[0] + " procedure (" + labels[1] + "): " << time_elap.count() << std::endl;
        out << "\nElapsed time of " + labels[0] + " procedure (" + labels[1] + "): " << time_elap.count() << std::endl;
        std::cout << "\nElapsed time of " + labels[0] + " procedure (" + labels[1] + "): " << time_elap.count() << std::endl;
        std::cout << "\nNormal termination. Congratulations!" << std::endl;
    }

} // !thermostat
} // !uovie