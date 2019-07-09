/* Langevin Dynamics */
#ifndef LANGEVIN_DYNAMICS_H_
#define LANGEVIN_DYNAMICS_H_

// standard C++ headers
#include <fstream>
#include <vector>
#include <cmath>
#include <random>

// Eigen matrix algebra library
#include <Eigen/Dense>

// uovie headers
#include "simu_para.h"
#include "phy_const.h"
#include "mol_geom.h"

namespace uovie {
namespace thermostat {
namespace ld {

    // shorthand for physical constants
    constexpr double h_bar = uovie::phy_const::red_Planck_const;
    constexpr double k = uovie::phy_const::Boltzmann_const;

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** LD Procedure (Side)                                ***/
    /*** ================================================== ***/

    class ld_procedure_side {
    public:
        ld_procedure_side() = default;
        ld_procedure_side(const Global::basic_simu_para& _bsp, const Global::system& _sys,
            const double _gamma) : bsp(_bsp), sys(_sys), gamma(_gamma) { }

        void implement();
        void implement(std::ofstream& out);

    private:
        const Global::basic_simu_para& bsp;
        const Global::system& sys;
        const double& gamma;    // friction coefficient

        const double& Dt = bsp.step_size;
        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const int dof = d * N;
        const double& T = sys.temperature;
        const double beta = 1 / (k * T);
        const double c1 = exp(-gamma * Dt / 2);
        const double c2 = sqrt(1 - pow(c1, 2));

        Eigen::ArrayXd m;          // mass
        Eigen::ArrayXd q;          // position
        Eigen::ArrayXd p;          // momentum
        Eigen::ArrayXd F;          // force
        Eigen::ArrayXd nrand;

        std::vector<std::mt19937> mtes;
        std::vector<std::normal_distribution<double>> nds;

        double kine_energy = 0;
        double pote_energy = 0;
        double ther_energy = 0;
        double cons_energy = 0;

        void initialize();
        void calc_physic_force();
        void upd_rand_num();

        void print_ld_procedure_title(std::ofstream& out);
        void print_ld_procedure_data(std::ofstream& out, double& t);
        void implement_one_step();
    };

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** LD Procedure (middle)                                ***/
    /*** ================================================== ***/

    class ld_procedure_middle {
    public:
        ld_procedure_middle () = default;
        ld_procedure_middle(const Global::basic_simu_para& _bsp, const Global::system& _sys,
            const double& _gamma) : bsp(_bsp), sys(_sys), gamma(_gamma) { }

        void implement();
        void implement(std::ofstream& out);

    private:
        const Global::basic_simu_para& bsp;
        const Global::system& sys;
        const double& gamma;    // friction coefficient

        const double& Dt = bsp.step_size;
        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const int dof = d * N;
        const double& T = sys.temperature;
        const double beta = 1 / (k * T);
        const double c1 = exp(-gamma * Dt);
        const double c2 = sqrt(1 - pow(c1, 2));

        Eigen::ArrayXd m;          // mass
        Eigen::ArrayXd q;          // position
        Eigen::ArrayXd p;          // momentum
        Eigen::ArrayXd F;          // force
        Eigen::ArrayXd nrand;

        std::vector<std::mt19937> mtes;
        std::vector<std::normal_distribution<double>> nds;

        double kine_energy = 0;
        double pote_energy = 0;
        double ther_energy = 0;
        double cons_energy = 0;

        void initialize();
        void calc_physic_force();
        void upd_rand_num();

        void print_ld_procedure_title(std::ofstream& out);
        void print_ld_procedure_data(std::ofstream& out, double& t);
        void implement_one_step();
    };

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** LD Procedure for PIMD (Middle)                     ***/
    /*** ================================================== ***/

    class ld_procedure_for_pimd {
    public:
        ld_procedure_for_pimd() = default;
        ld_procedure_for_pimd(const Global::basic_simu_para& _bsp, const Global::system& _sys,
            const double& _gamma, const int& _nbead) :
            bsp(_bsp), sys(_sys), gamma(_gamma), nbead(_nbead) { }

        void implement();
        void implement(std::ofstream& out);

    private:
        const Global::basic_simu_para& bsp;
        const Global::system& sys;
        const double& gamma;    // friction coefficient
        const int& nbead;

        const double& Dt = bsp.step_size;
        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const int dof = d * N;
        const double& T = sys.temperature;
        const double beta = 1 / (k * T);
        const double c1 = exp(-gamma * Dt);
        const double c2 = sqrt(1 - pow(c1, 2));
        const double fic_omega = sqrt(nbead) / (h_bar * beta);

        Eigen::ArrayXXd m;          // real mass
        Eigen::ArrayXXd q;          // real position
        Eigen::ArrayXXd m_tilde;    // transformed mass tilde
        Eigen::ArrayXXd m_bar;      // transformed mass bar
        Eigen::ArrayXXd r;          // transformed position
        Eigen::ArrayXXd s;          // fictitious momentum
        Eigen::ArrayXXd F;          // fictitious force
        Eigen::ArrayXXd nrand;

        std::vector<std::mt19937> mtes;
        std::vector<std::normal_distribution<double>> nds;

        Eigen::ArrayXXd kine_energy;
        Eigen::ArrayXXd pote_energy;

        double prim_kine_estor = 0;
        double prim_pote_estor = 0;
        double prim_pres_estor = 0;

        void initialize();
        void stag_trans();
        void inve_stag_trans();
        void calc_physic_force();
        void upd_rand_num();
        void calc_prim_estor();

        void print_ld_procedure_title(std::ofstream& out);
        void print_ld_procedure_data(std::ofstream& out, double& t);
        void implement_one_step();
    };

} // !ld
} // !thermostat
} // !uovie

#endif // !LANGEVIN_DYNAMICS_H_