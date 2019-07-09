/* Andersen Thermostat */
#ifndef ANDERSEN_THERMOSTAT_H_
#define ANDERSEN_THERMOSTAT_H_

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
namespace at {

    // shorthand for physical constants
    constexpr double h_bar = uovie::phy_const::red_Planck_const;
    constexpr double k = uovie::phy_const::Boltzmann_const;

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** AT Procedure (Side)                                ***/
    /*** ================================================== ***/

    class at_procedure_side {
    public:
        at_procedure_side() = default;
        at_procedure_side(const Global::basic_simu_para& _bsp, const Global::system& _sys,
            const double& _nu) : bsp(_bsp), sys(_sys), nu(_nu) { }

        void implement();
        void implement(std::ofstream& out);

    private:
        const Global::basic_simu_para& bsp;
        const Global::system& sys;
        const double& nu;

        const double& Dt = bsp.step_size;
        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const int dof = d * N;
        const double& T = sys.temperature;
        const double beta = 1 / (k * T);
        const double cri = 1 - exp(-nu * Dt / 2);

        Eigen::ArrayXd m;          // mass
        Eigen::ArrayXd q;          // position
        Eigen::ArrayXd p;          // momentum
        Eigen::ArrayXd F;          // force
        Eigen::ArrayXd nrand;
        Eigen::ArrayXd urand;

        std::vector<std::mt19937> nmtes;
        std::vector<std::normal_distribution<double>> nds;
        std::vector<std::mt19937> umtes;
        std::vector<std::uniform_real_distribution<double>> urds;

        double kine_energy = 0;
        double pote_energy = 0;
        double ther_energy = 0;
        double cons_energy = 0;

        void initialize();
        void calc_physic_force();
        void upd_uni_rand_num();

        void print_at_procedure_title(std::ofstream& out);
        void print_at_procedure_data(std::ofstream& out, double& t);
        void implement_one_step();
    };

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** AT Procedure (Middle)                              ***/
    /*** ================================================== ***/

    class at_procedure_middle {
    public:
        at_procedure_middle() = default;
        at_procedure_middle(const Global::basic_simu_para& _bsp, const Global::system& _sys,
            const double& _nu) : bsp(_bsp), sys(_sys), nu(_nu) { }

        void implement();
        void implement(std::ofstream& out);

    private:
        const Global::basic_simu_para& bsp;
        const Global::system& sys;
        const double& nu;

        const double& Dt = bsp.step_size;
        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const int dof = d * N;
        const double& T = sys.temperature;
        const double beta = 1 / (k * T);
        const double cri = 1 - exp(-nu * Dt);

        Eigen::ArrayXd m;          // mass
        Eigen::ArrayXd q;          // position
        Eigen::ArrayXd p;          // momentum
        Eigen::ArrayXd F;          // force
        Eigen::ArrayXd nrand;
        Eigen::ArrayXd urand;

        std::vector<std::mt19937> nmtes;
        std::vector<std::normal_distribution<double>> nds;
        std::vector<std::mt19937> umtes;
        std::vector<std::uniform_real_distribution<double>> urds;

        double kine_energy = 0;
        double pote_energy = 0;
        double ther_energy = 0;
        double cons_energy = 0;

        void initialize();
        void calc_physic_force();
        void upd_uni_rand_num();

        void print_at_procedure_title(std::ofstream& out);
        void print_at_procedure_data(std::ofstream& out, double& t);
        void implement_one_step();
    };

    //--------------------------------------------------------//

    /*** ================================================== ***/
    /*** AT Procedure for PIMD (Middle)                     ***/
    /*** ================================================== ***/

    class at_procedure_for_pimd {
    public:
        at_procedure_for_pimd() = default;
        at_procedure_for_pimd(const Global::basic_simu_para& _bsp, const Global::system& _sys,
            const double& _nu, const int& _nbead) :
            bsp(_bsp), sys(_sys), nu(_nu), nbead(_nbead) { }

        void implement();
        void implement(std::ofstream& out);

    private:
        const Global::basic_simu_para& bsp;
        const Global::system& sys;
        const double& nu;
        const int& nbead;

        const double& Dt = bsp.step_size;
        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const int dof = d * N;
        const double& T = sys.temperature;
        const double beta = 1 / (k * T);
        const double cri = 1 - exp(-nu * Dt);
        const double fic_omega = sqrt(nbead) / (h_bar * beta);

        Eigen::ArrayXXd m;          // real mass
        Eigen::ArrayXXd q;          // real position
        Eigen::ArrayXXd m_tilde;    // transformed mass tilde
        Eigen::ArrayXXd m_bar;      // transformed mass bar
        Eigen::ArrayXXd r;          // transformed position
        Eigen::ArrayXXd s;          // fictitious momentum
        Eigen::ArrayXXd F;          // fictitious force
        Eigen::ArrayXXd nrand;
        Eigen::ArrayXXd urand;

        std::vector<std::mt19937> nmtes;
        std::vector<std::normal_distribution<double>> nds;
        std::vector<std::mt19937> umtes;
        std::vector<std::uniform_real_distribution<double>> urds;

        Eigen::ArrayXXd kine_energy;
        Eigen::ArrayXXd pote_energy;

        double prim_kine_estor = 0;
        double prim_pote_estor = 0;
        double prim_pres_estor = 0;

        void initialize();
        void stag_trans();
        void inve_stag_trans();
        void calc_physic_force();
        void upd_uni_rand_num();
        void calc_prim_estor();

        void print_at_procedure_title(std::ofstream& out);
        void print_at_procedure_data(std::ofstream& out, double& t);
        void implement_one_step();
    };

} // !at
} // !thermostat
} // !uovie

#endif // !ANDERSEN_THERMOSTAT_H_