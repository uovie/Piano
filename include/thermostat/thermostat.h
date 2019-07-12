/* Thermostat */
#ifndef THERMOSTAT_H_
#define THERMOSTAT_H_

// standard C++ headers
#include <fstream>
#include <vector>
#include <chrono>

// Eigen matrix algebra library
#include <Eigen/Dense>

// uovie headers
#include "simu_para.h"
#include "phy_const.h"
#include "mol_geom.h"

namespace uovie {
namespace thermostat {

    // shorthand for physical constants
    constexpr double h_bar = uovie::phy_const::red_Planck_const;
    constexpr double k = uovie::phy_const::Boltzmann_const;

    /*** ================================================== ***/
    /*** Thermostat (Base)                                ***/
    /*** ================================================== ***/

    class thermostat_base {
    public:
        thermostat_base() = default;
        thermostat_base(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys) : fn_no_ex(_fn_no_ex), bsp(_bsp), sys(_sys) { }

        virtual void implement();

    protected:
        const std::string& fn_no_ex;
        const Global::basic_simu_para& bsp;
        const Global::system& sys;

        const double& Dt = bsp.step_size;
        const int& d = sys.dimension;
        const int& N = sys.num_part;
        const int dof = d * N;
        const double& T = sys.temperature;
        const double beta = 1 / (k * T);

        Eigen::ArrayXd m;          // mass
        Eigen::ArrayXd q;          // position
        Eigen::ArrayXd p;          // momentum
        Eigen::ArrayXd F;          // force

        double kin_ene = 0;
        double pot_ene = 0;

        virtual void initialize() = 0;
        void calc_physic_force();
        virtual void calc_physic_energy();
        virtual void implement_one_step() = 0;

        std::vector<std::string> labels;

        virtual void print_ther_proce_title(std::ofstream& chk, std::ofstream& out);
        virtual void print_ther_proce_data(std::ofstream& chk, std::ofstream& out, double& t);
        virtual void print_conclusion_info(std::ofstream& chk, std::ofstream& out,
            const std::chrono::duration<double>& time_elap);
    };


} // !thermostat
} // !uovie

#endif // !THERMOSTAT_H_