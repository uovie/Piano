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
#include "thermostat/thermostat.h"

namespace uovie {
namespace thermostat {
namespace ld {

    /*** ================================================== ***/
    /*** LD Procedure (Base)                                ***/
    /*** ================================================== ***/

    class ld_base : public thermostat_base {
    public:
        ld_base() = default;
        ld_base(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const double _gamma):
            thermostat_base(_fn_no_ex, _bsp, _sys), gamma(_gamma) { }

    protected:
        const double& gamma;    // friction coefficient

        Eigen::ArrayXd nrand;

        std::vector<std::mt19937> mtes;
        std::vector<std::normal_distribution<double>> nds;

        void initialize() override;
        void upd_rand_num();
    };
    
    /*** ================================================== ***/
    /*** LD Procedure (Side)                                ***/
    /*** ================================================== ***/

    class ld_side : public ld_base {
    public:
        ld_side() = default;
        ld_side(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const double _gamma):
            ld_base(_fn_no_ex, _bsp, _sys, _gamma)
        {
            labels.push_back("LD");
            labels.push_back("Side");
        }

    private:
        const double c1 = exp(-gamma * Dt / 2);
        const double c2 = sqrt(1 - pow(c1, 2));

        void implement_one_step() override;
    };

    /*** ================================================== ***/
    /*** LD Procedure (middle)                              ***/
    /*** ================================================== ***/

    class ld_middle : public ld_base {
    public:
        ld_middle () = default;
        ld_middle(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const double _gamma) :
            ld_base(_fn_no_ex, _bsp, _sys, _gamma)
        {
            labels.push_back("LD");
            labels.push_back("Middle");
        }

    private:
        const double c1 = exp(-gamma * Dt);
        const double c2 = sqrt(1 - pow(c1, 2));

        void implement_one_step() override;
    };

} // !ld
} // !thermostat
} // !uovie

#endif // !LANGEVIN_DYNAMICS_H_