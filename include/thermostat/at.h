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
#include "thermostat/thermostat.h"

namespace uovie {
namespace thermostat {
namespace at {

    /*** ================================================== ***/
    /*** AT Procedure (Base)                                ***/
    /*** ================================================== ***/

    class at_base : public thermostat_base {
    public:
        at_base() = default;
        at_base(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const double& _nu):
            thermostat_base(_fn_no_ex, _bsp, _sys), nu(_nu) { }

    protected:
        const double& nu;

        Eigen::ArrayXd nrand;
        Eigen::ArrayXd urand;

        std::vector<std::mt19937> nmtes;
        std::vector<std::normal_distribution<double>> nds;
        std::vector<std::mt19937> umtes;
        std::vector<std::uniform_real_distribution<double>> urds;

        void initialize() override;
        void upd_uni_rand_num();
    };

    /*** ================================================== ***/
    /*** AT Procedure (Side)                                ***/
    /*** ================================================== ***/

    class at_side : public at_base {
    public:
        at_side() = default;
        at_side(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const double& _nu):
            at_base(_fn_no_ex, _bsp, _sys, _nu)
        {
            labels.push_back("AT");
            labels.push_back("Side");
        }

    private:
        const double cri = 1 - exp(-nu * Dt / 2);
        void implement_one_step() override;
    };

    /*** ================================================== ***/
    /*** AT Procedure (Middle)                              ***/
    /*** ================================================== ***/

    class at_middle : public at_base {
    public:
        at_middle() = default;
        at_middle(const std::string& _fn_no_ex, const Global::basic_simu_para& _bsp,
            const Global::system& _sys, const double& _nu) :
            at_base(_fn_no_ex, _bsp, _sys, _nu)
        {
            labels.push_back("AT");
            labels.push_back("Middle");
        }

    private:
        const double cri = 1 - exp(-nu * Dt);
        void implement_one_step();
    };

} // !at
} // !thermostat
} // !uovie

#endif // !ANDERSEN_THERMOSTAT_H_