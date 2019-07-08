// Physical constant
#ifndef PHY_CONST_H_
#define PHY_CONST_H_

#include <string>

namespace uovie {
namespace phy_const {

    /*
     * CODATA Recommended Values of the Fundamental Physical Constants: 2014
     *
     * Comment format:
     *     (quantity, symbol, numerical value, unit, relative std. uncert.)
    */

    /* I
     * The CODATA recommended values of the fundamental constants of physics and chemistry based on the 2014 adjustment.
     */

    /*** ============================ ***/
    /*** Universal                    ***/
    /*** ============================ ***/

    constexpr auto light_speed_vac = 2.99792458e8; // ("speed of light in vacuum", "c, c0", 299792458, "m s−1", 0)
    constexpr auto magnet_const = 1.2566370614e-6; // ("magnetic constant", "mu0", 1.2566370614e-6, "N A−2", 0)
    constexpr auto elec_const = 8.854187817e-12; // ("electric constant", "epsilon0", 8.854187817e-12, "F m−1", 0)
    constexpr auto chara_imped_vac = 376.730313461; // ("characteristic impedance of vacuum", "Z0", 376.730313461, "Omega", 0)
    constexpr auto Newton_gravvi_const = 6.67408e-11; // ("Newtonian constant of gravitation", "G", 6.67408e-11, "m3 kg−1 s−2", 4.7e-5)
    constexpr auto Planck_const = 6.626070040e-34;// ("Planck constant", "h", 6.626070040e-34, "J s", 1.2e-8)
    constexpr auto red_Planck_const = 1/*1.054571800e-34*/; // ("reduced Planck constant", "h_bar", 1.054571800e-34, "J s", 1.2e-8)
    constexpr auto Planck_mass = 2.176470e-8; // ("Planck mass", "mp", 2.176470e-8, "kg", 2.3e-5)
    constexpr auto Planck_ener = 1.220910e19; // ("Planck energy", "mpc2", 1.220910e19, "GeV", 2.3e-5)
    constexpr auto Planck_temp = 1.416808e32; // ("Planck temperature", "Tp", 1.416808e32, "K", 2.3e-5)
    constexpr auto Planck_leng = 1.616229e-35; // ("Planck length", "lp", 1.616229e-35, "m", 2.3e-5)
    constexpr auto Planck_time = 5.39116e-44; // ("Planck time", "tp", 5.39116e-44, "s", 2.3e-5)

    /*** ============================ ***/
    /*** Electromagnetic              ***/
    /*** ============================ ***/



    /*** ============================ ***/
    /*** Physicochemistry             ***/
    /*** ============================ ***/

    constexpr auto Avogadro_const = 6.022140857e23; // ("Avogadro constant", "NA", 6.022140857e23, "mol-1", 1.2e-8)
    constexpr auto atomic_mass_const = 1.660539040e-27; // ("atomic mass constant", "mu", 1.660539040e-27, "kg", 1.2e-8)
    constexpr auto atomic_ener_const = 1.492418062e-10; // ("atomic energy constant", "mu", 1.492418062e-10, "J", 1.2e-8)
    constexpr auto Faraday_const = 96485.33289; // ("Faraday constant", "F", 96485.33289, "C mol-1" 6.2e-9)
    /* ... */
    constexpr auto Boltzmann_const = 1/*1.38064852e-23*/; // ("Boltzmann constant", "k", 1.38064852e-23, "J K-1", 5.7e-7);


    /* V
     * The values in SI units of some non-SI units based on the 2014 CODATA adjustment of the values of the constants.
     */

    /*** ========================================= ***/
    /*** Non-SI units accepted for use with the SI ***/
    /*** ========================================= ***/

    constexpr auto a_u_length = 0.52917721067e-10; // ("a.u. of length: Bohr radius (bohr)", "a0", 0.52917721067e-10, "m", 2.3e-10)

    /*** ========================================= ***/
    /*** Atomic units (a.u.)                       ***/
    /*** ========================================= ***/

    constexpr auto a_u_action = 1.054571800e-34; // ("a.u. of action: h/2pi", "h_bar", 1.054571800e-34, "J s", 1.2e-8)


} // !phy_const
} // !uovie

#endif
