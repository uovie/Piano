/* Langevin Dynamics */

// uovie headers
#include "thermostat/ld.h"
#include "model.h"

namespace uovie {
namespace thermostat {
namespace ld {

    /*** ===================================================== ***/
    /*** LD Procedure (Base) Class Member Functions            ***/
    /*** ===================================================== ***/

    // initialization
    void ld_base::initialize()
    {
        // resize arrays
        m.resize(dof);
        q.resize(dof);
        p.resize(dof);
        F.resize(dof);
        nrand.resize(dof);

        // initialize positions and masses
        int vi = 0;
        for (auto mi = 0; mi < sys.molecules.size(); mi++) {
            for (auto ai = 0; ai < sys.molecules[mi].atoms.size(); ai++) {
                for (auto di = 0; di < d; di++) {
                    m(vi) = sys.molecules[mi].atoms[ai].m;
                    q(vi) = sys.molecules[mi].atoms[ai].q[di];
                    vi++;
                }
            }
        }

        // initialize momenta
        std::mt19937 p_mte(36);
        for (auto ri = 0; ri < p.rows(); ri++) {
            std::normal_distribution<double> ndrm{ 0, sqrt(k * T * m(ri)) };
            p(ri) = ndrm(p_mte);
        }

        // set random-number generators
        std::random_device rd;
        std::mt19937 tmp_mte;
        std::normal_distribution<double> tmp_nd{ 0, 1 };
        for (int i = 0; i < dof; i++) {
            mtes.push_back(tmp_mte);
            mtes[i].seed(rd());
            nds.push_back(tmp_nd);
        }
    }

    void ld_base::upd_rand_num()
    {
        for (int i = 0; i < dof; i++)
            nrand(i) = nds[i](mtes[i]);
    }

    /*** ===================================================== ***/
    /*** LD Procedure (Side) Class Member Functions            ***/
    /*** ===================================================== ***/

    void ld_side::implement_one_step() {
        upd_rand_num();
        p = c1 * p + c2 * sqrt(1 / beta) * m.sqrt() * nrand;
        calc_physic_force();
        p += F * Dt / 2;
        q += p * Dt / m;
        calc_physic_force();
        p += F * Dt / 2;
        upd_rand_num();
        p = c1 * p + c2 * sqrt(1 / beta) * m.sqrt() * nrand;
    }

    /*** ===================================================== ***/
    /*** LD Procedure (middle) Class Member Functions          ***/
    /*** ===================================================== ***/

    void ld_middle::implement_one_step() {
        calc_physic_force();
        p += F * Dt / 2;
        q += p * Dt / (2 * m);
        upd_rand_num();
        p = c1 * p + c2 * sqrt(1 / beta) * m.sqrt() * nrand;
        q += p * Dt / (2 * m);
        calc_physic_force();
        p += F * Dt / 2;
    }

} // !ld
} // !thermostat
} // !uovie