/* Andersen Thermostat */

// uovie headers
#include "thermostat/at.h"

namespace uovie {
namespace thermostat {
namespace at {

    /*** ================================================== ***/
    /*** AT Procedure (Side) Class Member Functions         ***/
    /*** ================================================== ***/

    // initialization
    void at_base::initialize()
    {
        // resize arrays
        m.resize(dof);
        q.resize(dof);
        p.resize(dof);
        F.resize(dof);
        nrand.resize(dof);
        urand.resize(dof);

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

        // set normal-distributed random-number generators
        std::random_device rd;
        std::mt19937 tmp_nmte;
        std::normal_distribution<double> tmp_nd{ 0, 1 };
        for (int i = 0; i < dof; i++) {
            nmtes.push_back(tmp_nmte);
            nmtes[i].seed(rd());
            nds.push_back(tmp_nd);
        }

        // set uniform-distributed random-number generators
        std::mt19937 tmp_umte;
        std::uniform_real_distribution<double> tmp_ud{ 0, 1 };
        for (int i = 0; i < dof; i++) {
            umtes.push_back(tmp_umte);
            umtes[i].seed(rd());
            urds.push_back(tmp_ud);
        }

    }

    void at_base::upd_uni_rand_num()
    {
        for (int i = 0; i < dof; i++)
            urand(i) = urds[i](umtes[i]);
    }

    /*** ================================================== ***/
    /*** AT Procedure (Side) Class Member Functions         ***/
    /*** ================================================== ***/

    void at_side::implement_one_step()
    {
        upd_uni_rand_num();
        for (int i = 0; i < dof; i++) {
            if (urand(i) < cri) {
                nrand(i) = nds[i](nmtes[i]);
                p(i) = sqrt(1 / beta) * m(i) * nrand(i);
            }
        }
        
        calc_physic_force();
        p += F * Dt / 2;
        q += p * Dt / m;
        calc_physic_force();
        p += F * Dt / 2;

        upd_uni_rand_num();
        for (int i = 0; i < dof; i++) {
            if (urand(i) < cri) {
                nrand(i) = nds[i](nmtes[i]);
                p(i) = sqrt(1 / beta) * m(i) * nrand(i);
            }
        }
    }

    /*** ================================================== ***/
    /*** AT Procedure (Middle) Class Member Functions       ***/
    /*** ================================================== ***/

    void at_middle::implement_one_step()
    {
        calc_physic_force();
        p += F * Dt / 2;
        q += p * Dt / (2 * m);

        upd_uni_rand_num();
        for (int i = 0; i < dof; i++) {
            if (urand(i) < cri) {
                nrand(i) = nds[i](nmtes[i]);
                p(i) = sqrt(1 / beta) * m(i) * nrand(i);
            }
        }

        q += p * Dt / (2 * m);
        calc_physic_force();
        p += F * Dt / 2;
    }

} // !ld
} // !thermostat
} // !uovie