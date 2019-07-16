// standard C++ headers
#include <iostream>
#include <sstream>
#include <random>
#include <cstdlib>
#include <cassert>

// Piano headers
#include "process.h"
#include "phy_const.h"
#include "atom_data.h"

using namespace uovie::Global;

void process::open(const std::string& filename) {
    // check the extension
    if (filename.rfind(".vie") == std::string::npos)
        throw std::invalid_argument("Only .vie file is a valid input file.");

    in.open(filename);

    if (in.fail()) {
        std::cout << "Can not open the file " << filename << '.' << std::endl;
        exit(EXIT_FAILURE);
    }

    for (auto it = filename.begin(); it != filename.end() - 4; it++)
        fn_no_ex += *it;
}

void process::read() {

    constexpr double k = uovie::phy_const::Boltzmann_const;

    std::cout << "Read Infomation from " << fn_no_ex + ".vie" << std::endl;
    in >> job;

    {
        std::string tmp_data;
        if (job == "ld") {
            in >> tmp_data;
            des.push_back(tmp_data);
        }
        else if (job == "at") {
            in >> tmp_data;
            des.push_back(tmp_data);
        }
        else if (job == "nhc") {
            for (int i = 0; i < 2; i++) {
                in >> tmp_data;
                des.push_back(tmp_data);
            }
        }
        else if (job == "pimd") {
            for (int i = 0; i < 2; i++) {
                in >> tmp_data;
                des.push_back(tmp_data);
            }
        }
    }

    in >> bsp.run_time >> bsp.step_size >> bsp.data_coll_peri;

    bsp.run_time *= 1e-15 / uovie::phy_const::a_u_time;
    bsp.step_size *= 1e-15 / uovie::phy_const::a_u_time;

    in >> sys.dimension >> sys.volume >> sys.temperature >> sys.pressure;

    sys.temperature /= phy_const::a_u_energy;

    in >> sys.model_type;

    if (sys.model_type == "HO") {
        double tmp_data;
        in >> tmp_data;
        sys.model_para.push_back(tmp_data);
    }
    else if (sys.model_type == "LJ") {
        double tmp_data;
        for (int i = 0; i < 2; i++) {
            in >> tmp_data;
            sys.model_para.push_back(tmp_data);
        }
        // unit conversion
        sys.model_para[0] *= phy_const::energy_K / phy_const::a_u_energy;
        sys.model_para[1] /= phy_const::a_u_length * 1e10;
    }else
        throw "unsupported model";

    const double& T = sys.temperature;
    
    int mol_ctr = 0; int& mi = mol_ctr;
    std::mt19937 mte(36);
    std::string line;
    molecule tmp_mole;
    
    while (getline(in, line)) {
        std::string str;
        if (line.size() >= 5) {
            for (int a = 0; a < 5; a++)
                str.push_back(line[a]);
        }
        if (str == "*****") {
            int natom;
            in >> natom;
            sys.molecules.push_back(tmp_mole);
            for (int ai = 0; ai < natom; ai++) {
                atom tmp_atom;
                sys.molecules[mi].atoms.push_back(tmp_atom);
                in >> sys.molecules[mi].atoms[ai].symbol;
                for (int di = 0; di < sys.dimension; di++) {
                    double tmp_q;
                    sys.molecules[mi].atoms[ai].q.push_back(tmp_q);
                    in >> sys.molecules[mi].atoms[ai].q[di];
                    sys.molecules[mi].atoms[ai].q[di] /= phy_const::a_u_length * 1e10;
                }
            }
            // identify atomic numbers and assign atomic masses, momentums, forces
            for (int ai = 0; ai < natom; ai++) {
                for (int ei = 0; ei < atom_data::element.size(); ei++) {
                    if (sys.molecules[mi].atoms[ai].symbol == atom_data::element[ei]) {
                        sys.molecules[mi].atoms[ai].atomic_number = ei + 1;
                        sys.molecules[mi].atoms[ai].m = atom_data::atomic_mass_kg[ei] / phy_const::a_u_mass;
                        break;
                    }
                }
                std::normal_distribution<double> ndrm{ 0, sqrt(k * T * sys.molecules[mi].atoms[ai].m) };
                for (int di = 0; di < sys.dimension; di++) {
                    sys.molecules[mi].atoms[ai].p.push_back(ndrm(mte));
                    sys.molecules[mi].atoms[ai].F.push_back(0);
                }
            }
            mol_ctr++;
        }
    }

    in.close(); // close files

    sys.num_part = 0;
    for (auto mi = 0; mi < sys.molecules.size(); mi++)
        sys.num_part += sys.molecules[mi].atoms.size();

}