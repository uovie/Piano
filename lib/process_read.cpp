// file class member function "read()"
#include <iostream>
#include <sstream>
#include <cassert>
#include "../include/iostream.h"
#include "../include/phy_const.h"
#include "../include/atom_data.h"

using namespace uovie::Global;

void process::read() {

    constexpr double angstrom_to_bohr = 1e-10 / Phy_Const::a_u_length;

    std::cout << "Read Infomation from " << filename << std::endl;

    in >> sys.dimension;

    std::string line;
    
    int mol_ctr = 0;
    molecule tmp_mole;
    
    while (getline(in, line)) {
        std::string str;
        if (line.size() >= 5) {
            for (int a = 0; a < 5; a++)
                str.push_back(line[a]);
        }
        if (str == "*****") {
            int num_atom;
            in >> num_atom;
            sys.molecules.push_back(tmp_mole);
            for (int i = 0; i < num_atom; i++) {
                atom tmp_atom;
                sys.molecules[mol_ctr].atoms.push_back(tmp_atom);
                in >> sys.molecules[mol_ctr].atoms[i].symbol;
                for (int d = 0; d < sys.dimension; d++) {
                    double tmp_q, tmp_p, tmp_F;
                    sys.molecules[mol_ctr].atoms[i].q.push_back(tmp_q);
                    sys.molecules[mol_ctr].atoms[i].p.push_back(tmp_p);
                    sys.molecules[mol_ctr].atoms[i].F.push_back(tmp_F);
                    in >> sys.molecules[mol_ctr].atoms[i].q[d];
                    sys.molecules[mol_ctr].atoms[i].p[d] = 1;
                    sys.molecules[mol_ctr].atoms[i].F[d] = 0;
                }
            }
            // identify atomic numbers and assign atomic masses
            for (auto atom_it = sys.molecules[mol_ctr].atoms.begin(); atom_it != sys.molecules[mol_ctr].atoms.end(); atom_it++) {
                for (auto ele_it = Atom_Data::element.begin(); ele_it != Atom_Data::element.end(); ele_it++) {
                    if ((*atom_it).symbol == *ele_it) {
                        (*atom_it).atomic_number = static_cast<int>(ele_it - Atom_Data::element.begin()) + 1;
                        (*atom_it).m = 1; // Atom_Data::atomic_mass[(*atom_it).atomic_number - 1];
                        break;
                    }
                }
            }
        }
    }

    for (int mi = 0; mi < sys.molecules.size(); mi++) {
        for (int ai = 0; ai < sys.molecules[mi].atoms.size(); ai++) {
            std::cout << "\n" << sys.molecules[mi].atoms[ai].symbol;
            for (int di = 0; di < sys.dimension; di++) {
                std::cout << "\t" << sys.molecules[mi].atoms[ai].q[di];
            }
            std::cout << "\n";
        }
    }
}