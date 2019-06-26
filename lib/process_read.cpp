// file class member function "read()"
#include <iostream>
#include <cassert>
#include "../include/iostream.h"
#include "../include/phy_const.h"
#include "../include/atom_data.h"

using namespace uovie::Global;

void process::read() {

    constexpr double angstrom_to_bohr = 1e-10 / Phy_Const::a_u_length;

    std::cout << "Read Infomation from " << filename << std::endl;

    std::string temp_data;
    int mol_ctr = 0;
    molecule mole;
    while (!in.eof()) {
        
        do {
            getline(in, temp_data);             // read redundant data
        } while (!in.eof() && temp_data[0] != '*');
        if (temp_data[0] != '*')
            break;

        sys.molecules.push_back(mole);
        int num_atom;
        in >> num_atom;
        
        atom atom_data;
        double temp_q[3];
        for (int i = 0; i < num_atom; i++) {
            in >> atom_data.symbol >> temp_q[0] >> temp_q[1] >> temp_q[2];

            atom_data.q[0] = temp_q[0] * angstrom_to_bohr;
            atom_data.q[1] = temp_q[1] * angstrom_to_bohr;
            atom_data.q[2] = temp_q[2] * angstrom_to_bohr;

            sys.molecules[mol_ctr].atoms.push_back(atom_data);
        }

        // identify atomic numbers and assign atomic masses
        for (auto atom_it = sys.molecules[mol_ctr].atoms.begin(); atom_it != sys.molecules[mol_ctr].atoms.end(); atom_it++) {
            for (auto ele_it = Atom_Data::element.begin(); ele_it != Atom_Data::element.end(); ele_it++) {
                if ((*atom_it).symbol == *ele_it) {
                    (*atom_it).atomic_number = (int)(ele_it - Atom_Data::element.begin()) + 1;
                    (*atom_it).m = 1; // Atom_Data::atomic_mass[(*atom_it).atomic_number - 1];
                    break;
                }
            }
        }

        mol_ctr++;
    }
}