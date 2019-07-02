/* Harmonic Bead-Spring Model */
#ifndef HARMONIC_BEAD_SPRING_H_
#define HARMONIC_BEAD_SPRING_H_

#include <vector>

namespace uoive {
namespace model {

    class hbs {

        int dimension;
        int num_part;
        double volume;
        double temperature;
        double pressure;
        std::vector<molecule> molecules;
    };


} // !model
} // !uovie


#endif // !HARMONIC_BEAD_SPRING_H_