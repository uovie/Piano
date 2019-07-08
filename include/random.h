#ifndef RANDOM_H_
#define RANDOM_H_

#include <random>

namespace uovie {
    // random-number engines
    std::mt19937 nhc_tmv_mte(27);
    std::mt19937 car_mom_mte(36);
}

#endif // !RANDOM_H_