// Simulation Parameters
#ifndef SIMU_PARA_H_
#define SIMU_PARA_H_

namespace uovie {
namespace Global {

    class basic_simu_para {
    public:
        double run_time;                    // run time
        double time_step_size;              // step size
        int data_coll_peri;                 // data collection period
    };
    
}   // Global
}   // uovie
#endif