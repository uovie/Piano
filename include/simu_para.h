// Simulation Parameters
#ifndef SIMU_PARA_H_
#define SIMU_PARA_H_

namespace uovie {
namespace Global {

    class basic_simu_para {
    public:
        double run_time;        // run time
        double step_size;       // time step size
        int data_coll_peri;     // data collection period
    };
    
} // !Global
} // !uovie
#endif // !SIMU_PARA_H_