// Simulation Parameters
#ifndef SIMU_PARA_H_
#define SIMU_PARA_H_


namespace uovie {
namespace Global {

    class basic_simu_para {
    public:
        basic_simu_para() = default;
        basic_simu_para(double rt, double ss, int dcp):
            run_time(rt), step_size(ss), data_coll_peri(dcp) { }

        double run_time;                    // run time
        double step_size;                   // step size
        int data_coll_peri;                 // data collection period
    };
    
}   // Global
}   // uovie

#endif