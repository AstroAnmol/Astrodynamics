#ifndef IND_H
#define IND_H
#include <eigen-3.3.7/Eigen/Dense>
#include "orbit.h"

class IndirectMethod{
    public:
        /* IndirectMethod(); */
        IndirectMethod(Eigen::VectorXd param);
        double getFitness();
        Eigen::Vector3d Thrust(Eigen::VectorXd x);//x is R,V,m,Lambda_R, Lambda_V, lambda_m        
        Eigen::VectorXd differential(Eigen::VectorXd x);//x is R,V,m,Lambda_R, Lambda_V, lambda_m
        void propagate();//
        // to save the propagated trajectory
        void save(std::string name);
        // print some values
        void print(std::string value);
        // void set_GAparam(Eigen::VectorXd param);
    private:
        //Eigen::VectorXd GAparam; //t0, m0, V_inf_0_x, V_inf_0_y, Lambda_R0, tf
        double t0, tf;
        double m0;
        double lambda_m0=1;
        double mLEO;
        Eigen::Vector3d R0, V0, Lambda_R0, Lambda_V0;
        //Eigen::VectorXd x0;

        double v_inf_0;
        Eigen::Vector3d V_inf_0;

        double ve, vc, del_v; // Escape Velocity, Circular orbit Velocity

        //double T0, c, c_dash, epsilon;
        const double T0=0.02; // Thrust at 1 AU
        const double c=1;   // exhaust velocity
        const double c_dash=0.0955;
        const double epsilon=0.06;

        
        //orbit mars
        Eigen::Vector3d RM0, VM0; // Mars vector on 01/01/2000
        orbit Mars;
        Eigen::Vector3d RMtf, VMtf;

        //orbit earth
        Eigen::Vector3d RE0, VE0; // Earth vector on 01/01/2000
        orbit Earth;
        Eigen::Vector3d REt0, VEt0;

        double fitness;
};


#endif