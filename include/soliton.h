#ifndef SOLITON_H
#define SOLITON_H

#include <eigen-5.0.0/Eigen/Dense>

class Soliton {
public:
    // constructor
    Soliton();
    Soliton(double angle, double height, double vel_multiplier, Eigen::Vector3d pos, Eigen::Vector3d vel);

    // set functions
    void set_params(double angle, double height, double vel_multiplier);
    void set_debris_state(Eigen::Vector3d pos, Eigen::Vector3d vel);

    // get functions
    Eigen::Vector3d get_velocity();
    double get_time_to_reach_cone_base();

private:
    // soliton parameters
    double cone_angle; // in radians
    double cone_height; // in km
    double sol_vel_multiplier;
    double time_to_reach_cone_base; // in seconds
    Eigen::Vector3d soliton_velocity; // in km/s

    Eigen::Vector3d debris_pos; // debris position when soliton is generated
    Eigen::Vector3d debris_vel; // debris velocity when soliton is generated
};

#endif // SOLITON_H