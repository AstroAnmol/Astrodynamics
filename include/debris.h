#ifndef DEBRIS_H
#define DEBRIS_H

#include <eigen-5.0.0/Eigen/Dense>
#include "orbit.h"

class Debris {
public:
    // constructor

    // set functions
    void set_state(Eigen::Vector3d pos, Eigen::Vector3d vel);
    void set_soliton_cone(double angle, double height);

    void check_detection();


    // get functions
    double get_time_to_reach_cone_base();
    double get_soliton_vel();

    // read satellite orbit
    void read_sat_orbit(std::string name);
    // Read a CSV produced by orbit::propagate_2BP (or a compatible numeric CSV).
    void read_propagation_csv(const std::string &filename);

private:
    // mass and size
    double mass; // in kg
    double size; // in meters

    // state vectors
    Eigen::Vector3d position; // in km
    Eigen::Vector3d velocity; // in km/s

    // soliton parameters
    Eigen::Vector3d soliton_vel;

    // cone parameters for soliton
    double cone_angle; // in radians
    double cone_height; // in km
    double time_to_reach_cone_base; // in seconds

    // detection satellite orbit
    Eigen::ArrayXXd sat_orbit; // time, position (x,y,z), velocity (vx,vy,vz)

    // check within cone
    bool within_cone(Eigen::Vector3d pos);

    // check within spherical detection range
    bool within_spherical_range(Eigen::Vector3d pos, double time);

};

#endif // DEBRIS_H