#ifndef SATELLITE_H
#define SATELLITE_H

#include <eigen-5.0.0/Eigen/Dense>
#include "orbit.h"

class Satellite {
public:
    // constructor
    Satellite();

    // set functions

private:
    // Body frame: x: velocity, y: right, z: down
    // size (in body frame)
    // dimensions of the satellite
    double sat_x_size;
    double sat_y_size;
    double sat_z_size;

    // boom length for sensors
    double boom_length;

    // sensor positions in body frame (mounted at corners of front yz face)
    Eigen::Vector3d sensor_1_BF, sensor_2_BF, sensor_3_BF, sensor_4_BF;
    Eigen::Vector3d corner_1_BF, corner_2_BF, corner_3_BF, corner_4_BF;


    orbit sat_orbit; // satellite orbit

    // current state vectors
    Eigen::Vector3d position; // in km
    Eigen::Vector3d velocity; // in km/s

    // Body frame in ECI frame

    Eigen::Matrix3d R_BF_to_ECI;

};

#endif // SATELLITE_H