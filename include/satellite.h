#ifndef SATELLITE_H
#define SATELLITE_H

#include <eigen-5.0.0/Eigen/Dense>
#include "orbit.h"
#include "soliton.h"

class Satellite {
public:
    // constructor
    Satellite();

    // set functions
    void set_soliton_state(Eigen::Vector3d pos, Eigen::Vector3d vol);

    // detection functions
    bool detect_soliton();
    bool detect_soliton_over_time();

private:
    // Body frame: x: velocity, y: right, z: down
    // size (in body frame)
    // dimensions of the satellite
    double sat_x_size;
    double sat_y_size;
    double sat_z_size;

    // boom length for sensors
    double boom_length;
    double detection_freq; // in Hz

    // sensor positions in body frame (mounted at corners of front yz face)
    Eigen::Vector3d sensor_1_BF, sensor_2_BF, sensor_3_BF, sensor_4_BF;
    Eigen::Vector3d corner_1_BF, corner_2_BF, corner_3_BF, corner_4_BF;

    // sensor positions in ECI frame
    Eigen::Vector3d sensor_1_ECI, sensor_2_ECI, sensor_3_ECI, sensor_4_ECI;


    Orbit sat_orbit; // satellite orbit

    // current state vectors
    double time; // in seconds
    Eigen::Vector3d position; // in km
    Eigen::Vector3d velocity; // in km/s

    // Future state vectors
    Eigen::ArrayXXd future_state; // time, position (x,y,z), velocity (vx,vy,vz)
    
    // Body frame in ECI frame
    Eigen::Matrix3d R_BF_to_ECI;

    // soliton
    Soliton soliton;
    double time_to_reach_cone_base;
    Eigen::Vector3d debris_position, debris_velocity;

    //csv read
    void read_propagation_csv(const std::string &filename);
    void read_future_state(std::string name);

    // wake parameters
    double plane_angle; // in radians

    // plane parameters (BF)
    Eigen::Vector3d back_plane_normal_BF, plane_normal_BF_1, plane_normal_BF_2, plane_normal_BF_3, plane_normal_BF_4;
    double d_back_BF, d_1_BF, d_2_BF, d_3_BF, d_4_BF;
    
    // Detection functions
    double detections;
    // check if a point is in the wake of satellite
    bool within_wake(Eigen::Vector3d pos);

};

#endif // SATELLITE_H