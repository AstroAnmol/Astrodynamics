#include "satellite.h"

// constructor
Satellite::Satellite() {
    // Define default satellite parameters

    // Body frame: x: velocity, y: right, z: down
    // size (in body frame)
    sat_x_size = 3.0*0.0001; // 30 centimeters (3U)
    sat_y_size = 2.0*0.0001;   // 20 centimeters (2U)
    sat_z_size = 2.0*0.0001;   // 20 centimeters (2U)

    // boom length for sensors
    boom_length = 0.001; // 1 meter

    // define a satellite body frame (x: velocity, y: right, z: down)
    corner_1_BF = Eigen::Vector3d(sat_x_size/2, sat_y_size/2, sat_z_size/2);
    corner_2_BF = Eigen::Vector3d(sat_x_size/2, -sat_y_size/2, sat_z_size/2);
    corner_3_BF = Eigen::Vector3d(sat_x_size/2, -sat_y_size/2, -sat_z_size/2);
    corner_4_BF = Eigen::Vector3d(sat_x_size/2, sat_y_size/2, -sat_z_size/2);

    sensor_1_BF = corner_1_BF + boom_length*Eigen::Vector3d(1, 1, 1).normalized();
    sensor_2_BF = corner_2_BF + boom_length*Eigen::Vector3d(1, -1, 1).normalized();
    sensor_3_BF = corner_3_BF + boom_length*Eigen::Vector3d(1, -1, -1).normalized();
    sensor_4_BF = corner_4_BF + boom_length*Eigen::Vector3d(1, 1, -1).normalized();

    // satellite orbit
    Eigen::VectorXd OE_sat(6);
    double a, e, i, omega, Omega, theta;
    a =         750 + 6371; // km;
    e =         0.063;
    i =         45;
    omega =     0;
    Omega =     0;
    theta =     0;

    OE_sat << a, e, i, omega, Omega, theta;
    sat_orbit.set_OE(OE_sat);

    // current state vectors
    Eigen::VectorXd sat_state = sat_orbit.get_cartesian();
    position = sat_state.segment(0,3); // in km
    velocity = sat_state.segment(3,3); // in km/s

    // Body frame in ECI frame
    // construct rotation matrix from body frame to ECI frame
    Eigen::Vector3d sat_x_BF = velocity.normalized(); // x axis along velocity
    Eigen::Vector3d sat_z_BF = -position.normalized(); // z axis points to Earth center
    Eigen::Vector3d sat_y_BF = sat_z_BF.cross(sat_x_BF).normalized(); // y axis completes right-handed system

    R_BF_to_ECI.col(0) = sat_x_BF;
    R_BF_to_ECI.col(1) = sat_y_BF;
    R_BF_to_ECI.col(2) = sat_z_BF;
}