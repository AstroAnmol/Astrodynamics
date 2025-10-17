#include "debris.h"
#include "orbit.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

// Read a CSV file produced by orbit::propagate_2BP (header + numeric rows)
void Debris::read_propagation_csv(const std::string &filename) {
    std::ifstream ifs(filename.c_str());
    if (!ifs.is_open()){
        std::cerr << "Debris::read_propagation_csv: failed to open '" << filename << "'\n";
        return;
    }

    std::string line;
    // read header
    if (!std::getline(ifs, line)){
        std::cerr << "Debris::read_propagation_csv: file empty: '" << filename << "'\n";
        return;
    }

    std::vector<std::vector<double>> rows;

    while (std::getline(ifs, line)){
        if (line.size() == 0) continue;
        std::vector<double> values;
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ',')){
            // trim
            size_t start = cell.find_first_not_of(" \t\r\n");
            size_t end = cell.find_last_not_of(" \t\r\n");
            if (start == std::string::npos) { cell = ""; }
            else cell = cell.substr(start, end - start + 1);

            if (cell.empty()) { values.push_back(0.0); continue; }

            try {
                double v = std::stod(cell);
                values.push_back(v);
            } catch (...) {
                // filter non-numeric characters
                std::string filtered;
                for (char c: cell) if ((c>='0' && c<='9') || c=='-' || c=='+' || c=='.' || c=='e' || c=='E') filtered.push_back(c);
                if (!filtered.empty()){
                    try { values.push_back(std::stod(filtered)); }
                    catch(...) { values.push_back(0.0); }
                } else {
                    values.push_back(0.0);
                }
            }
        }
        if (!values.empty()) rows.push_back(values);
    }

    if (rows.empty()) {
        std::cerr << "Debris::read_propagation_csv: no data rows found in '" << filename << "'\n";
        return;
    }

    size_t cols = 0;
    for (auto &r: rows) if (r.size() > cols) cols = r.size();

    Eigen::ArrayXXd out(rows.size(), cols);
    out.setZero();
    for (size_t i=0;i<rows.size();++i){
        for (size_t j=0;j<rows[i].size();++j) out(i,j) = rows[i][j];
    }

    sat_orbit = out;
    return;
}

// Read satellite propagation file (tries "name_file.csv" then "name") and
// populate internal `sat_orbit` matrix using read_propagation_csv.
void Debris::read_sat_orbit(std::string name) {
    std::string fn1 = name + "_file.csv";
    std::string fn2 = name;

    // try first filename
    read_propagation_csv(fn1);
    if (sat_orbit.size() == 0) {
        // try second filename
        read_propagation_csv(fn2);
    }

    if (sat_orbit.size() == 0) {
        std::cerr << "Debris::read_sat_orbit: failed to read '" << fn1 << "' or '" << fn2 << "'\n";
        return;
    }

    // Basic validation
    if (sat_orbit.rows() < 1 || sat_orbit.cols() < 14) {
        std::cerr << "Debris::read_sat_orbit: unexpected CSV layout (rows=" << sat_orbit.rows() << ", cols=" << sat_orbit.cols() << ")\n";
        // still keep sat_orbit as read, but caller should handle it
    }

    return;
}


// set state
void Debris::set_state(Eigen::Vector3d pos, Eigen::Vector3d vel) {
    position = pos;
    velocity = vel;

    soliton_vel = vel * 1.2; // assuming soliton moves with 1.2 times the debris initially
}

// set soliton cone parameters
void Debris::set_soliton_cone(double angle, double height) {
    cone_angle = angle;
    cone_height = height;

    time_to_reach_cone_base = cone_height / soliton_vel.norm(); // time to reach the apex of the cone
}



void Debris::check_detection() {
    // Function to check if the satellite is within the soliton cone
    // This function would use the satellite's position from its orbit
    // and compare it with the debris position and soliton parameters.
    // Implementation would depend on the specific detection criteria.
    bool detected_global = false;
    for (int i = 0; i < sat_orbit.rows(); ++i) {
        // Time (sec),Semi-Major Axis (km),Eccentricity,Inclination (deg),RAAN (deg),Argument of Periapsis (deg),True Anomaly (deg),Orbital Energy (km^2/sec^2),Radius_1 (km),Radius_2 (km),Radius_3 (km),Velocity_1 (km/s),Velocity_2 (km/s),Velocity_3 (km/s),Acceleration_1 (km/s^2),Acceleration_2 (km/s^2),Acceleration_3 (km/s^2),Angular_Momemntum_1 (km^2/s),Angular_Momemntum_2 (km^2/s),Angular_Momemntum_3 (km^2/s)
        Eigen::Vector3d sat_pos = sat_orbit.row(i).segment<3>(8);
        Eigen::Vector3d sat_vel = sat_orbit.row(i).segment<3>(11);
        double time = sat_orbit(i,0);
        
        // define 4 sensor locations (in km from satellite center) in body frame
        // define a satellite body frame (x: velocity, y: right, z: down)
        Eigen::Vector3d sensor_1_BF, sensor_2_BF, sensor_3_BF, sensor_4_BF;
        Eigen::Vector3d corner_1_BF, corner_2_BF, corner_3_BF, corner_4_BF;

        double boom_length = 0.001; // 1 meters boom length
        double sat_x_size = 3.0*0.0001; // 30 centimeters
        double sat_yz_size = 2.0*0.0001; // 20 centimeters
        corner_1_BF = Eigen::Vector3d(sat_x_size/2, sat_yz_size/2, sat_yz_size/2);
        corner_2_BF = Eigen::Vector3d(sat_x_size/2, -sat_yz_size/2, sat_yz_size/2);
        corner_3_BF = Eigen::Vector3d(sat_x_size/2, -sat_yz_size/2, -sat_yz_size/2);
        corner_4_BF = Eigen::Vector3d(sat_x_size/2, sat_yz_size/2, -sat_yz_size/2);

        sensor_1_BF = corner_1_BF + boom_length*Eigen::Vector3d(1, 1, 1).normalized();
        sensor_2_BF = corner_2_BF + boom_length*Eigen::Vector3d(1, -1, 1).normalized();
        sensor_3_BF = corner_3_BF + boom_length*Eigen::Vector3d(1, -1, -1).normalized();
        sensor_4_BF = corner_4_BF + boom_length*Eigen::Vector3d(1, 1, -1).normalized();

        // construct rotation matrix from body frame to ECI frame
        Eigen::Vector3d sat_x_BF = sat_vel.normalized(); // x axis along velocity
        Eigen::Vector3d sat_z_BF = -sat_pos.normalized(); // z axis points to Earth center
        Eigen::Vector3d sat_y_BF = sat_z_BF.cross(sat_x_BF).normalized(); // y axis completes right-handed system

        Eigen::Matrix3d R_BF_to_ECI;
        R_BF_to_ECI.col(0) = sat_x_BF;
        R_BF_to_ECI.col(1) = sat_y_BF;
        R_BF_to_ECI.col(2) = sat_z_BF;

        // sensor positions in ECI frame
        Eigen::Vector3d sensor_1_ECI = sat_pos + R_BF_to_ECI * sensor_1_BF;
        Eigen::Vector3d sensor_2_ECI = sat_pos + R_BF_to_ECI * sensor_2_BF;
        Eigen::Vector3d sensor_3_ECI = sat_pos + R_BF_to_ECI * sensor_3_BF;
        Eigen::Vector3d sensor_4_ECI = sat_pos + R_BF_to_ECI * sensor_4_BF;

        // check if any sensor is within the cone and within spherical detection range
        static bool detected = false; // to avoid multiple prints
        if (within_cone(sensor_1_ECI) && within_spherical_range(sensor_1_ECI, time)){
            detected = true;
            detected_global = true;
            // print detection event
            std::cout << "Detection at time " << time << " sec by sensor 1\n";
        } 
        if (within_cone(sensor_2_ECI) && within_spherical_range(sensor_2_ECI, time)){
            detected = true;
            detected_global = true;
            std::cout << "Detection at time " << time << " sec by sensor 2\n";
        }
        if (within_cone(sensor_3_ECI) && within_spherical_range(sensor_3_ECI, time)){
            detected = true;
            detected_global = true;
            std::cout << "Detection at time " << time << " sec by sensor 3\n";
        }
        if (within_cone(sensor_4_ECI) && within_spherical_range(sensor_4_ECI, time)){
            detected = true;
            detected_global = true;
            std::cout << "Detection at time " << time << " sec by sensor 4\n";
        }
    }
    // if never detected, you can print no detection
    if (!detected_global){
        std::cout << "No detection";
    }
}

bool Debris::within_cone(Eigen::Vector3d pos) {
    // debris position is the cone apex
    Eigen::Vector3d cone_apex = position;

    // position vector from apex
    Eigen::Vector3d debris_to_pos = pos - position;

    bool within_height = (debris_to_pos.norm() <= cone_height);
    bool within_angle = (debris_to_pos.dot(soliton_vel) / (debris_to_pos.norm() * soliton_vel.norm()) >= cos(cone_angle));
    bool same_side = (debris_to_pos.dot(soliton_vel) >= 0);

    return (within_height && within_angle && same_side);
}

bool Debris::within_spherical_range(Eigen::Vector3d pos, double time) {
    Eigen::Vector3d debris_to_pos = pos - position;
    double range = time * soliton_vel.norm();
    double del_range = 0.001 * soliton_vel.norm();
    return (debris_to_pos.norm() >= range - del_range && debris_to_pos.norm() <= range + del_range);
}

// get time to reach cone base
double Debris::get_time_to_reach_cone_base() {
    return time_to_reach_cone_base;
}

// get soliton velocity magnitude
double Debris::get_soliton_vel() {
    return soliton_vel.norm();
}