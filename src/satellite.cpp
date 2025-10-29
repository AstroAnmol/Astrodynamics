#include "satellite.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>


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
    detection_freq = 10000; // 1000 Hz
    detections = 0;

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
    time = 0.0; // initial time

    // Body frame in ECI frame
    // construct rotation matrix from body frame to ECI frame
    Eigen::Vector3d sat_x_BF = velocity.normalized(); // x axis along velocity
    Eigen::Vector3d sat_z_BF = -position.normalized(); // z axis points to Earth center
    Eigen::Vector3d sat_y_BF = sat_z_BF.cross(sat_x_BF).normalized(); // y axis completes right-handed system

    R_BF_to_ECI.col(0) = sat_x_BF;
    R_BF_to_ECI.col(1) = sat_y_BF;
    R_BF_to_ECI.col(2) = sat_z_BF;

    // sensor positions in ECI frame
    sensor_1_ECI = position + R_BF_to_ECI * sensor_1_BF;
    sensor_2_ECI = position + R_BF_to_ECI * sensor_2_BF;
    sensor_3_ECI = position + R_BF_to_ECI * sensor_3_BF;
    sensor_4_ECI = position + R_BF_to_ECI * sensor_4_BF;

    // wake parameters
    plane_angle = 20.0 * M_PI / 180; // radians
    
    back_plane_normal_BF << -1, 0, 0;
    d_back_BF = sat_x_size/2;

    plane_normal_BF_1 << -std::sin(plane_angle), 0, -std::cos(plane_angle);
    d_1_BF = sat_x_size/2 *std::sin(plane_angle) - sat_z_size/2 *std::cos(plane_angle);

    plane_normal_BF_2 << -std::sin(plane_angle), 0, std::cos(plane_angle);
    d_2_BF = sat_x_size/2 *std::sin(plane_angle) - sat_z_size/2 *std::cos(plane_angle);

    plane_normal_BF_3 << -std::sin(plane_angle), std::cos(plane_angle), 0;
    d_3_BF = sat_x_size/2 *std::sin(plane_angle) - sat_y_size/2 *std::cos(plane_angle);

    plane_normal_BF_4 << -std::sin(plane_angle), -std::cos(plane_angle), 0;
    d_4_BF = sat_x_size/2 *std::sin(plane_angle) - sat_y_size/2 *std::cos(plane_angle);

    // soliton 
    soliton.set_params(45.0 * M_PI / 180, 10, 1.2); // 45 deg cone angle, 10 km height, 1.2x debris velocity
}

// Set parameters
void Satellite::set_soliton_state(Eigen::Vector3d pos, Eigen::Vector3d vel) {
    debris_position = pos;
    debris_velocity = vel;
    soliton.set_debris_state(debris_position, debris_velocity);
    time_to_reach_cone_base = soliton.get_time_to_reach_cone_base();
    std::cout<< "Soliton Velocity: " << soliton.get_velocity().transpose() << " km/s\n";
    sat_orbit.propagate_2BP(1/detection_freq, time_to_reach_cone_base, 0, "satellite_propagation");

    read_future_state("satellite_propagation");
    detections = 0;
}

// check if a point is in the satellite wake
bool Satellite::within_wake(Eigen::Vector3d pos) {
    // position vector from satellite center
    Eigen::Vector3d sat_to_pos = pos - position;
    
    bool behind_back_plane, behind_plane_1, behind_plane_2, behind_plane_3, behind_plane_4;
    behind_back_plane = (back_plane_normal_BF.dot(sat_to_pos) - d_back_BF) >= 0;
    behind_plane_1 = (plane_normal_BF_1.dot(sat_to_pos) - d_1_BF) >= 0;
    behind_plane_2 = (plane_normal_BF_2.dot(sat_to_pos) - d_2_BF) >= 0;
    behind_plane_3 = (plane_normal_BF_3.dot(sat_to_pos) - d_3_BF) >= 0;
    behind_plane_4 = (plane_normal_BF_4.dot(sat_to_pos) - d_4_BF) >= 0;

    return behind_back_plane && behind_plane_1 && behind_plane_2 && behind_plane_3 && behind_plane_4;
}

// check if the soliton is detected at current time
bool Satellite::detect_soliton() {
    bool detected = false;
    bool debris_within_wake = within_wake(debris_position);
    if (debris_within_wake) {
        // std::cout<< "Debris is within satellite wake at time " << time << " seconds.\n";
        // return false;
    }
    else{
        // check each sensor
        bool sensor_1_detected = soliton.within_cone(sensor_1_ECI) && soliton.within_spherical_range(sensor_1_ECI, time, detection_freq);
        bool sensor_2_detected = soliton.within_cone(sensor_2_ECI) && soliton.within_spherical_range(sensor_2_ECI, time, detection_freq);
        bool sensor_3_detected = soliton.within_cone(sensor_3_ECI) && soliton.within_spherical_range(sensor_3_ECI, time, detection_freq);
        bool sensor_4_detected = soliton.within_cone(sensor_4_ECI) && soliton.within_spherical_range(sensor_4_ECI, time, detection_freq);

        if (sensor_1_detected) {
            std::cout << "Sensor 1 detected soliton at time " << time << " seconds.\n";
            detections += 1;
        }
        if (sensor_2_detected) {
            std::cout << "Sensor 2 detected soliton at time " << time << " seconds.\n";
            detections += 1;
        }
        if (sensor_3_detected) {
            std::cout << "Sensor 3 detected soliton at time " << time << " seconds.\n";
            detections += 1;
        }
        if (sensor_4_detected) {
            std::cout << "Sensor 4 detected soliton at time " << time << " seconds.\n";
            detections += 1;
        }
        detected = sensor_1_detected || sensor_2_detected || sensor_3_detected || sensor_4_detected;
    }
    return detected;
}

// check if soliton is detected at any time in future_state
bool Satellite::detect_soliton_over_time() {
    
    for (int i = 0; i < future_state.rows(); ++i) {
        time = future_state(i, 0);
        Eigen::Vector3d pos= future_state.row(i).segment<3>(8);
        Eigen::Vector3d vel= future_state.row(i).segment<3>(11);
        position = pos;
        velocity = vel;
        // update sensor positions in ECI frame
        
        // construct rotation matrix from body frame to ECI frame
        Eigen::Vector3d sat_x_BF = velocity.normalized(); // x axis along velocity
        Eigen::Vector3d sat_z_BF = -position.normalized(); // z axis points to Earth center
        Eigen::Vector3d sat_y_BF = sat_z_BF.cross(sat_x_BF).normalized(); // y axis completes right-handed system

        R_BF_to_ECI.col(0) = sat_x_BF;
        R_BF_to_ECI.col(1) = sat_y_BF;
        R_BF_to_ECI.col(2) = sat_z_BF;

        // sensor positions in ECI frame
        sensor_1_ECI = position + R_BF_to_ECI * sensor_1_BF;
        sensor_2_ECI = position + R_BF_to_ECI * sensor_2_BF;
        sensor_3_ECI = position + R_BF_to_ECI * sensor_3_BF;
        sensor_4_ECI = position + R_BF_to_ECI * sensor_4_BF;

        bool local_detected = detect_soliton();
    }
    std::cout << "Total Detections: " << detections << "\n";
    if (detections>=2){return true;}
    else {return false;}
}

// Read a CSV file produced by orbit::propagate_2BP (header + numeric rows)
void Satellite::read_propagation_csv(const std::string &filename) {
    std::ifstream ifs(filename.c_str());
    if (!ifs.is_open()){
        std::cerr << "Satellite::read_propagation_csv: failed to open '" << filename << "'\n";
        return;
    }

    std::string line;
    // read header
    if (!std::getline(ifs, line)){
        std::cerr << "Satellite::read_propagation_csv: file empty: '" << filename << "'\n";
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

    future_state = out;
    return;
}

// Read satellite propagation file (tries "name_file.csv" then "name") and
// populate internal `future_state` matrix using read_propagation_csv.
void Satellite::read_future_state(std::string name) {
    std::string fn1 = "Results/" + name + "_file.csv";
    std::string fn2 = "Results/" + name;

    // try first filename
    read_propagation_csv(fn1);
    if (future_state.size() == 0) {
        // try second filename
        read_propagation_csv(fn2);
    }

    if (future_state.size() == 0) {
        std::cerr << "Satellite::read_future_state: failed to read '" << fn1 << "' or '" << fn2 << "'\n";
        return;
    }

    // Basic validation
    if (future_state.rows() < 1 || future_state.cols() < 14) {
        std::cerr << "Satellite::read_future_state: unexpected CSV layout (rows=" << future_state.rows() << ", cols=" << future_state.cols() << ")\n";
        // still keep future_state as read, but caller should handle it
    }

    return;
}