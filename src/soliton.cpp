#include "soliton.h"

// constructor
Soliton::Soliton(double angle, double height, double vel_multiplier, Eigen::Vector3d pos, Eigen::Vector3d vel) {
    cone_angle = angle;
    cone_height = height;
    sol_vel_multiplier = vel_multiplier;
    debris_pos = pos;
    debris_vel = vel;

    soliton_velocity = sol_vel_multiplier * debris_vel;
    time_to_reach_cone_base = cone_height / soliton_velocity.norm(); // time to reach the apex of the cone
}

// set functions
void Soliton::set_params(double angle, double height, double vel_multiplier){
    cone_angle = angle;
    cone_height = height;
    sol_vel_multiplier = vel_multiplier;
    
    soliton_velocity = sol_vel_multiplier * debris_vel;
    time_to_reach_cone_base = cone_height / soliton_velocity.norm(); // time to reach the apex of the cone
}

void Soliton::set_debris_state(Eigen::Vector3d pos, Eigen::Vector3d vel) {
    debris_pos = pos;
    debris_vel = vel;
}

// get functions
Eigen::Vector3d Soliton::get_velocity() {
    return soliton_velocity;
}

double Soliton::get_time_to_reach_cone_base() {
    return time_to_reach_cone_base;
}
