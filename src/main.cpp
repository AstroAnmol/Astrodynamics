#define _USE_MATH_DEFINES

#include <iostream>
#include <fstream>
#include <eigen-3.3.7/Eigen/Dense>
#include <cmath>
#include <queue>
#include <string>
#include <vector>
#include <random>
#include "orbit.h"
#include "cr3bp.h"
#include "indirect_method.h"
#include "GA.h"


int main(){
        srand(time(0));
        GA trial1;
        Eigen::VectorXd BestGeneT1;
        BestGeneT1=trial1.get_BestGene();;
        std::cout<<BestGeneT1.transpose()<<std::endl;
        std::cout<<trial1.get_BestFitness()<<std::endl;
        IndirectMethod Best(BestGeneT1);
        Best.save("Trial_1");
}