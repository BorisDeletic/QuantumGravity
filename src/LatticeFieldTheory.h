//
// Created by Boris Deletic on 15/10/2022.
//

#ifndef QUANTUMGRAVITY_LATTICEFIELDTHEORY_H
#define QUANTUMGRAVITY_LATTICEFIELDTHEORY_H

#include <vector>
#include <random>
#include "Point.h"

using namespace std;

class LatticeFieldTheory {
public:
    LatticeFieldTheory(int n, double lambda, double delta);

    double potential(double x);
    double magnetisation();

    double fieldValue(Point& p) {return field[p.point[0]][p.point[1]][p.point[2]][p.point[3]]; };
    void setField(Point& p, double newField) { field[p.point[0]][p.point[1]][p.point[2]][p.point[3]] = newField; }

    double neighbourSum(Point& p);
    double scalarActionDifference(Point& p, double newField, double kappa);

    int stepMH(double kappa);
    int runStepsMH(double kappa, int steps);

    vector<vector<vector<vector<double> > > > field;

    int n;
    double lambda;
    double delta;

    random_device rd;  // Will be used to obtain a seed for the random number engine
    mt19937 gen; // Standard mersenne_twister_engine seeded with rd()
    uniform_int_distribution<> pointRNG; //random numbers in lattice
    uniform_real_distribution<double> uniformRNG;
};


#endif //QUANTUMGRAVITY_LATTICEFIELDTHEORY_H
