//
// Created by Boris Deletic on 15/10/2022.
//

#ifndef QUANTUMGRAVITY_POINT_H
#define QUANTUMGRAVITY_POINT_H

#include <vector>
#include <iostream>

using namespace std;

class Point {
public:
    Point(vector<int> point)
            : point(point) {}

    vector<Point> getNeighbours(int n); //n is lattice size

    vector<int> point;

    friend ostream& operator << (ostream &o, const Point &l);
};


#endif //QUANTUMGRAVITY_POINT_H
