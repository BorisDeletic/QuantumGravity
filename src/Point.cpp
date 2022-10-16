//
// Created by Boris Deletic on 15/10/2022.
//

#include "Point.h"

using namespace std;

vector<Point> Point::getNeighbours(int n) {
    vector<Point> neighbours;

    for (int i = 0; i < point.size(); i++) {
        Point above(point);
        Point below(point);

        above.point[i] = point[i] + 1 < n ? point[i] + 1 : 0;
        below.point[i] = point[i] - 1 < 0 ? n - 1 : point[i] - 1;

        neighbours.push_back(above);
        neighbours.push_back(below);
    }

    return neighbours;
}


ostream &operator<<(ostream &o, const Point &l) {
    o << "(";
    for (int p : l.point) {
        o << p << ", ";
    }
    o << ")" << endl;
    return o;
}
