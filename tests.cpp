//
// Created by Boris Deletic on 16/10/2022.
//

#include <iostream>
#include "src/Point.h"
#include "src/LatticeFieldTheory.h"


void testNeighbours()
{
    Point testP({0,1,1,1});

    vector<Point> neighbours = testP.getNeighbours(10);
    for (Point neighbour : neighbours) {
        cout << neighbour;
    }
}

