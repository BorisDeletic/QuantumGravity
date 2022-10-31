//
// Created by Boris Deletic on 15/10/2022.
//

#include "LatticeFieldTheory.h"


LatticeFieldTheory::LatticeFieldTheory(int n, double lambda, double delta) :
        n(n),
        lambda(lambda),
        delta(delta),
        gen(rd()),
        pointRNG(0, n - 1),
        uniformRNG(0,1)
{
    field = vector<vector<vector<vector<double> > > >(n,
                vector<vector<vector<double> > >(n,
                     vector<vector<double> >(n,
                           vector<double>(n, 0)
                           )
                    )
            );
}

double LatticeFieldTheory::potential(double x) {
    return lambda*(x*x-1) * (x*x-1) + x*x;
}


double LatticeFieldTheory::neighbourSum(Point& p) {
    // compute the sum of all neighbour field values at a point

    vector<Point> neighbours = p.getNeighbours(n);

    double fieldSum = 0;
    for (Point p : neighbours) {
        fieldSum += fieldValue(p);
    }

    return fieldSum;
}

double LatticeFieldTheory::scalarActionDifference(Point &p, double newField, double kappa) {
    // compute change in the action for new field value at a point
    double oldField = fieldValue(p);
    double potentialChange = potential(oldField) - potential(newField);

    double actionChange = 2 * kappa * neighbourSum(p) * (oldField - newField) - potentialChange;

    return actionChange;
}

int LatticeFieldTheory::stepMH(double kappa) {
    vector<int> rands;
    rands.push_back(pointRNG(gen));
    rands.push_back(pointRNG(gen));
    rands.push_back(pointRNG(gen));
    rands.push_back(pointRNG(gen));

    Point p = Point(rands);

    uniform_real_distribution<double> distr(-delta, delta);
    double rngField = distr(gen);

    double newField = fieldValue(p) + rngField;
    double actionChange = scalarActionDifference(p, newField, kappa);

    double r = uniformRNG(gen); // random in (0,1)

    if ((actionChange < 0) || r < exp(-actionChange)) {
        setField(p, newField);
        return 1;
    }
    return 0;
}

int LatticeFieldTheory::runStepsMH(double kappa, int steps) {
    int accepted = 0;

    for (int i = 0; i < steps; i++) {
        accepted += stepMH(kappa);
    }
    return accepted;
}

double LatticeFieldTheory::magnetisation() {
    // return mean field value in lattice
    double mag = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                for (int l = 0; l < n; l++) {
                    mag += field[i][j][k][l];
                }
            }
        }
    }

    return mag;
}

