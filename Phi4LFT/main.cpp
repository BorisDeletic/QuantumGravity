#include <iostream>
#include "src/Point.h"
#include "src/LatticeFieldTheory.h"
#include "tests.h"

using namespace std;

void meanMag() {
    int n = 4;

    int latticeSize = n*n*n*n;
    int equilSteps = latticeSize * 1000;
    int batches = 2000;

    double lambda = 1.5;
    double delta = 1.5;
    vector<double> kappas;
    for (double k = 0.08; k < 0.18; k+= 0.01) kappas.push_back(k);

    for (double k : kappas) {
        LatticeFieldTheory lattice(n, lambda, delta);
        lattice.runStepsMH(k, equilSteps);

        double meanMagnetisation = 0;

        for (int i = 0; i < batches; i++) {
            lattice.runStepsMH(k, latticeSize*2);
            meanMagnetisation += lattice.magnetisation();
        }

        meanMagnetisation = abs(meanMagnetisation / batches) / latticeSize;

        cout << k << "\t" << meanMagnetisation << "\n";
    }
}

int main() {
  //  testNeighbours({1,1,1,0});
    meanMag();

    return 0;
}
