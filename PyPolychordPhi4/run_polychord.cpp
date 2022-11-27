#include "interfaces.hpp"
#include "interfaces.h"

#include <random>
#include <iomanip>
#include <iostream>
#include <unistd.h>
#include <cstring>
#include <memory>
#include <math.h>
#include <functional>

// #include "Integrator.hpp"
#include "Integrands.hpp"

using namespace std;
using namespace SamplingForStripper;

// Wrapper around random number generator for convience.
// Mimics the structure in STRIPPER.
struct rand
{
    static ranlux24 gen;
    static uniform_real_distribution<> dist;
    static void setseed(unsigned i) { gen.seed(i); };
    static double getrand() { return dist(gen); };
};

ranlux24 rand::gen;
uniform_real_distribution<> rand::dist(0., 1.);

shared_ptr<Function> fptr = nullptr;

double loglikelihood(double theta[], int nDims, double phi[], int nDerived)
{
    return fptr->logl({theta, theta + nDims});
}

// Prior function
//
// Either write your prior code directly into this function, or call an
// external library from it. This should transform a coordinate in the unit hypercube
// stored in cube (of size nDims) to a coordinate in the physical system stored in theta
//
// This function is called from likelihoods/fortran_cpp_wrapper.f90
// If you would like to adjust the signature of this call, then you should adjust it there,
// as well as in likelihoods/my_cpp_likelihood.hpp
//
void prior(double cube[], double theta[], int nDims)
{
    //============================================================
    // insert prior code here
    //
    //
    //============================================================
    for (int i = 0; i < nDims; i++)
        theta[i] = cube[i];
}

// Dumper function
//
// This function gives you runtime access to variables, every time the live
// points are compressed by a factor settings.compression_factor.
//
// To use the arrays, subscript by following this example:
//
//    for (int i_dead=0;i_dead<ndead;i_dead++)
//    {
//        for (int j_par=0;j_par<npars;j_par++)
//            std::cout << dead[npars*i_dead+j_par] << " ";
//        std::cout << std::endl;
//    }
//
// in the live and dead arrays, the rows contain the physical and derived
// parameters for each point, followed by the birth contour, then the
// loglikelihood contour
//
// logweights are posterior weights
//
void dumper(int ndead, int nlive, int npars, double *live, double *dead, double *logweights, double logZ, double logZerr)
{
    cout << "Current integral estimate  = " << exp(logZ) << " +/- " << exp(logZ + logZerr) - exp(logZ) << endl;
    // cout << "Number of evaluations = "<< nlive << endl;
}

int main(int argc, char *argv[])
{

    /****************************************************************************
     *                                                                          *
     * parse command line                                                       *
     *                                                                          *
     ****************************************************************************/

    if (argc == 1)
    {
        const string space(strlen(argv[0]), ' ');
        cout << "usage: " << argv[0] << " [ -i intergrator] [ -n points ] [ -d dim ]\n"
             << "       " << space << " [ -f function ] [ -s seed ] [ -o filename ]\n"
             << "       " << space << " [ -l filename ] [ -p points ]\n"
             << "       This routine performs a Monte Carlo integration of a selected predefiend function\n"
             << "       modelling a cross section. This is to validate new integrators and to measure performance.\n"
             << "       Optional parameters are:\n\n"
             << "       -i integrator: 0 - Vanilla MC (default)\n"
             << "                      1 - MySampler\n"
             << "       -f function  : 0 - d-dim. gaussian with random mean & width\n\n"
             << "                      1 - same as 0 but with random theta peaks\n"
             //    << "       -p number of point to integrate during optimization (default 0)\n\n"
             << "       -n number of live points (default 10^2)\n\n"
             //    << "       -x Max number of dead point, early termination (default = -1)\n\n"
             << "       -d integration space dimensions (default 1)\n\n"
             << "       -o filename  : chains file location\n\n"
             //    << "       -l filename  : load optimization dump from specified file\n\n"
             << "       -s random number genertor seed\n"
             << "       -p boost prior samples by factor (default -1)\n";
        return 0;
    }

    unsigned integrator = 0;
    unsigned function = 0;
    unsigned long long nlive = 100;
    unsigned long long prior_factor = 10;
    // unsigned long long maxndead = -1;
    // unsigned long long nevents_opt = 0;
    unsigned dimensions = 1;
    unsigned seed = 42;
    string filename = "chains";

    int c;
    opterr = 0;
    while ((c = getopt(argc, argv, "i:f:n:d:o:s:p:")) != -1)
        switch (c)
        {
        case 'i':
            integrator = atoi(optarg);
            break;
        case 'f':
            function = atoi(optarg);
            break;
        case 'n':
            nlive = atoi(optarg);
            break;
        case 'p':
            prior_factor = atoi(optarg);
            break;
        case 'd':
            dimensions = atoi(optarg);
            break;
        // case 'x':
        //     maxndead = atoi(optarg);
        //     break;
        case 's':
            seed = atoi(optarg);
            break;
        case 'o':
            filename = optarg;
            break;
        case '?':
            if (optopt == 'i' || optopt == 'f' || optopt == 'n' ||
                optopt == 'd' || optopt == 'o' || optopt == 's' || optopt == 'p') //|| optopt == 'x')
                cerr << "Missing argument for option -"
                     << static_cast<char>(optopt) << endl;
            else
                cerr << "Unknown option -"
                     << static_cast<char>(optopt) << endl;
        default:
            return 1;
        }

    cout << "# Configuration finished:" << endl
         << "# Integrator : " << integrator << endl
         << "# Function   : " << function << endl
         << "# Live points     : " << nlive << endl
         << "# Dimensions : " << dimensions << endl
         << "# Seed       : " << seed << endl
         << "# Number prior samples       : " << prior_factor*nlive << endl;

    // << "# Max Dead       : " << maxndead << endl;

    /****************************************************************************
     *                                                                          *
     * Construct integrand                                                      *
     *                                                                          *
     ****************************************************************************/

    switch (function)
    {
    case 0:
        fptr = make_shared<RandomGaussian>(dimensions);
        break;
    case 1:
    {
        auto newfunc = make_shared<RandomGaussianPeaked>(dimensions);
        vector<double> pos_min;
        vector<double> pos_max;

        double delta = pow(1e-3, 1. / static_cast<double>(dimensions));
        for (unsigned i = 0; i < dimensions; i++)
        {
            pos_min.push_back(0.1 - delta / 2.);
            pos_max.push_back(0.1 + delta / 2.);
        }
        newfunc->addPeak(1., pos_min, pos_max);
        fptr = newfunc;
    };
    break;
    default:
        cerr << "Function not implemented";
        return -1;
    }

    /****************************************************************************
     *                                                                          *
     * Construct integrator                                                     *
     *                                                                          *
     ****************************************************************************/

    // initialize random number generator
    // rand::setseed(seed);

    // shared_ptr<Integrator> iptr = nullptr;

    // switch (integrator)
    // {
    //     case 0: iptr = make_shared<Integrator>(&rand::getrand,dimensions);
    //             break;
    //     case 1:
    //     {
    //         unsigned par = dimensions;
    //         iptr = make_shared<MySampler>(&rand::getrand,par);
    //         break;
    //     }
    //     default:
    //     cerr << "Integrator not implemented"; return -1;
    // }

    /*
    Polychord specific settings, copied from $POLYCHORD_DIR/src/drivers/polychord_CC.cpp
    */

    int nDerived;
    nDerived = 0;

    Settings settings(dimensions, nDerived);

    settings.nlive = nlive;
    settings.num_repeats = settings.nDims * 5;
    settings.do_clustering = false;

    settings.precision_criterion = 1e-3;
    settings.logzero = -1e30;

    settings.base_dir = filename;
    settings.file_root = "samples";

    settings.write_resume = false;
    settings.read_resume = false;
    settings.write_live = true;
    settings.write_dead = true;
    settings.write_stats = true;

    // settings.max_ndead=maxndead;
    settings.seed = seed;
    settings.equals = false;
    settings.posteriors = true;
    settings.cluster_posteriors = false;
    settings.maximise = false;

    settings.feedback = 1;
    settings.compression_factor = 0.36787944117144233;
    settings.synchronous = true;

    settings.nprior = prior_factor * settings.nlive;
    
    settings.boost_posterior = dimensions * settings.num_repeats;

    run_polychord(loglikelihood, prior, dumper, settings);
    cout << "Analytic integral value  = " << fptr->integral() << endl;
}
