import pypolychord as pc
import time
from pypolychord.priors import UniformPrior,GaussianPrior
from phi4 import ScalarPhi4Action
import jax.numpy as jnp
import os
import anesthetic
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
# plt.style.use("prb")

M2 = -4.0
lam = 1.0
#3x3 lattice for now
L=8
nf=1
#change alpha to match figure
alpha=0.0
lattice_shape = (nf,L,L)

phi4_action = ScalarPhi4Action(M2=M2, lam=lam, lat=lattice_shape,alpha=alpha)

class lattice_prior(object):
    #Lattice prior is ~trivial for now but set up as a class
    def __init__(self):
        self.lattice_prior=GaussianPrior(mu=0,sigma=1)
    def __call__(self,theta):
        return self.lattice_prior(theta)

prior=lattice_prior()

#wrap the action up and for convenience save the magnetization as a derived param
def wrapped_like(theta):
    mag=phi4_action.mag(theta)
    return float(-(phi4_action(theta))),[mag]

nDims=jnp.ones(lattice_shape).flatten().shape[0]

settings = pc.settings.PolyChordSettings(nDims=nDims,nDerived=1)
settings.nlive=200
settings.read_resume = False
settings.write_resume = False
settings.base_dir="phi4"
settings.file_root="alpha_%s"%str(alpha)
settings.num_repeats = 100
settings.boost_posterior=settings.num_repeats
settings.nprior = settings.nlive * 10

pc.run_polychord(wrapped_like, prior=prior, nDims=nDims, nDerived=1, settings=settings)


def draw_magnetization(samples,filename):
    #
    print("saving png")
    mags=samples.to_numpy()[...,-4]
    f, a = plt.subplots()
    a.hist(mags,bins=30,weights=samples.get_weights())
    a.set_xlabel(r"$\bar{\phi}$")
    a.set_ylabel(r"Density")
    f.suptitle(r"$\phi^4$ theory, $V(\phi)=\frac{\mu}{2}\phi^2 + \lambda \phi^4 + \alpha\phi\,, \quad \lambda=%s, \mu=%s$"%(str(lam),str(M2)))
    f.savefig(filename+".png")

if rank==0:
    ns=anesthetic.read_chains(root=os.path.join(settings.base_dir,settings.file_root))
    draw_magnetization(ns,os.path.join(settings.base_dir,settings.file_root))

