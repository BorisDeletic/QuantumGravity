import jax.numpy as jnp

class ScalarPhi4Action:
    def __init__(self, M2, lam,lat,alpha):
        self.M2 = M2
        self.lam = lam
        self.lat=lat
        self.alpha=alpha
        
    def __call__(self, theta):
        # alpha=theta[-1]
        # theta=theta[:-1]
        theta=theta.reshape(self.lat)
        # potential term
        action_density = self.M2/2 * theta ** 2 + self.lam * theta ** 4 + self.alpha/(self.lat[-1]*self.lat[-2]) * theta
        # kinetic term (discrete Laplacian)
        Nd = len(theta.shape) - 1
        dims = range(1,Nd+1)
        for mu in dims:
            action_density += 2 * theta ** 2
            action_density -= theta * jnp.roll(theta, -1, mu)
            action_density -= theta * jnp.roll(theta, 1, mu)
        return jnp.sum(action_density, axis=tuple(dims))
        
    def mag(self,theta):
        return jnp.sum(theta)/(self.lat[-1]*self.lat[-2])