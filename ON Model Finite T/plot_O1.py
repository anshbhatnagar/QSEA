from cProfile import label
import csv
from math import factorial
import numpy as np
import matplotlib.pyplot as plt
from cosmoTransitions import generic_potential
import json

def pot(phi, mPhi2, alphaPhi, lPhi):

    freePhi4 = (mPhi2/2)*phi**2+(alphaPhi/6)*phi**3+(lPhi/24)*phi**4

    return freePhi4

class model1(generic_potential.generic_potential):
    """
    A sample model which makes use of the *generic_potential* class.

    This model doesn't have any physical significance. Instead, it is chosen
    to highlight some of the features of the *generic_potential* class.
    It consists of two scalar fields labeled *phi1* and *phi2*, plus a mixing
    term and an extra boson whose mass depends on both fields.
    It has low-temperature, mid-temperature, and high-temperature phases, all
    of which are found from the *getPhases()* function.
    """
    def init(self, vev, gamma):
        """
          m1 - tree-level mass of first singlet when mu = 0.
          m2 - tree-level mass of second singlet when mu = 0.
          mu - mass coefficient for the mixing term.
          Y1 - Coupling of the extra boson to the two scalars individually
          Y2 - Coupling to the two scalars together: m^2 = Y2*s1*s2
          n - degrees of freedom of the boson that is coupling.
        """
        # The init method is called by the generic_potential class, after it
        # already does some of its own initialization in the default __init__()
        # method. This is necessary for all subclasses to implement.

        # This first line is absolutely essential in all subclasses.
        # It specifies the number of field-dimensions in the theory.
        self.Ndim = 1

        # This next block sets all of the parameters that go into the potential
        # and the masses. This will obviously need to be changed for different
        # models.
        
        self.v = vev
        sf = gamma
        self.mPhi2 = -2*sf**4/self.v**2
        self.lPhi = 12*(sf/self.v)**4

        self.alphaPhi = 0

        self.Tmax = 10.


        # self.renormScaleSq is the renormalization scale used in the
        # Coleman-Weinberg potential.
        self.renormScaleSq = (4)**2

    def forbidPhaseCrit(self, X):
        """
        forbidPhaseCrit is useful to set if there is, for example, a Z2 symmetry
        in the theory and you don't want to double-count all of the phases. In
        this case, we're throwing away all phases whose zeroth (since python
        starts arrays at 0) field component of the vev goes below -5. Note that
        we don't want to set this to just going below zero, since we are
        interested in phases with vevs exactly at 0, and floating point numbers
        will never be accurate enough to ensure that these aren't slightly
        negative.
        """
        return (np.array([X])[...,0] < -5.0).any()

    def V0(self, X):
        """
        This method defines the tree-level potential. It should generally be
        subclassed. (You could also subclass Vtot() directly, and put in all of
        quantum corrections yourself).
        """
        # X is the input field array. It is helpful to ensure that it is a
        # numpy array before splitting it into its components.
        X = np.asanyarray(X)
        # x and y are the two fields that make up the input. The array should
        # always be defined such that the very last axis contains the different
        # fields, hence the ellipses.
        # (For example, X can be an array of N two dimensional points and have
        # shape (N,2), but it should NOT be a series of two arrays of length N
        # and have shape (2,N).)
        phi = X[...,0]

        return pot(phi, self.mPhi2, self.alphaPhi, self.lPhi)

    def boson_massSq(self, X, T):
        X = np.array(X)
        phi = X[...,0]

        # We need to define the field-dependent boson masses. This is obviously
        # model-dependent.
        # Note that these can also include temperature-dependent corrections.
        mPhiSq = self.mPhi2+(self.lPhi/2)*phi**2+self.alphaPhi*phi
        M = np.array([mPhiSq])

        # At this point, we have an array of boson masses, but each entry might
        # be an array itself. This happens if the input X is an array of points.
        # The generic_potential class requires that the output of this function
        # have the different masses lie along the last axis, just like the
        # different fields lie along the last axis of X, so we need to reorder
        # the axes. The next line does this, and should probably be included in
        # all subclasses.

        M = np.rollaxis(M, 0, len(M.shape))

        # The number of degrees of freedom for the masses. This should be a
        # one-dimensional array with the same number of entries as there are
        # masses.
        dof = np.array([1])

        # c is a constant for each particle used in the Coleman-Weinberg
        # potential using MS-bar renormalization. It equals 1.5 for all scalars
        # and the longitudinal polarizations of the gauge bosons, and 0.5 for
        # transverse gauge bosons.
        c = np.array([1.5])

        return M, dof, c

    def approxZeroTMin(self):
        # There are generically two minima at zero temperature in this model,
        # and we want to include both of them.
        return [np.array([self.v]), np.array([-self.v])]



colors = ['blue', 'green', 'orange', 'red', 'purple', 'black']

with open('data.json', 'r') as f:
  data = json.load(f)

with open('qsea-params.json', 'r') as f:
  params = json.load(f)

plots = data['runs']
phi = data['phi']
treeLevel = data['treeLevel']

Gamma = params['Gamma']
vev = params['v']

shift=0

for run in plots:
    ctr=0
    for potVal in run['pot']:
        if ctr == int(0.5*len(run['pot'])):
            shift=float(potVal)
        ctr+=1

    i = 0
    for potVal in run['pot']:
        run['pot'][i] = float(potVal-shift)
        i+=1

    shift =0

showError = False

if not showError:
    i=0
    for run in plots:
        plt.plot(phi,run['pot'], color=colors[i],label=r'$T={0}$'.format(run['T']))
        i+=1

    plt.plot(phi,treeLevel, label = 'Tree')


m = model1(vev, Gamma)
X = np.array([phi]).T

if not showError:
    plt.plot(X[...,0],m.V0(X)-m.V0([0]), label = 'Tree')
    i=0
    for run in plots:
        T = run['T']
        plt.plot(X[...,0],m.Vtot(X,T)-m.Vtot([0],T), linestyle ='dashed', color = colors[i], label=r'$T={0}$'.format(T))
        i+=1
else:
    i=0
    for run in plots:
        T = run['T']
        plt.plot(X[...,0],np.log10(np.abs((m.Vtot(X,T)-m.Vtot([0],T)-run['pot'])/run['pot'])), linestyle ='dashed', color = colors[i], label=r'$T={0}$'.format(T))
        i+=1

plt.grid(True, which='both')

plt.title(r'$m^2={0}$, $\lambda={1}$'.format(m.mPhi2, m.lPhi))
plt.xlim([-2,2])
plt.ylim([-0.04,0.02])
plt.legend(loc='upper right')
plt.show()
