
import numpy as np
import numpy.linalg as la
from collections import namedtuple

FIELDS = ['xraw', 'yraw', 'a', 'b',
          'eps', 'center', 'nform',
          'popt', 'paxes', 'phi0']

FNAN = len(FIELDS)*[None]


class Ellipse(object):

    def __init__(self, x, y):
        
        self.fdata = namedtuple('fdata', FIELDS )
        self.fdata.xraw = x
        self.fdata.yraw = y
        
        w, v = self._fit_algebraic()
        self.fdata.popt = v.flatten()
        self.fdata.paxes = w

        # TODO - use property decorator
        self.normal_form()
        self.params()

    def _fit_algebraic(self):
        """Returns the design matrix for an algebraic fit of an ellipsis
        """

        SC = np.array([])
        for x, y in zip(self.fdata.xraw, self.fdata.yraw):
            SC = np.append(SC, np.array([x**2, y**2, 2.*x*y, x, y, 1.]))
        SC = np.matrix(SC.reshape(-1, 6))

        # design matrix
        DM = SC.T*SC
        w, v = la.eig(DM)
        w, v = self._lambda_min(w, v)

        return w, v

    def _lambda_min(self, w, v):
        """Returns the minimum eigenvalue and the corresponding eigenvector"""

        ind = np.where(w==w.min())[0]
        
        return np.array(w[ind]), np.array( v[:, ind])   

    def normal_form(self):
          
        par = self.fdata.popt
        A = np.matrix([[par[0], par[2]], [par[2], par[1]]])
        val, vec = la.eig(A)
        
        b = np.matrix(par[3:5]).T
        t = (-1./2.*b.T*A.I).T
        c = t.T*A*t + b.T*t + par[-1]
           
        
        self.fdata.center = np.array(t).flatten()
        self.fdata.nform = np.append(val, c)
        self.fdata.paxes = vec
        
        if c > 0.:
            self.fdata.nform *= -1.0

        return self.fdata.nform

    
    def paxes(self):

        pa = self.fdata.paxes
        c = self.fdata.center.T.flatten()

        scale = max((self.fdata.a, self.fdata.b))*1.1
        py = np.array(( [-1*pa[0,0], pa[0,0]], [-1*pa[1,0], pa[1,0]] )) 
        px = np.array(( [-1*pa[0,1], pa[0,1]], [-1*pa[1,1], pa[1,1]] )) 

        return px*scale+c[0], py*scale+c[1]


    def params(self):
        
        nform = self.fdata.nform
        pa = self.fdata.paxes

        a = np.sqrt(-nform[2]/nform[0])
        b = np.sqrt(-nform[2]/nform[1])

        if a > b:
            eps = np.sqrt(1.0-(b/a)**2)
            phi0 = np.arctan(pa[1, 0]/pa[0, 0])
        else:
            eps = np.sqrt(1.0-(a/b)**2)
            phi0 = np.arctan(pa[1, 1]/pa[0, 1])

        self.fdata.a = a
        self.fdata.b = b
        self.fdata.eps = eps
        self.fdata.phi0 = phi0


    def ellipse(self, n=300, rawfit=False):    

        phi = np.linspace(0, 2*np.pi, n)

        if self.fdata.a > self.fdata.b:
            p = self.fdata.b
        else:
            p = self.fdata.a
            
        rad = p/np.sqrt(1-(self.fdata.eps*np.cos(phi))**2)
        x, y = pol2cat(phi, rad, deg=False)

        if rawfit:
            x, y = rotate(x, y, self.fdata.phi0)
            x += self.fdata.center[0]
            y += self.fdata.center[1]
            
        return x, y


def rotate(x, y, phi):
        
    xr = x*np.cos(phi) - y*np.sin(phi)
    yr = x*np.sin(phi) + y*np.cos(phi)

    return xr, yr

def pol2cat(phi, rad, deg=True):

    if deg:
        phi = np.radians(phi)

    return rad * np.cos(phi), rad * np.sin(phi)
        

def ellipse_polar(phi, b, e, t):    
    
    rad = b/np.sqrt(1-(e*np.cos(phi-t))**2)
    
    return rad
