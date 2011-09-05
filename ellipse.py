
import numpy as np
import numpy.linalg as la
from collections import namedtuple
from matplotlib import patches
from matplotlib import lines
from matplotlib.pyplot import arrow
from scipy.optimize import leastsq

FIELDS = ['xraw', 'yraw', 'a', 'b',
          'eps', 'center', 'nform',
          'popt', 'paxes', 'phi0']

FNAN = len(FIELDS)*[None]


class EllipsePatch(patches.Ellipse):

    def __init__(self, xy, width, height, angle=0.0, **kwargs):
        
        if width < height:
            width, height = height, width

        super(EllipsePatch, self).__init__(xy, width, height, angle, **kwargs)

        self.angle = angle
        self.width = width
        self.height = height
        self.center = xy

    def paxes(self, *args, **kwargs): 

        cx, cy = self.center                
        sc = (self.width + self.height)*0.5
        ssin = np.sin(np.radians(self.angle))*sc
        scos = np.cos(np.radians(self.angle))*sc

        xline = lines.Line2D( (-scos+cx, scos+cx),
                              (-ssin+cy, ssin+cy), *args, **kwargs)
        
        yline = lines.Line2D( (-ssin+cx, ssin+cx),
                              (scos+cy, -scos+cy), *args, **kwargs)
        
        return xline, yline
        
    def text(self):
            
        txt = 'Parameters:\n'
        txt += 'cx: %.3f\n' %self.center[0]
        txt += 'cy: %.3f\n' %self.center[1]
        txt += 'a: %.3f\n' %(self.width/2.0)
        txt += 'b: %.3f\n' %(self.height/2.0)
        txt += 'eps: %.3f' %self.eps

        return txt

    @property
    def eps(self):

        if self.width < self.height:
            eps = np.sqrt(1- self.width**2/self.height**2)
        else:
            eps = np.sqrt(1- self.height**2/self.width**2)

        return eps


class Ellipse(object):



    def __init__(self, x, y):
        
        self.fdata = namedtuple('fdata', FIELDS )
        self.fdata.xraw = x
        self.fdata.yraw = y
        
        w, v = self._fit_algebraic()
        self.fdata.popt = v.flatten()
        self.fdata.paxes = w

        # TODO - use property decorator
        self._normal_form()
        self._params()


    def __getattr__(self, attr):
        
        if hasattr(self.fdata, attr) and not attr.startswith('_'):
            return getattr(self.fdata, attr)
        else:
            return getattr(self, attr)

    def _fit_algebraic(self):
        """Returns the design matrix for an algebraic fit of an ellipsis
        """

        SC = np.array([])
        for x, y in zip(self.fdata.xraw, self.fdata.yraw):
            SC = np.append(SC, np.array([x**2, y**2, 2.*x*y, x, y, 1.]))
        SC = np.matrix(SC.reshape(-1, 6))

        # design matrix
        DM = SC.T*SC
        import pdb; pdb.set_trace()
        w, v = la.eig(DM)
        w, v = self._lambda_min(w, v)

        return w, v

    def _lambda_min(self, w, v):
        """Returns the minimum eigenvalue and the corresponding eigenvector"""

        ind = np.where(w==w.min())[0]
        
        return np.array(w[ind]), np.array( v[:, ind])   

    def _normal_form(self):
          
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


    def _params(self):
        
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


class EllipseGeometric(object):
    
    def __init__(self, x, y, p0=None):

        self.x = x
        self.y = y
        self._p0 = p0
        self.fit()
        
    def mfunc(self, p, x, y):
            
        phi = self.phi() - p[4]
        # cx, cy, a, b, alpha, phi
        xbar = p[2]*np.cos(phi)
        ybar = p[3]*np.sin(phi)
#        rotm = self.rotmat(p[4])
        rotx = np.cos(p[4])*xbar - np.sin(p[4])*ybar
        roty = np.sin(p[4])*xbar + np.cos(p[4])*ybar

        print x.shape
        print rotx.shape
        print y.shape
        print roty.shape

#        import pdb; pdb.set_trace()

       # di = np.sqrt((p[0]-x)**2 + (p[1]-y)**2) - np.sqrt(rotx**2 + roty**2)
        di = (x - p[0])**2 + (y -p[1])**2  - (rotx**2 + roty**2)
        return di

    def func(self, p, x, y):

        return (x-p[0])**2/p[2]**2 + (y-p[1])**2/p[3]**2 - 1

    def chi(self, x, y):

        return lambda p: self.func(p, x, y).flatten()

    def xbar(self, a, b, phi):

        return a*np.cos(phi)
        
    def rotmat(self, alpha):

        return np.matrix([np.cos(alpha), -np.sin(alpha)],
                         [np.sin(alpha), np.cos(alpha)])

    def fit(self):
    
        efunc = self.chi(self.x, self.y)
        p0 = self.p0
        popt,err = leastsq(efunc, p0, full_output=False)
#        import pdb; pdb.set_trace()

        self.center = popt[:2]
        self.a = popt[2]
        self.b = popt[3]
        self.phi0 = np.degrees(popt[4])

    @property
    def p0(self):
        
        if not self._p0:
        #        p0 = np.append(np.array([0.0, 0.0, 111, 112, 0.0]), phi0)
            return np.append( [0.0, 0.0, 111, 112, 0.0], self.phi)
        else:
            return self._p0
#           return np.append(self._p0, self.phi())

    def phi(self, x0=0.0, y0=0.0):
        
        return np.arctan2(self.y-y0, self.x-x0)

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
