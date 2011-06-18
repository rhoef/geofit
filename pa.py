
import numpy as np
import numpy.linalg as la

def matA(par):
    
    A = np.matrix([ [ par[0], par[2] ],
                    [ par[2], par[1]]])
    
    return A
               
def par_diag(par, val, vec):
    
    npar = np.append(val, 0.0,)
    npar = np.append(npar, par[3:5]*vec.T )
    npar = np.append(npar, par[-1]) 

    return npar

def normal_form(opar, val, vec):

    par = par_diag(opar, val, vec)
    c0 = par[-1] - par[3]**2/(2*par[0])**2 - par[4]**2/(2*par[1])**2
    print np.array([np.sqrt(abs(1/par[0])), np.sqrt(abs(1/par[1])), c0])
    return np.array([np.sqrt(abs(c0/par[0])), np.sqrt(abs(c0/par[1])), c0])



def onb(par):
    
    A = matA(par)
    val, vec = la.eig(A)
    D = vec*A*vec.T

    return vec, val, D


def xlimits(a, b, c, d, e, f):
    
    p = (c*e - 2*d*b)/(2*(c**2 - b * a))
    q = (e**2 - 4 * b * f)/(4*(c**2 - b * a))
    
    xmin = -p/2 - np.sqrt((p/2)**2 - q)
    xmax = -p/2 + np.sqrt((p/2)**2 - q)

    return xmin*0.999, xmax*0.999


def ellipse(par, n=300):

    (a, b, c, d, e, f) = par.flatten()
    xlim = xlimits(a, b, c, d, e, f)
    x = np.linspace(xlim[0], xlim[1], n)
    
    p = (2*c*x+e)/b
    q = (a*x**2 + d*x +f)/b    
    y = -p/2 + np.sqrt((p/2)**2 - q) # upper half
    y = np.append(y, (-p/2 - np.sqrt((p/2)**2 - q))[::-1]) # lower half

    x = np.append(x, x[::-1])

    return x, y

def ellipse_polar(phi, b, e, t):    

    rad = b/np.sqrt(1-(e*np.cos(phi-t))**2)
    
    return rad


def lambda_min(w, v):
    
    ind = np.where(w==w.min())[0]
    
    return np.array(w[ind]), np.array( v[:, ind])   


def pol2cat(rad, phi, deg=True):

    if deg:
        phi = np.radians(phi)

    return rad * np.cos(phi), rad * np.sin(phi)

def scatter_matrix(points):
    """Returns the design matrix for an algebraic fit of an ellipsis
    """
    
    DM = np.array([])
    for (x, y) in points:
        DM = np.append(DM, np.array([x**2, y**2, 2*x*y, x, y, 1]))
#        DM = np.append(DM, np.array([x**2, y**2, x*y, x, y, 1]))
    DM = np.matrix(DM.reshape(-1, 6))
        
    return DM

def mag(x, axis=None):

    return np.sqrt( (x**2).sum(axis=axis))
