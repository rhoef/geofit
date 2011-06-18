#!/usr/bin/python

import numpy as np
import pylab as pl
from numpy.random import randn, rand
from matplotlib.pyplot import figure, show, draw

from pdb import set_trace
from numpy import linalg as LA
from matplotlib.patches import Circle

from pa import onb
from pa import ellipse
from pa import scatter_matrix
from pa import ellipse_polar
from pa import pol2cat
from pa import lambda_min
from pa import normal_form
from pa import mag
from pa import par_diag

def fit_ellipse(x, y):

    points = [(i, j) for i, j in zip(x, y)]
    S = scatter_matrix(points)
    DM = S.T*S
    w, v = LA.eig(DM)
    w, v = lambda_min(w, v)

    xf, yf = ellipse(v.flatten())
    pa, pl, D = onb(v.flatten())

    npar = par_diag(v.flatten(), pl, pa) 
    nform = normal_form(v.flatten(), pl, pa) 
    print v.flatten()
    print npar
    print nform

    fig = figure(1)
    ax = fig.add_subplot(111, aspect='equal') 
    ax.plot(x, y, 'bo')
    ax.plot(xf, yf, '-r', lw=2)
    ax.axhline(y=0.0, color='k')
    ax.axvline(x=0.0, color='k')
    
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
        
    r = np.sqrt((x*x + y*y))
    sc = r.mean()*1.3
        
    # gca().add_patch(Circle((0,0), r.mean(), fill=False, color='b', lw=1))
    ax.plot(np.array([-1*pa[0,0], pa[0,0]])*sc, np.array([-1*pa[1,0], pa[1,0]])*sc, 
            '-k', lw=2, ls='--')
    ax.plot(np.array([-1*pa[0,1], pa[0,1]])*sc, np.array([-1*pa[1,1], pa[1,1]])*sc, 
            '-k', lw=2, ls='--')

    if abs(xlim[0]- xlim[1]) > abs(ylim[0]-ylim[1]):
        ax.set_ylim(xlim)
    else:
        ax.set_xlim(ylim)

        
    ax.set_title('Deviation')
    draw()
       
    show()
      
if __name__ == '__main__':

    sens = 1.0
    nmeas = 1000
    err = randn(nmeas)*0.6
    phi = np.linspace(0, 2*np.pi, nmeas)  # angles of the channels

    # first line of data is the inital postion
    # data = np.loadtxt('dipoletest.csv', delimiter=',')
    # y = data[1:,0] - data[0,0]
    # x = data[1:,1] - data[0,1]

    r0 = ellipse_polar(phi, 170, rand(), rand()*np.pi)
    x, y = pol2cat(r0+err, phi, deg=False)
    
    fit_ellipse(x, y)
    
    show()
  
