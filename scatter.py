#!/usr/bin/python

import numpy as np
import pylab as pl
from numpy.random import randn, rand
from matplotlib.pyplot import figure, show, draw, gca

from pdb import set_trace
from numpy import linalg as LA
from matplotlib.patches import Circle

from ellipse import Ellipse
from ellipse import EllipsePatch
from ellipse import ellipse_polar
from ellipse import pol2cat


def fit_ellipse(x, y):

    el  = Ellipse(x, y)
    el_hl1 = EllipsePatch( (0, 0), 2*el.a, 2*el.b, fill=False, color='g')
    el_raw = EllipsePatch( tuple(el.center), 2*el.a, 2*el.b, 
                           angle=np.degrees(el.phi0), fill=False, color='r', )


    print 'Angle: ', np.degrees(el.phi0)
    print 'a:', el.a
    print 'b:', el.b
    print 'Excentricity:', el.eps


    fig = figure(1)
    ax = fig.add_subplot(111, aspect='equal') 
    ax.plot(x, y, 'bo')
    ax.axhline(y=0.0, color='k')
    ax.axvline(x=0.0, color='k')
    
    ax = gca()
    ax.add_artist(el_hl1)
    ax.add_artist(el_raw)
    xline, yline = el_raw.paxes(ls='--', color='k')
    ax.add_artist(xline)
    ax.add_artist(yline)

    
    lmax = max(list(ax.get_xlim())+list(ax.get_ylim()))
    ax.set_ylim((-lmax, lmax))
    ax.set_xlim(ax.get_ylim())
    ax.set_title('Deviation')
       
    show()
      

if __name__ == '__main__':

  # sens = 1.0
  # nmeas = 100
  # err = randn(nmeas)*0.1
  # phi = np.linspace(0, 2*np.pi, nmeas)  # angles of the channels
  # r0 = ellipse_polar(phi, 4, rand(), rand()*np.pi)
  # x, y = pol2cat(phi, r0+err, deg=False)

    data= np.loadtxt('dipoletest.csv', delimiter=',')
    x = data[1:,0] - data[0,0]
    y = data[1:,1] - data[0,1]
    
    
    fit_ellipse(x+1, y+2)
    
    show()
  
