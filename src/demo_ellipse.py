#!/usr/bin/python
"""
demo_ellipse.py

Fit an ellipse using different methods
"""

__author__ = 'rudolf.hoefler@gmail.com'
__copyright__ = 'WTFL'
__svn_id__ = '$Id$'

import numpy as np
from numpy.random import randn, rand
from matplotlib.pyplot import figure
from matplotlib.pyplot import show

from ellipse import Ellipse
from ellipse import EllipseBookstein
from ellipse import EllipseTrace
from ellipse import EllipseLevenberg
#from ellipse import EllipseGaussNewton
from ellipse import EllipsePatch
from ellipse import ellipse_polar
from ellipse import pol2cat

from matplotlib.font_manager import FontProperties
fprops =  FontProperties(size="small")

def fit_ellipse(x, y):
    """Compare different fit routines

    1) Ellipse
    2) EllipseBookstein
    3) EllipseTrace
    4) EllipseLevenberg
    """

    # fit raw data
    el  = Ellipse(x, y)
    elb = EllipseBookstein(x, y)
    elt = EllipseTrace(x, y)
    ell = EllipseLevenberg(x, y, p0=(elt.center[0], elt.center[1],
                                     elt.a, elt.b, elt.phi0))

    # plot fit in normal form
    el_hl1 = EllipsePatch( (0, 0), 2*el.a, 2*el.b, fill=False, color='g',
                           label="No Constraint")
    # plot fit of raw data
    el_raw = EllipsePatch( tuple(el.center), 2*el.a, 2*el.b,
                           angle=np.degrees(el.phi0), fill=False, color='r',
                           lw=1, label="No Constraint rotated")

    el_book = EllipsePatch( tuple(elb.center), 2*elb.a, 2*elb.b,
                            angle=np.degrees(elb.phi0), fill=False, color='y',
                            lw=1, label="Bookstein")
    el_trace = EllipsePatch( tuple(elt.center), 2*elt.a, 2*elt.b,
                         angle=np.degrees(elt.phi0), fill=False, color='m',
                         lw=1, label="Trace")
    el_lev = EllipsePatch(tuple(ell.center), 2*ell.a, 2*ell.b,
                          angle=np.degrees(ell.phi0), fill=False, color='c',
                          lw=1, label="Levenberg")

    fig, ax = figAxes()
    ax.plot(x, y, 'bo')
    ax.add_artist(el_hl1)
    ax.add_artist(el_book)
    ax.add_artist(el_raw)
    ax.add_artist(el_trace)
    ax.add_artist(el_lev)

    xline, yline = el_raw.paxes(ls='--', color='k')
    ax.add_artist(xline)
    ax.add_artist(yline)


    leg = ax.legend((el_raw, el_hl1, el_trace, el_book, el_lev),
                    ("No Constraint ", "No Constraint rotated",
                     "Trace", "Bookstein", "Levenberg"),
                    frameon=True, fancybox=True, shadow=False,
                    bbox_to_anchor=(1.0, 1.01), loc=2,
                    prop=fprops)

    leg.draggable(True, use_blit=True)
    leg.get_frame().set_alpha(0.9)
    fig.text(0.05, 0.40, el_lev.text(), fontsize=11)
    limits(ax)

    print_params("No Constraint", el)
    print_params("Trace", elt)
    print_params("Bookstein", elb)
    print_params("Levenberg", ell)


def print_params(title, el):
    print 30*'-'
    print title
    print 'Center: %8.3f, %8.3f' %tuple(el.center)
    print '  a, b: %8.3f, %8.3f' %(el.a, el.b)
    print ' Alpha: %8.3f' %el.phi0
    print '   Eps: %8.3f' %el.eps

def limits(ax):
    lmax = max(list(ax.get_xlim())+list(ax.get_ylim()))
    ax.set_ylim((-lmax, lmax))
    ax.set_xlim(ax.get_ylim())

def figAxes():
    fig = figure(figsize=(10, 6))
    ax = fig.add_subplot(111, aspect='equal')
    ax.axhline(y=0.0, color='k')
    ax.axvline(x=0.0, color='k')
    return fig, ax


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description=
                                     'Fit and plot ellipses from 2d data')
    parser.add_argument('-f', '--file', dest="file",
                        help='csv file', default=None)
    parser.add_argument('-n', '--noise', dest="noise", type=float,
                        help='Random noise', default=0.3)
    parser.add_argument('-p', '--points', dest="points", type=int,
                        help='Number of data points', default=100)
    args = parser.parse_args()

    if args.file:
        data= np.loadtxt(args.file, delimiter=',')
        x = data[1:, 0] - data[0, 0]
        y = data[1:, 1] - data[0, 1]
    else:
        err = randn(args.points)*args.noise
        phi = np.linspace(0, 2*np.pi, args.points)  # angles of the channels
        r0 = ellipse_polar(phi, 4, rand(), rand()*np.pi)
        x, y = pol2cat(phi, r0+err, deg=False)

    fit_ellipse(x, y)
    show()

