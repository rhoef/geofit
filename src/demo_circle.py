#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
demo_circle.py

Fit an ellipse using different methods.
"""

__author__ = 'rudolf.hoefler@gmail.com'
__copyright__ = 'WTFL'

from matplotlib import use
use("GtkAgg", warn=False)

import numpy as np
from numpy.random import randn
from matplotlib.pyplot import show

from ellipse import pol2cat
from scatter import figAxes
from scatter import limits
from scatter import fprops
from circle import CircleAlgebraic
from circle import CircleLevenberg
from circle import CirclePatch


def fit_circle(x, y):
    """Compare different fit routines

    1) Circle algebraic fit
    2) Circle geomeotric fit
    """

    # fit raw data
    alg = CircleAlgebraic(x, y)
    geo = CircleLevenberg(x, y, p0=(alg.center[0], alg.center[1],
                                    alg.radius))

    # plot fit in normal form
    circ_alg = CirclePatch(alg.center, alg.radius,fill=False, color='g',
                           label="Algebraic fit")
    # plot fit of raw data
    circ_geo = CirclePatch(geo.center, geo.radius, fill=False, color='r',
                         lw=1, label="Levenberg-Marquardt")

    # mpl figure
    fig, ax = figAxes()
    ax.plot(x, y, 'bo')
    ax.add_artist(circ_alg)
    ax.add_artist(circ_geo)
    ax.axvline(geo.center[0], ls='--', color='k')
    ax.axhline(geo.center[1], ls='--', color='k')

    leg = ax.legend((circ_alg, circ_geo),
                    ("Algebraic fit", "Levenberg-Marquardt"),
                    frameon=True, fancybox=True, shadow=False,
                    bbox_to_anchor=(1.0, 1.01), loc=2,
                    prop=fprops)

    leg.draggable(True, use_blit=True)
    leg.get_frame().set_alpha(0.9)
    fig.text(0.05, 0.40, geo.text(), fontsize=11)
    print_params("Algebraic fit", alg)
    print_params("Geomeotric fit", geo)

def print_params(title, alg):
    print 30*'-'
    print title
    print 'Center: %8.3f, %8.3f' %tuple(alg.center)
    print 'Radius: %8.3f' %(alg.radius)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=
                                     'Fit and plot ellipses from 2d data')
    parser.add_argument('-f', '--file', dest="file",
                        help='csv file', default=None)
    parser.add_argument('-r', '--radius', default=1.0, type=float,
                        help="Radius of the circle")
    parser.add_argument('-x' '--deltaX', dest="x", type=float,
                        default=0.0, help="Center position in x")
    parser.add_argument('-y' '--deltaY', dest="y", type=float,
                        default=0.0, help="Center position in y")
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
        x, y = pol2cat(phi, args.radius+err, deg=False)

    fit_circle(x+args.x, y+args.y)
    show()
