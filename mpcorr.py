#!/usr/bin/python

from numpy import ones, cos, matrix, pi, arange, degrees, zeros
from numpy import degrees, radians, array, round
from numpy.random import randn
from pylab import plot, show, xlim, figure,title, polar, gca
from pdb import set_trace
from numpy import linalg as LA

# suffixes
# dp->dipole 
# ch->channel of the mulitpole

def chVoltage(dp, phi):
    
    return cos(dp - phi)

sens = 1.0
npoles = 8.0
nmeas = 8.0
eamp = 0.2
err = eamp*randn(nmeas)
r0 = ones(npoles)*10

# ei = matrix(cos(phi))
# C = matrix(d0).T*ei

phi = arange(npoles)/npoles*pi  # angles of the channels
dpole = arange(nmeas)/nmeas*2*pi  # angles of the dipole fields

npoles = int(npoles)
nmeas = int(nmeas)

mx = array([chVoltage(dpole[i], phi[j]) for i in range(nmeas)
            for j in range(npoles)]).reshape(nmeas, npoles)

phi_ch = matrix(dpole).T*matrix(ones(npoles))
phi_dp = matrix(ones(nmeas)).T*matrix(phi)

A = matrix(round(cos(phi_dp-phi_ch), 5))
wadiag = abs(A).sum()/npoles
Aw = abs(A)/wadiag
w, v = LA.eig(A)

# x = matrix(err)*A.I
figure(1)
polar(phi, r0, '-b', lw=3)
polar(phi, r0 + err, 'ob')
gca().set_rmax(r0.max()*1.05)

#set_trace()

#polar(phi, d0 + x, 'or')
title('Deviation')

show()



