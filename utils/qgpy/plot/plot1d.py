from pylab import *
import numpy as np

from atmos import twoDslices

from qgpy.data import rdqg#,qgtools

def zonalMean(istart,istop,type='q',ilevel=0,superpose=False):
    """ Plot Zonal mean (x-average) of quantity type """
    ndump=istop-istart+1
    for idump in range(istart,istop+1):
        print "file number: ",idump
        data=rdqg.qgData(idump)
        var2d=data.getVar(type,ilevel)
        if idump==istart:
            y=data.y ; ny=y.size
            var1d=np.zeros(ny)
        var1d=var1d+np.mean(var2d,0)    
    var1d=var1d/ndump

    # Plot results
    if not(superpose):
        figure()
    plot(y,var1d)
    

def plotMeanVel(istart,istop):
    """ Plot upper and lower level zonal mean velocities """
    y,qU,psiU,uU,vU=qgtools.getMean(istart,istop,ilevel=0)
    y,qL,psiL,uL,vL=qgtools.getMean(istart,istop,ilevel=1)
    figure()
    plot(y/1.e3,uU,'k-',linewidth=3)
    plot(y/1.e3,uL,'k--',linewidth=3)
    plot(y/1.e3,40.*exp(-y**2/2.5e6**2),'k')
    xlim(-1.5e4,+1.5e4)
    ylim(-10,+50)
    xlabel('y (1000 km)')
    ylabel('U (m/s)')
