import numpy as np
from pylab import *

from qgpy.data import rdqg

def getXi(data,zcut):
    "Return relative vorticity at level zcut"
    dx=data.x[1]-data.x[0] ; dy=data.y[1]-data.y[0]
    u=data.getVar('u',ilevel=zcut) ; v=data.getVar('v',ilevel=zcut)
    nx,ny=u.shape
    psibar=data.psibar[:,:,zcut]
    ubar=-(psibar[1:nx+1,2:ny+2]-psibar[1:nx+1:,0:ny])/dy/2.
    u=u-ubar
    xi=(np.roll(u,1,axis=1)-np.roll(u,-1,axis=1))/dy-\
       (np.roll(v,1,axis=0)-np.roll(v,-1,axis=0))/dx
    return xi[:,:]/2.e0

def getMeanStress(istart,istop,ilevel=0):
    """ Get Time average of uprime*vprime at ilevel"""
    ndump=istop-istart+1
    y,qbar,psibar,ubar,vbar=getMean(istart,istop,ilevel)
    stress=np.zeros(y.size)
    for idump in range(istart,istop+1):
        print idump
        data=rdqg.qgData(idump)
        u=data.getVar('u',ilevel)
        v=data.getVar('v',ilevel)
        uprime=getVarPrime(u,ubar)
        vprime=getVarPrime(v,vbar)
        stress=stress+np.mean(uprime*vprime,0)
    stress=stress/ndump
    #Plot result
    figure()
    plot(y/1.e3,stress)

def getHeatFlux(istart,istop):
    """ Get Time average of uprime*vprime at ilevel"""
    ndump=istop-istart+1
    y,qU,psiU,uU,vU=getMean(istart,istop,ilevel=0)
    y,qL,psiL,uL,vL=getMean(istart,istop,ilevel=1)
    stress=np.zeros(y.size)
    for idump in range(istart,istop+1):
        print idump
        data=rdqg.qgData(idump)
        u=data.getVar('u',ilevel=0)
        v=data.getVar('v',ilevel=0)
        psiprimeU=getVarPrime(data.getVar('psi',ilevel=0),psiU)
        psiprimeL=getVarPrime(data.getVar('psi',ilevel=1),psiL)
        uprime=getVarPrime(u,uU)
        vprime=getVarPrime(v,vU)
        stress=stress+np.mean(vprime*(psiprimeU-psiprimeL)/7.e5,0)
    stress=stress/ndump
    #Plot result
    figure()
    plot(y/1.e3,stress)

def getMean(istart,istop,ilevel=0):
    """ Get Time average of basic variables at ilevel"""
    ndump=istop-istart+1
    data=rdqg.qgData()
    for idump in range(istart,istop+1):
        print idump
        data.load(idump)
        if idump==istart:
            nx=data.x.size ; ny=data.y.size
            y=data.y
            qbar  =np.zeros((ny))
            psibar=np.zeros((ny))
            ubar  =np.zeros((ny))
            vbar  =np.zeros((ny))
        qbar  =  qbar+np.mean(data.getVar('q'  ,ilevel),0)
        psibar=psibar+np.mean(data.getVar('psi',ilevel),0)
        ubar  =  ubar+np.mean(data.getVar('u'  ,ilevel),0)
        vbar  =  vbar+np.mean(data.getVar('v'  ,ilevel),0)
    qbar=qbar/ndump ; psibar=psibar/ndump
    ubar=ubar/ndump ; vbar=vbar/ndump
    return y,qbar,psibar,ubar,vbar

def getVarPrime(u,ubar):
    """ Return 2D map of fluctuations for given variable stored in u """
    nx,ny=u.shape
    ubarXY=np.repeat(ubar[newaxis,:],nx).reshape(ny,nx)
    return transpose(transpose(u)-ubarXY)
