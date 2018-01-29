from pylab import *

import numpy as np

from qgpy.data import rdqg
from atmos import twoDslices
from atmos.utils import str_suffix

def makePlot(y):
    """ Add information about Topo """
    dxB0=0.68e0
    x0=4.1125e0 ; x1=12.3375e0
    
    ymin=y.min() ; ymax=y.max()
    plot([x0,x0],[ymin,ymax],'k')
    plot([x0+dxB0,x0+dxB0],[ymin,ymax],'k--')
    plot([x0-dxB0,x0-dxB0],[ymin,ymax],'k--')
    plot([x1,x1],[ymin,ymax],'k')
    plot([x1+dxB0,x1+dxB0],[ymin,ymax],'k--')
    plot([x1-dxB0,x1-dxB0],[ymin,ymax],'k--')

def uMean(istart,istop,type='u',superpose=False):
    """ Plot time evolution of u at midchannel """
    umean=np.zeros(istop-istart+1)
    time=np.zeros(istop-istart+1)
    uabsmean=0.
    for idump in range(istart,istop+1):
        print "File number: ",idump
        data=rdqg.qgData(idump=idump)
        nx,ny=data.getVar(type).shape
#        umean[idump-istart]=np.mean(data.getVar(type)[:,ny/2],0)
        umean[idump-istart]=np.mean(np.mean(data.getVar(type)[:,:],0),0)
        #umean[idump-istart]=data.u[0,ny/2,0]
        uabsmean=uabsmean+np.mean(np.abs(data.u))
        #uabsmean=uabsmean+np.mean(data.getVar(type))
        time[idump-istart]=data.time
    if not(superpose):
        figure()
    plot(time,umean)
    uabsmean=uabsmean/(istop-istart+1)
    return uabsmean                         

def psiMean2D(istart,istop):
    """ Plot time evolution of u at midchannel """
    for idump in range(istart,istop+1):
        print "File number: ",idump
        data=rdqg.qgData(idump=idump)
        if (idump==istart):
            nx,ny,nlayer=data.psi.shape ; nx=nx-2 ; ny=ny-2
            psiMean=np.zeros((nx,ny))
        psiMean=psiMean+data.getVar(type='psi',ilevel=0)
    psiMean=psiMean/(istop-istart+1)

    #Define ppplot object
    plotObj=twoDslices.initObj()
    plotObj.x=data.x ; plotObj.y=data.y
    plotObj.f=psiMean
    #Axis labels
    plotObj.xlabel='X'
    plotObj.ylabel='Y'
    #Plot array
    twoDslices.plot2d(plotObj,figAspect=False)

def polarCont(istart,istop,ilevel=0,threshold=None,type='Z'):
    """ Make polar contour plot of psi """

    # Get data
    icount=0
    for idump in range(istart,istop+1):
        go=0
        print "File number: ",idump
        data=rdqg.qgData(idump=idump)
        if (idump==istart):
            nx,ny,nlayer=data.psi.shape
            psi=np.zeros((ny,nx))
        if (threshold is not(None)):
            if type=='Z':
                if np.mean(data.u[:,:,ilevel])>threshold:
                    icount=icount+1
                    go=1
            if type=='B':
                if np.mean(data.u[:,:,ilevel])<threshold:
                    icount=icount+1
                    go=1
        else:
            go=1
            icount=icount+1
        psi=psi+go*(data.psi[:,:,ilevel]+data.psibar[:,:,ilevel]).transpose()
    print icount
    psi=psi/icount
    psi=np.flipud(psi)

    ydim,xdim=psi.shape

    # extend of the box
    rMin=1. ; rMax=4.
    extent=[rMin,rMax,-np.pi/2.,3./2.*np.pi]

    # Compute plot coordinates
    xp,yp,ax=getCoord(xdim,ydim,extent)
    
    # Make plot
    #map1 = ax.pcolormesh(xp,yp,psi)
    map1 = ax.contour(xp,yp,psi,20,colors='k',linestyles='solid')

    # Add inner & outer circles
    addCircles(ax,rMin,rMax)

    # Add contour of topography
    addTopo(ax,extent,xp,yp,ydim,xdim)

    return

def getCoord(xdim,ydim,extent):
    """ Compute plot coordinates for polar plots """

    # Define array of figures
    fig, ax = subplots(1, 1, figsize=(8., 8.))
    subplots_adjust(bottom=0.15, left=0.15, wspace=0.5)

    rtile  = np.ones(xdim)
    ttile  = np.ones(ydim)
    dtheta = (extent[3] - extent[2])/xdim
    r      = np.linspace(extent[0], extent[1], ydim)
    theta  = np.linspace(extent[2], extent[3]  \
                             , xdim)
    theta  = np.linspace(extent[2] - dtheta/2., extent[3] + dtheta/2. \
                             , xdim)
    r      = r[:, np.newaxis]*rtile
    theta  = theta*ttile[:, np.newaxis]
    xp     = r*np.cos(theta)
    yp     = r*np.sin(theta)

    return xp,yp,ax

def addCircles(ax,rMin,rMax):
    """ Add inner & outer circles """
    eps=1.e-6
    x0=arange(-rMin,rMin+eps,2.*rMin/1000.)
    ax.plot(x0,np.sqrt((rMin+eps)**2-x0**2),'k',linewidth=3)
    ax.plot(x0,-np.sqrt((rMin+eps)**2-x0**2),'k',linewidth=3)
    x0=arange(-rMax,rMax+eps,2.*rMax/1000.)
    ax.plot(x0,np.sqrt((rMax+eps)**2-x0**2),'k',linewidth=3)
    ax.plot(x0,-np.sqrt((rMax+eps)**2-x0**2),'k',linewidth=3)
    return

def addTopo(ax,extent,xp,yp,nr,ntheta):
    """ Add topography shaded in grey """
    r      = np.linspace(extent[0],extent[1],nr)
    theta  = np.linspace(extent[2],extent[3],ntheta)
    topo=np.zeros((nr,ntheta))
    angleTopo=72. ; phiTopo=angleTopo/2./180.*pi
    for itheta in range(ntheta):
        if (abs(theta[itheta])<phiTopo):
            topo[:,itheta]=1.
        if (abs(theta[itheta]-pi)<phiTopo):
            topo[:,itheta]=1.
    ax.contourf(xp,yp,topo,1,colors=['white','lightgrey'])
    return

def addExtrema(ax,psi):
    """ Add letters at location of extrema """
    return

def savePNG(istart,istop):
    """ Save images to png file """
    for idump in range(istart,istop+1):
        polarCont(idump,idump)
        savefig('pngPsi/psi_'+str_suffix(idump,4),)
        clf()
        close()

   
#
#
#
def poissoncheck(idump=0):
    data=rdqg.qgData(idump)
    x=data.x ; y=data.y
    dx=data.x[1]-data.x[0]
    dy=data.y[1]-data.y[0]
    psi=data.psibar[:,:,0]
    nx,ny=psi.shape ; nx=nx-2 ; ny=ny-2
    qstar=data.qbar[1:nx+1,1:ny+1,0]
    q=np.zeros((nx,ny))
    for j in range(1,ny+1):
        for i in range(1,nx+1):
            q[i-1,j-1]=(psi[i+1,j]-2.*psi[i,j]+psi[i-1,j])/dx**2+\
                       (psi[i,j+1]-2.*psi[i,j]+psi[i,j-1])/dy**2
#            q[i-1,j-1]=(psi[i,j+1]-2.*psi[i,j]+psi[i,j-1])/dy**2
    
    figure()
    plot(y,qstar[1,:],'b')
    plot(y,q[1,:],'r')
    figure()
    plot(x,qstar[:,5],'b')
    plot(x,q[:,5],'r')

#    plot(y,qstar[1,:]-q[1,:])
#    figure()
#    plot(y,psi[0,1:ny+1])
