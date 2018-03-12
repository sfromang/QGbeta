from pylab import *
import matplotlib.pyplot as plt
import os

from atmos import twoDslices
from utils.utils import str_suffix

from qgpy.data import rdqg
from qgpy.plot import qgtools

def cut2dXY(idump,type='q',ilevel=0,minVal=None,maxVal=None,km=None):
    """2D plot of cut in XY plane"""

    # Define constants
    if (km is None):
        km=1.e3

    # Read data
    data=rdqg.qgData()
    data.load(idump)

    #Define ppplot object
    plotObj=twoDslices.initObj()

    #Load arrays
    plotObj.x=data.x/km ; plotObj.y=data.y/km
    plotObj.f=data.getVar(type=type,ilevel=ilevel)

    #Axis labels
    plotObj.xlabel='X'
    plotObj.ylabel='Y'

    #Plot boundary
    if (minVal is not(None)):
        plotObj.vmin=minVal
    if (maxVal is not(None)):
        plotObj.vmax=maxVal

    #Plot array
    twoDslices.plot2d(plotObj,figAspect=False)

def xiContour(idump,type='q',zcut=0,pltMotion=False):
    """2D plot of cut in XY plane"""

    # Define constants
    km=1.e3
    km=1.e0

    # Get relative vorticity
    data=rdqg.qgData(idump)
#    xi=qgtools.getArray(data,type='xi',zcut=zcut)# - xi0
    xi=data.getVar(type=type,ilevel=0)
    print xi.shape,data.x.shape
    print 'Time (in h)=',data.time/3600.
    print xi.min(),xi.max()

    # window with correct size
    x=data.x/km ; y=data.y/km
    xsize=np.max(x)-np.min(x)
    ysize=np.max(y)-np.min(y)
    w,h=figaspect(float(ysize/xsize))
    figure(figsize=(w,h))

    #Plot contour
#    levels=(-5.e-5,-2.5e-5,-1.5e-5,1.3e-5,2.5e-5,5.e-5,8.e-5)
# Next 4 lines for Gwendal stuff
#    levels=(-1.e-4,-7.5e-5,-5.e-5,-2.5e-5,2.5e-5,5.e-5,7.5e-5,1.e-4)
#    xiCont=plt.contour(x,y,transpose(xi),levels,colors='k')
#    plot([x.min(),x.max()],[0.,0.],color='k',linewidth=0.5)
#    plot([x.min(),x.max()],[-1000.,-1000.],':',color='k',linewidth=0.5)
    matplotlib.rcParams['contour.negative_linestyle']='dashed'
#    levels=(-0.4,-0.35,-0.3,-0.20,-0.15,-0.10,0.0,+0.05,+0.1,+0.15,+0.2)
#    levels=np.arange(xi.min(),xi.max(),(xi.max()-xi.min())/20.)
#    print levels
#    xiCont=plt.contour(x,y[abs(y)<1.57],transpose(xi[:,abs(y)<1.57]),levels,colors='k')
    xiCont=plt.contour(x,y,transpose(xi),15,colors='k')
#    plot([x.min(),x.max()],[0.,0.],color='k',linewidth=0.5)
#    plot([x.min(),x.max()],[-1000.,-1000.],':',color='k',linewidth=0.5)
    xlabel('X (in km)')
    ylabel('Y (in km)')
    if pltMotion:
        xvort,yvort=getMaxVorticity(0,idump,zcut)
        plot(xvort,yvort,linewidth=2,color='k')
    print xi.min(),xi.max()

def meanContour(start,stop,type='psi'):
    """ Plot contour of quantity type averaged in timpe """

    # Compute averaged quantity
    for idump in range(start,stop+1):
        data=rdqg.qgData(idump)
        if (idump==start):
            var=data.getVar(type=type,ilevel=0)
        else:
            var=var+data.getVar(type=type,ilevel=0)
    var=var/(stop-start+1)

    # window with correct size
    x=data.x ; y=data.y
    xsize=np.max(x)-np.min(x)
    ysize=np.max(y)-np.min(y)
    w,h=figaspect(float(ysize/xsize))
    figure(figsize=(w,h))

    matplotlib.rcParams['contour.negative_linestyle']='dashed'
    xiCont=plt.contour(x,y,transpose(var),15,colors='k')
    xlabel('X')
    ylabel('Y')

def getMaxVorticity(start,stop,zcut):

    # Get maximum vorticity coordinates
    xmax=np.zeros(stop-start+1)
    ymax=np.zeros(stop-start+1)
    for idump in range(start,stop+1):
        data=rdqg.qgData(idump)
        if (idump==start):
            nx=data.x.size ; ny=data.y.size
        xi=qgtools.getArray(data,type='xi',zcut=zcut) #- xi0
        xmax[idump-start]=data.x[np.argmax(xi)/ny]/1.e3
        ymax[idump-start]=data.y[np.mod(np.argmax(xi),ny)]/1.e3
    return xmax,ymax

def animq(istart,istop,saveImage=False,minVal=None,maxVal=None):
    for idump in range(istart,istop+1):
        print idump
        cut2dXY(idump,type='u',ilevel=0,minVal=minVal,maxVal=maxVal)
        # data=rdqg.qgData(idump,filedir='./')
        # vortex.x=data.x/km ; vortex.y=data.y/km
        # vortex.f=data.xi-xi0
        # if (idump==istart):
        #     vortex.vmin=vortex.f.min()
        #     vortex.vmax=vortex.f.max()
        # twoDslices.plot2d(jet,figAspect=True)
        # twoDslices.plot2d(vortex,figAspect=True,superpose=True)
        draw()
        if saveImage:
            savefig('png/anim'+str_suffix(idump,length=4)+'.png')
            clf()
            close()

