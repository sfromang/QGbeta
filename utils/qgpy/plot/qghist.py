import matplotlib.pyplot as plt
from pylab import *
import numpy as np
import os

from qgpy.data import rdqg

class history:

    def __init__(self,dirname='./',filename='history.txt',verbose=True):
        """ Define class to read history file """
        self.dirname=dirname
        self.filename=filename
        self.verbose=verbose
        return

    def load(self):
        """ Read history data """
        if (self.verbose):
            print '   Reading history data...'
        filename=self.dirname+self.filename
        self.time,self.dt,self.uprobe,self.umask,self.umean=\
            np.loadtxt(filename,unpack='true',usecols=(0,1,2,3,4))

    def plot(self,pltName=None,clPlt=False,superpose=False,\
                savePlt=False,saveFileName='history.png',Trot=1.):
        """ Plot history data """
        if not(superpose):
            figure()
        plt.plot(self.time/Trot,self.umean)
        plt.plot(self.time/Trot,smooth(self.umean,100))
        plt.xlabel('Time')
        plt.ylabel('Mean zonal velocity')
        if pltName is not None:
            plt.title(pltName)
        saveAndclean(savePlt=savePlt,filename=self.dirname+saveFileName,clPlt=clPlt)
        return

    def plotSmooth(self,pltName=None,clPlt=False,superpose=False,\
                       savePlt=False,saveFileName='history.png',Trot=1.):
        """ Plot history data """
        if not(superpose):
            figure()
        plt.plot(self.time/Trot,smooth(self.umean,100))
        plt.xlabel('Time')
        plt.ylabel('Mean zonal velocity')
        if pltName is not None:
            plt.title(pltName)
        saveAndclean(savePlt=savePlt,filename=self.dirname+saveFileName,clPlt=clPlt)
        return

def saveAndclean(savePlt=False,filename='file.png',clPlt=True):
    """ Small little function to save and clean a figure """
    if (savePlt):
        plt.savefig(filename)
    if (clPlt):
        plt.clf()
        plt.close()
    return    

def getDir():
    """ Return a list with the subdirectories contained 
    in the current directory """
    return [f for f in os.listdir('.') if os.path.isdir(f)]

def windHist(istart,istop,superpose=False,Tday=86400,wLog=False):
    """ Plot mean zonal wind time history """
    ndump=istop-istart+1 ; 
    time=np.zeros(ndump) ; vHist=np.zeros(ndump)
    for idump in range(istart,istop+1):
        print idump
        #time[idump-istart],vHist[idump-istart]=getMaxV(idump)
        #time[idump-istart],vHist[idump-istart]=getMeanV(idump)
        time[idump-istart],vHist[idump-istart]=getVpos(idump,9.59,0.)
    if not(superpose):
        figure()
    if wLog:
        semilogy(time[1:]/Tday,vHist[1:])
    else:
        plot(time[1:]/Tday,vHist[1:])


def getSigma(istart,istop,verbose=False):
    """ Return growth rate """
    ndump=istop-istart+1 ; Tday=86400. ; km=1.e3
    time=np.zeros(ndump) ; vHist=np.zeros(ndump)
    for idump in range(istart,istop+1):
        if (verbose):
            print idump
        time[idump-istart],vHist[idump-istart]=getMaxV(idump)
    sigma=(np.log(vHist)-np.roll(np.log(vHist),-1))/(time-np.roll(time,-1))
    return np.mean(sigma[1:])

def getMeanZonalWind(idump):
    """ Compute mean zonal wind at y=0 from output data """
    data=rdqg.qgData(idump)
    ny=data.y.size-2
    return data.time,np.mean(np.mean(data.getVar('u'),0)[ny/2:ny/2+2])

def getMaxV(idump):
    """ Compute maximum of y-vel from output data """
    data=rdqg.qgData(idump)
    return data.time,np.max(data.getVar('v'))

def getMeanV(idump):
    """ Compute maximum of x-vel from output data """
    data=rdqg.qgData(idump)
    return data.time,np.mean(data.getVar('u'))

def getVpos(idump,x,y):
    """ Compute maximum of x-vel from output data at given location """    
    data=rdqg.qgData(idump)
    i=np.int((x-data.x.min())/(data.x[1]-data.x[0]))
    j=np.int((y-data.y.min())/(data.y[1]-data.y[0]))
    return data.time,data.getVar('u',ilevel=0)[i,j]
    #return data.time,(data.getVar('u',ilevel=0)[i,j]+data.getVar('u',ilevel=1)[i,j])/2.

def annotatePlt():
    """ Annotate history plot """
    xlabel('t/Trot')
    ylabel('Velocity')

    Trot=4.*pi
    xblock=np.array([17000.,32000.,38000.,48000.]) ; xblock=xblock/Trot
    xplt=np.array([35500.,43000.])/Trot
    ylines=[-0.1,0.3]

    # Mark limits of blocked flows
    plot([xblock[0],xblock[0]],ylines,':k',linewidth=1)
    plot([xblock[1],xblock[1]],ylines,':k',linewidth=1)
    fill_between([xblock[0],xblock[1]],[0.3,0.3],[-0.1,-0.1],color='lightgray')

    plot([xblock[2],xblock[2]],ylines,':k',linewidth=1)
    plot([xblock[3],xblock[3]],ylines,':k',linewidth=1)
    fill_between([xblock[2],xblock[3]],[0.3,0.3],[-0.1,-0.1],color='lightgray')

    plot([xplt[0],xplt[0]],ylines,'--k',linewidth=2)
    plot([xplt[1],xplt[1]],ylines,'--k',linewidth=2)

    xlim(0.,6500.)
    ylim(-0.1,0.3)

def getHist(Trot=1.,tmin=0.,dirname='./',verbose=True,\
                doPlt=True,pltName=None,clPlt=False,superpose=False,\
                savePlt=False,saveFileName='history.png'):
    """ Read and mean velocity from history file. 
    Options available to plot and save data """

    if (verbose):
        print '   Reading history data...'
    filename=dirname+'history.txt'
    time,dt,uprobe,umask,umean=\
            np.loadtxt(filename,unpack='true',usecols=(0,1,2,3,4))
    if (verbose):
        print '   Mean velocity:',np.mean(umean[time>tmin])
    
    if (doPlt):
        if not(superpose):
            figure()
        plt.plot(time/Trot,umean)
        plt.plot(time/Trot,smooth(umean,100))
        plt.xlabel('Time')
        plt.ylabel('Mean zonal velocity')
        if pltName is not None:
            plt.title(pltName)
    saveAndclean(savePlt,dirname+saveFileName,clPlt)

    return umean[time>tmin]

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def regimeStat():
    A=[1.2,1.4,1.45,1.5,1.6,1.7,1.9]
    relLow=[1.,0.957600,0.801951,0.625173,0.19029,0.0114,0.]
    figure()
    plot(A,relLow,'s')
    plot([1.,1.37],[1.,1.],'k--')
    plot([1.37,1.67],[1.,0.],'k--')
    plot([1.67,2.],[0.,0.],'k--')
    xlim(1.,2.)
    ylim(-0.1,1.2)
    xlabel('Forcing Amplitude (proxy for Ro number)')
    ylabel('Time spent in blocked state')
    return

