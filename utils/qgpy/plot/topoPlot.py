import numpy as np
import matplotlib.pyplot as plt

from qgpy.plot import qghist,histoFit

def analyseSeries(xlim=(0.,0.15),ylim=(0.,100.),tmin=0):
    """ Make a global analysis of the series """
    mydirectories=qghist.getDir() ; ndir=len(mydirectories)
    idir=0
    for currdir in mydirectories:
        print 'Directory:',currdir
        print '   Number of remaining directories: ',ndir-idir ; idir=idir+1
        pltName='A='+currdir[1:5]+' & noiseAmp='+currdir[13:]
        history=qghist.history(dirname=currdir+'/',verbose=False)
        history.load() ; umean=history.umean[history.time>tmin]
        history.plot(clPlt=True,savePlt=True,\
                         pltName=pltName,saveFileName=currdir+'_history.png')
        plotHisto(umean,pltHistoFit=True,dirname=currdir+'/',\
                      xlim=xlim,ylim=ylim,clPlt=True,savePlt=True,\
                      pltName=pltName,saveFileName=currdir+'_histogram.png')
    return

def plotHisto(umean,pltHistoFit=False,dirname='./',xlim=None,ylim=None,clPlt=False,
              pltName=None,savePlt=False,saveFileName='histogram.png'):
    plt.figure()
    plt.hist(umean,bins=30,normed=True)
    plt.xlabel('Mean zonal velocity')
    plt.ylabel('Relative frequency')
    if (pltHistoFit):
        histoFit.histoFit(umean,doPlt=True)
    if (xlim is not None):
        plt.xlim(xlim)
    if (ylim is not None):
        plt.ylim(ylim)
    if (pltName is not None):
        plt.title(pltName)
    qghist.saveAndclean(savePlt,dirname+saveFileName,clPlt)
    return

def getRegime(A=None,eps0=None,ucrit=1.,tmin=0.,superpose=False):
    """ Compute regime prevalence based on history data """
    mydirectories=qghist.getDir() ; ndir=len(mydirectories)
    idir=0
    if (eps0 is not None):
        Avals=np.zeros(1)
    for currdir in mydirectories:
        ACondition=(np.float(currdir[1:5])==A) and (A is not None)
        noiseCondition=(np.float(currdir[13:])==eps0) and (eps0 is not None)
        if (ACondition or noiseCondition):
            print 'Directory:',currdir
            print '   Reading history data...'
            history=qghist.history(dirname=currdir+'/',verbose=False)
            history.load() ; umean=history.umean[history.time>tmin]
            A0,A1,x0,x1,dx=histoFit(umean,doPlt=False)
            if (x0+dx>x1-dx):
                if (umean.mean()<ucrit):
                    blockedFrac=1.
                else:
                    blockedFrac=0.
            else:
                blockedFrac=A0/(A0+A1)
                
            print '      Relative time in blocked state= ',blockedFrac
            print '      Relative time in zonal state= ',1.-blockedFrac
            if (eps0 is not None):
                if idir==0:
                    Avals=np.float(currdir[1:5])
                    blocked=blockedFrac
                    idir=idir+1
                else:
                    Avals=np.append(Avals,np.float(currdir[1:5]))
                    blocked=np.append(blocked,blockedFrac)
    if not(superpose):
        plt.figure()
    plt.plot(Avals[np.argsort(Avals)],blocked[np.argsort(Avals)],'o--')
    plt.xlabel('Forcing amplitude A')
    plt.ylabel('Relative time spent in blocked state')
    plt.ylim(-0.1,1.2)
    return Avals[np.argsort(Avals)],blocked[np.argsort(Avals)]

def meanRegime(epsMin=0.6,epsMax=1.1,epsInt=0.1):
    """ Get mean regime over all experiments """
    nbExp=np.int(np.round((epsMax-epsMin)/epsInt))+1
    epsVals=np.linspace(epsMin,epsMax,nbExp)
    superpose=False ; meanBlocked=None
    for eps in epsVals:
        Avals,blockFreq=getRegime(eps0=eps,ucrit=0.05,tmin=1.e5,superpose=superpose)
        if not(isinstance(meanBlocked,np.ndarray)):
            meanBlocked=np.zeros(blockFreq.size)
        meanBlocked=meanBlocked+blockFreq
        superpose=True
    meanBlocked=meanBlocked/nbExp
    plt.plot(Avals,meanBlocked,'k',linewidth=2)

