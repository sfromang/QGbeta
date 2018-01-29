import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize

from qgpy.plot import qghist

def histoFit(umean,doPlt=True):
    """ Fit histogram of umean frequencies by a double Gaussian """
    hist,binEdges=np.histogram(umean,bins=30,normed=True)
    bins=[ (binEdges[i]+binEdges[i+1])/2. for i in range(binEdges.size-1) ]
    # First try single Gaussian fit
    pCoeff,pCov=optimize.curve_fit(fOneGaussian,bins,hist,p0=[50.,umean.mean()**0.5,7.e-3])
    A0I=pCoeff[0] ; A1I=0. ; x0I=pCoeff[1] ; x1I=0. ; dxI=pCoeff[2]        
    R2I=getRsq(bins,hist,A0I,A1I,x0I,x1I,dxI)
    # Perform fit with two Gaussians
    pCoeff,pCov=optimize.curve_fit(f,bins,hist,p0=[50.,50.,0.2,0.3,0.005])
    if (pCoeff[2]<pCoeff[3]):
        A0=pCoeff[0] ; A1=pCoeff[1] ; dx=pCoeff[4]
        x0=pCoeff[2]
        x1=pCoeff[3]
    else:
        A1=pCoeff[0] ; A0=pCoeff[1] ; dx=pCoeff[4]
        x1=pCoeff[2]
        x0=pCoeff[3]
    R2=getRsq(bins,hist,A0,A1,x0,x1,dx)
    if ((R2>0.) and (R2>R2I) and (np.abs(R2-R2I)>1.e-4)):
        print '   Best fit with two Gaussians: ',R2                
    else:
        A0=A0I ; A1=A1I ; x0=x0I ; x1=x1I ; dx=dxI
        print '   Best fit with one Gaussian: ',R2I
    if (doPlt):
        xvals=np.linspace(0.,np.max(bins),100.)
        #qghist.plotHisto(umean)
        plt.plot(xvals,[f(x,A0,A1,x0,x1,dx) for x in xvals],'r',linewidth=2)
    return A0,A1,x0**2,x1**2,np.abs(dx)

def f(x,A0,A1,x0,x1,dx):
    """ Function to fit p=[A0,A1,x0,x1,dx] """
    f0=A0*np.exp(-(x-x0**2)**2/2./dx**2)
    f1=A1*np.exp(-(x-x1**2)**2/2./dx**2)
    return f0+f1

def fOneGaussian(x,A0,x0,dx):
    """ Function to fit p=[A0,A1,x0,x1,dx] """
    return A0*np.exp(-(x-x0**2)**2/2./dx**2)

def getRsq(x,y,A0,A1,x0,x1,dx):
    """ Compute R^2 fit goodness """
    return 1-RSS(x,y,A0,A1,x0,x1,dx)/TSS(x,y)

def TSS(x,y):
    """ Compute variance of the data """
    return np.sum((y-y.mean())**2)

def RSS(x,y,A0,A1,x0,x1,dx):
    """ Compute residual variance of the data """
    return np.sum([(y[i]-f(x[i],A0,A1,x0,x1,dx))**2 for i in range(y.size)])
