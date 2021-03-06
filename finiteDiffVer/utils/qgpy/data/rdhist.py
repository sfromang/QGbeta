
import matplotlib.pyplot as plt
import numpy as np

from utils.utils import smooth
from utils.plt_utils import saveAndclean

class hstData:

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
                savePlt=False,saveFileName='history.png',Trot=1.,smooth=False):
        """ Plot history data """
        if not(superpose):
            plt.figure()
        plt.plot(self.time/Trot,self.umean)
        if (smooth):
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

