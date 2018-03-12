import tables
import numpy as np
import netCDF4 as ncdf

class qgData:
    def __init__(self,idump=1,filedir='./'):
        #Dump directory name
        dir=filedir+'output_'+str_suffix(idump)+'/'

        #Read first output to get basic info
        filename=dir+'qgout.000000'
        data=qgFile(filename)
        
        # Get time & size
        self.time=data.time
        self.nx=data.nx
        self.ny=data.ny

        # Define output arrays
        self.x=data.x[1:self.nx+1]
        self.y=data.y[1:self.ny+1]
        self.q=data.q
        self.psi=data.psi
        self.u=data.u
        self.v=data.v
        self.qbar  =data.qbar
        self.psibar=data.psibar

    def getVar(self,type='q',ilevel=0,ghost=False):
        """ Return 2D variable of type at ilevel """
        nx=self.nx ; ny=self.ny
        imin=1 ; imax=nx+1
        jmin=1 ; jmax=ny+1
        if (ghost):
            imin=0 ; imax=nx+2
            jmin=0 ; jmax=ny+2
        dicVar={'q'  : self.q[imin:imax,jmin:jmax,ilevel]+self.qbar[imin:imax,jmin:jmax,ilevel],\
                'psi': self.psi[imin:imax,jmin:jmax,ilevel]+self.psibar[imin:imax,jmin:jmax,ilevel],\
                'u'  : self.u[imin:imax,jmin:jmax,ilevel],\
                'v'  : self.v[imin:imax,jmin:jmax,ilevel],\
                'qbar'  : self.qbar[imin:imax,jmin:jmax,ilevel],\
                'psibar': self.psibar[imin:imax,jmin:jmax,ilevel],\
                'qprime'  : self.q[imin:imax,jmin:jmax,ilevel],\
                'psiprime': self.psi[imin:imax,jmin:jmax,ilevel],\
                'psiS': np.sum(self.psi[imin:imax,jmin:jmax,:]+self.psibar[imin:imax,jmin:jmax,:],2)}
        return dicVar[type]

class qgFile:

    def __init__(self,filename='qgout.000000',nvar=8,nbuf=3):

        #Open file
        f=open(filename,'rb')

        #Get simulations parameters from file
        self.time=get_array(f,1,'f8')
        dim=get_array(f,3,'i4')  ; [nx,ny,nlayers]=dim
        self.nx=nx ; self.ny=ny ; self.nlayers=nlayers
        params=get_array(f,2,'f8') ; [beta,Rrossby]=params
        self.beta=beta ; self.Rrossby=Rrossby

        #Get ndim and allocate arrays
        shape=(nlayers,ny+2,nx+2)
                
        #Get grid position and state variables
        ntot=(nx+2)*(ny+2)*nlayers
        self.x     =get_array(f,nx+2,'f8')
        self.y     =get_array(f,ny+2,'f8')
        self.q     =np.transpose(get_array(f,ntot,'f8').reshape(shape),(2,1,0))
        self.psi   =np.transpose(get_array(f,ntot,'f8').reshape(shape),(2,1,0))
        self.u     =np.transpose(get_array(f,ntot,'f8').reshape(shape),(2,1,0))
        self.v     =np.transpose(get_array(f,ntot,'f8').reshape(shape),(2,1,0))
        self.qbar  =np.transpose(get_array(f,ntot,'f8').reshape(shape),(2,1,0))
        self.psibar=np.transpose(get_array(f,ntot,'f8').reshape(shape),(2,1,0))

        #Close file
        f.close()

def get_array(fid,nb_count,type):
    bits_nb='i4'
    pad=np.fromfile(fid,count=1,dtype=bits_nb)
    array=np.fromfile(fid,count=nb_count,dtype=type)
    pad=np.fromfile(fid,count=1,dtype=bits_nb)
    return array

def GetSubArray3d(array,ivar,nx,ny,nz,nbuf,ndim):
    if ndim==1:
        SubArray3d=np.transpose(array[ivar,0:nz,0:ny,nbuf:nx+nbuf],(2,1,0))  
    if ndim==2:
        SubArray3d=np.transpose(array[ivar,0:nz,nbuf:ny+nbuf,nbuf:nx+nbuf],(2,1,0))  
    if ndim==3:
        SubArray3d=np.transpose(array[ivar,nbuf:nz+nbuf,nbuf:ny+nbuf,nbuf:nx+nbuf],(2,1,0))  
    
    return SubArray3d

def GetRoot(D,M,E,gamma,guess):
    #solve Q-P-E=0 with NR to get Q E=E-D Q=Q-D
    x=guess ; epsilon=1.
    while (epsilon>1.e-6):
        root=Getf(x,D,M,E,gamma)
        droot=Getfprime(x,D,M,E,gamma)
        x=x-root/droot
        epsilon=abs(root/droot/x)
    return x

def Getf(x,D,M,E,gamma):
    u2=M**2/((x+D)**2-M**2) ; lor=(1+u2)**(0.5)
    Xsi=(x-u2/(lor+1)*D)/lor**2
    p=(gamma-1.0)/gamma*Xsi
    return x-p-E

def Getfprime(x,D,M,E,gamma):
    u2=M**2/((x+D)**2-M**2) ; lor=(1+u2)**(0.5)
    dpdx=(gamma-1.)/gamma
    dpdx=dpdx*(1.+M**2/(x+D)**2*(1.-D*lor/(x+D)))
    return 1.-dpdx

def str_suffix(n,length=6):
    fewzero=''
    for i in range(length-len(str(n))):
        fewzero=fewzero+'0'
    return fewzero+str(n)

