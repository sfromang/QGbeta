!###########################################################
!###########################################################
!###########################################################
subroutine condinit
  use user_params
  use params
  use variables
  implicit none
  integer :: i,j,iseed,ilayer
  real rvalue
  real(dp) :: amp
  namelist /init_params/Umax,sigma,tau,tauM,nuH

  if (verbose) write (*,*) "Entering condinit.f90..."

  ! Default parameters for the problem
  ! Default parameters are based on Zurita-Gotor et al. (2014)
  Umax=40.d0         ! Jet max velocity
  sigma=2.5e6        ! Jet width
  tau=20.d0*86400.d0 ! Forcing time
  tauM=3.d0*86400.d0 ! Dissipation time
  nuH=5.d15          ! Hyperdiffusion coefficient

  ! Read problem input data
  open(unit=1,file='input',status='old')
  read(1,init_params)
  close(1)

  ! Compute forcing
  do j=0,ny+1
     psibar(:,j,1)=-Umax*sqrt(pi)/2.d0*sigma*erf(y(j)/sigma)
  end do
  call getQ(psibar,qbar) ; call computeBCzeroGrad(qbar,nx,ny,nlayers)

  iseed=0 ; amp=1.e-10

  ! Background state
  psi(:,:,1)=0.d0 ; psi(:,:,2)=0.d0
  do j=0,ny+1
     do i=0,nx+1
        do ilayer=1,nlayers
           call ran2(iseed,rvalue)
           psi(i,j,ilayer)=psibar(i,j,ilayer)+psibar(i,j,ilayer)*amp*(rvalue-0.5)
        end do
     end do
  end do
  !call computeBCzeroGrad(psi,nx,ny,nlayers)

  ! Compute PV
  !call getQ(psi,q) ; call computeBCzeroGrad(q,nx,ny,nlayers)

  ! Compute FFT and store variable in spectral space
  call var2svar(psibar,spsibar,nx,ny,nlayers)
  call var2svar(q,sq,nx,ny,nlayers)
  call var2svar(psi,spsi,nx,ny,nlayers)
  call svar2var(spsi,psi,nx,ny,nlayers)
  call spsi2sq(spsi,sq,nx,ny,nlayers,LambdaInvSq)
  call svar2var(sq,q,nx,ny,nlayers)
  call computeBCzeroGrad(q,nx,ny,nlayers)
  call computeBCzeroGrad(psi,nx,ny,nlayers)

  ! Compute other variables (qm1,u,v) for saving purposes.
  sqm1=sq
  call spsi2su(spsi,su,nx,ny,nlayers,2.d0*pi/(ymax-ymin)) ! Compute u in spectral space
  call su2u(su,u,nx,ny,nlayers)                           ! Transform u back to real space
  call spsi2sv(spsi,sv,nx,ny,nlayers,2.d0*pi/(xmax-xmin)) ! Compute v in spectral space
  call svar2var(sv,v,nx,ny,nlayers)                       ! Transform v back to real space

  return
end subroutine condinit
!====================================================================
!  numerical recipes random number generator ran2
!    requires input seed value=iseed
!    returns real random number=rvalue
!    Also updates iseed for next call 
!
      subroutine ran2(iseed,rvalue)
     
      integer iseed
      real rvalue
      integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real AM,EPS,RNMX
      parameter (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,&
               & IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,&
               & IR2=3791,NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      integer idum2,jj,kk,iv(NTAB),iy
      data idum2/123456789/, iv/NTAB*0/, iy/0/
!
      idum=iseed
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 jj=NTAB+8,1,-1
          kk=idum/IQ1
          idum=IA1*(idum-kk*IQ1)-kk*IR1
          if (idum.lt.0) idum=idum+IM1
          if (jj.le.NTAB) iv(jj)=idum
11      continue
        iy=iv(1)
      endif
      kk=idum/IQ1
      idum=IA1*(idum-kk*IQ1)-kk*IR1
      if (idum.lt.0) idum=idum+IM1
      kk=idum2/IQ2
      idum2=IA2*(idum2-kk*IQ2)-kk*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      jj=1+iy/NDIV
      iy=iv(jj)-idum2
      iv(jj)=idum
      if(iy.lt.1)iy=iy+IMM1
      rvalue=min(AM*iy,RNMX)
      iseed=idum
      return
      end
