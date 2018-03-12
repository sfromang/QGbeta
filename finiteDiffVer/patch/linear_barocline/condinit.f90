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
  real(dp) :: amp,maxPsibar,kx,xi
  real(dp) :: Utop,Lx
  logical :: randomInit
  namelist /init_params/Utop,nuH,Lx,randomInit

  if (verbose) write (*,*) "Entering condinit.f90..."

  ! Default parameters for the problem
  Utop=30.d0   ! Upper layer velocity velocity
  Lx=4.e6      ! Wavelength of perturbation
  nuH=5.d15    ! Hyperdiffusion coefficient
  randomInit=.false. ! True if initial conditions consists in random noise

  ! Read problem input data
  open(unit=1,file='input',status='old')
  read(1,init_params)
  close(1)

  ! Compute variables useful for eigenmodes computations
  kx=2.*pi/Lx
  xi=(2.d0*LambdaInvSq-kx**2)/(2.d0*LambdaInvSq+kx**2)
  
  ! Compute initial steady state psibar and PV qbar (including ghost zones)
  do j=0,ny+1
     psibar(:,j,1)=-Utop*y(j)
  end do
  qbar(:,:,1)=-(psibar(:,:,1)-psibar(:,:,2))*LambdaInvSq        
  qbar(:,:,2)=+(psibar(:,:,1)-psibar(:,:,2))*LambdaInvSq

  ! Perturbed streamfunction
  iseed=0 ; amp=1.e-10 ; maxPsibar=maxval(psibar)
  psi(:,:,1)=0.d0 ; psi(:,:,2)=0.d0
  if (randomInit) then
     ! Random perturbations on perturbed streamfunction
     do i=0,nx+1
        do ilayer=1,2
           call ran2(iseed,rvalue)
           psi(i,j,ilayer)=maxPsibar*amp*(rvalue-0.5)
        end do
     end do
  else
     ! Sinusoidal perturbations on perturbed streamfunction
     do i=0,nx+1
        psi(i,:,1)=cos(kx*x(i))-xi**0.5*sin(kx*x(i))
        psi(i,:,2)=cos(kx*x(i))+xi**0.5*sin(kx*x(i))
     end do
     psi=maxPsibar*amp*psi
  endif

  ! Compute PV
  call getQ(psi,q) ; call computeBC(q,nx,ny,nlayers)

  ! Compute other variables (qm1,u,v) for saving purposes.
  qm1=q
  call getVel(psi+psibar,u,v)
  call computeBC(u,nx,ny,nlayers)
  call computeBC(v,nx,ny,nlayers)

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
