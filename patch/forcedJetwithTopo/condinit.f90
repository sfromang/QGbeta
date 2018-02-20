!###########################################################
!###########################################################
!###########################################################
subroutine condinit
  use user_params
  use noise_params
  use params
  use variables
  implicit none
  integer :: i,j,iseed,ilayer
  real(dp) :: rvalue
  real(dp) :: amp,f1,f2,y0,y1,factor,scaling
  namelist /init_params/forcingType,varForcing,A,Amin,Amax,Ainc,kappa,\
                        topoType,h0,amp,noiseAmp,tauM,w0,dw0,\
                        varNoise,noiseMinAmp,noiseMaxAmp,noiseInc,\
                        redTime,redSpace,scaling

  if (verbose) write (*,*) "Entering condinit.f90..."

  ! Default parameters for the problem
  varForcing=.false. ! cst vs. varying amplitude of the large scale flow forcing
  forcingType='tian'
  A=2.d0
  Amin=2.d0
  Amax=2.d0
  Ainc=1.d0
  kappa=5.d-3
  !nuH=5.d-5         ! Hyperdiffusion coefficient
  topoType='tianJFM'
  tauM=1.d1         ! Bottom layer friction
  h0=0.1d0          ! Bottom topography height
  ! Noise default parameters
  amp=0.d0           ! Initial random noise amplitude on streamfunction
  noiseAmp=0.d0      ! Amplitude of noise forcing in vorticity tendency equation
  varNoise=.false.   ! cst vs. varying amplitude of the noise amplitude
  noiseMinAmp=0.d0   ! noise minimum amplitude in case of time varying noise
  noiseMaxAmp=0.d0   ! noise maximum amplitude in case of time varying noise
  noiseInc=0.d0      ! noise amplitude increment
  w0=0.1d0           ! Wavenumber of the forcing (in 1D, the period of the noise is Nx*|0.5-w0| cells)
  dw0=0.01d0         ! Half-width of the forcing in k-space (same unit as w0)
  redTime=.false.    ! Noise correlated in time
  redSpace=.false.   ! Noise correlated in space
  scaling=-1.d0      ! Scaling factor of the reference PV profile obtained from equilibrium states

  ! Read problem input data
  open(unit=1,file='input',status='old')
  read(1,init_params)
  close(1)

  ! Compute bottom layer topography
  call computeTopography

  ! Compute vorticity forcing
  call getForcingAmp(time)
  call computeForcing

  ! Initial flow
  call initRandom
  psi(:,:,1:nlayers)=0.d0
  iseed=0
  do j=1,ny
     do i=1,nx
        do ilayer=1,nlayers
           call random_number(rvalue)
           psi(i,j,ilayer)=psibar(i,j,ilayer)+psibar(i,j,ilayer)*amp*(rvalue-0.5)
        end do
     end do
  end do
  ! Streamfunction BC in x and y
  psi(0,:,1)=psi(nx,:,1) ; psi(nx+1,:,1)=psi(1,:,1)
  psi(:,0,1)=sum(psi(1:nx,1,1))/nx ; psi(:,ny+1,1)=sum(psi(1:nx,ny,1))/nx

  ! Compute PV
  call getQ(psi,q) ; call computeBC(q,nx,ny,nlayers)
  q=qbar
  ! do j=1,ny
  !    do i=1,nx
  !       do ilayer=1,nlayers
  !          call random_number(rvalue)
  !          qbar(i,j,ilayer)=cos(2.*twopi/(xmax-xmin)*x(i))*cos(2*pi/(ymax-ymin)*y(j))
  !       end do
  !    end do
  ! end do
  ! call var2svar(qbar,sqbar,nx,ny,nlayers)
  ! do j=1,ny
  !    write (*,*) j,kx2(2,j),ky2(2,j),sqbar(2,j,1)
  ! end do
  ! call sq2spsi(sqbar,spsibar,nx,ny,nlayers)
  ! call spsi2sq(spsibar,sqbar,nx,ny,nlayers)
  ! call svar2var(sqbar,qbar,nx,ny,nlayers)
  ! call svar2var(spsibar,psibar,nx,ny,nlayers)
  ! q=qbar ; psi=psibar

  !sqbar=0.d0
  !sqbar(2,2,1)=cmplx(1.d0,0.d0)
  !call svar2var(sqbar,qbar,nx,ny,nlayers) ; call computeBC(q,nx,ny,nlayers)

  
  ! Initial flow from existing equilibrium state (stored in files zonal.bin and blocked.bin)
  if (scaling>=0) then
     call computeInitFlow(scaling)
     call computeBC(q,nx,ny,nlayers)
  endif

  ! Compute other variables (qm1,u,v) for saving purposes.
  qm1=q
  !call getVel(psi+psibar,u,v)
  call getVel(psi,u,v)

  ! Compute red noise
  !call getRedNoise(noise)
  !u(1:nx,1:ny,1)=noise
  !call output
  !stop

  ! Compute FFT and store variable in spectral space
  !call var2svar(psibar,spsibar,nx,ny,nlayers)
  call var2svar(q,sq,nx,ny,nlayers)
  call var2svar(psi,spsi,nx,ny,nlayers)
  !call svar2var(spsi,psi,nx,ny,nlayers)
  !call spsi2sq(spsi,sq,nx,ny,nlayers,LambdaInvSq)
  !call svar2var(sq,q,nx,ny,nlayers)
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
!===============================================================================
!> Compute forcing amplitude as a function of time
!===============================================================================
subroutine getForcingAmp(time)
  use params , only : tlim
  use user_params , only : varForcing,A,Amin,Amax,Ainc
  use precision
  implicit none
  real(dp), intent(in) :: time
  integer :: Nstep,tstep 

  if (varForcing) then
     if (Amax>Amin) then
        Nstep=int((Amax-Amin)/Ainc)
        tstep=int(tlim/Nstep)
        A=Amin+Ainc*int(time/tstep)
     else
        Nstep=int((Amin-Amax)/Ainc)
        tstep=int(tlim/Nstep)
        A=Amin-Ainc*int(time/tstep)
     endif
  endif

  write (*,*) "         Forcing amplitude: ",A
  
  return
end subroutine getForcingAmp
!===============================================================================
!> Compute initial condition for the flow (q & psi) 
!===============================================================================
subroutine computeInitFlow(scaling)
  use params
  use variables
  real(dp), intent(in) :: scaling
  real(dp) :: dummy
  real(dp), dimension(0:nx+1,0:ny+1,1:nlayers) :: qzonal,qblocked !< PV array of zonal & blocked flows

  ! Read PV of zonal flow
  open(unit=10,file='../zonal.bin',status='old',form='unformatted')
  read(10) dummy
  read(10) nx,ny,nlayers
  read(10) beta,lambda
  read(10) x
  read(10) y
  read(10) qzonal
  close(10)

  ! Read PV of blocked flow
  open(unit=10,file='../blocked.bin',status='old',form='unformatted')
  read(10) dummy
  read(10) nx,ny,nlayers
  read(10) beta,lambda
  read(10) x
  read(10) y
  read(10) qblocked
  close(10)
  
  ! Compute PV & streamfunction
  ! do i=0,nx+1
  !    q(i,:,:)=scaling*half*sum((qzonal+qblocked),dim=1)/nx
  ! end do
  q=scaling*qzonal+(1.d0-scaling)*qblocked
  call poissonfft(q(:,:,1),psi(:,:,1),.false.,.false.,.false.,.true.) 
  psi(0,:,1)=psi(nx,:,1) ; psi(nx+1,:,1)=psi(1,:,1)
  psi(:,0,1)=sum(psi(1:nx,1,1))/nx ; psi(:,ny+1,1)=sum(psi(1:nx,ny,1))/nx
  
  return
end subroutine computeInitFlow
!===============================================================================
!> Compute forcing term (qbar & psibar)
!===============================================================================
subroutine computeForcing
  use params
  use variables , only : qbar,psibar,sqbar,spsibar,y
  use user_params , only : forcingType,A
  implicit none
  real(dp) :: y0,y1,factor,f1,f2
  integer :: j

  ! Compute forcing
  if (forcingType=='tian') then
     do j=0,ny+1
        y0=pi/32.d0 ; y1=pi/4.d0 ; factor=ymax/(pi/2.d0)
        y0=y0*factor ; y1=y1*factor
        f1=-(y(j)-y1)**2/y0**2
        f2=-(y(j)+y1)**2/y0**2
        qbar(0:nx+1,j,1)=A*(exp(f1)-7.d0/13.d0*exp(f2))
     end do
  endif
  if (forcingType=='charney') then
     do j=0,ny+1
        qbar(0:nx+1,j,1)=A*sqrt(2.d0)*sin(y(j))
     end do
  endif

  ! Copy same mean state in both layers in case of two layers
  if (nlayers==2) then
     qbar  (:,:,2)=  qbar(:,:,1)
     psibar(:,:,2)=psibar(:,:,1)
  endif

  ! Dirichlet lateral BC (with psi=0)
  call var2svar(qbar,sqbar,nx,ny,nlayers)
  call dealiasing(sqbar,nx,ny,nlayers)
  call svar2var(sqbar,qbar,nx,ny,nlayers)
  call sq2spsi(sqbar,spsibar,nx,ny,nlayers)
  call dealiasing(spsibar,nx,ny,nlayers)
  call svar2var(spsibar,psibar,nx,ny,nlayers)

  return
end subroutine computeForcing
!===============================================================================
!> Compute bottom layer topography
!===============================================================================
subroutine computeTopography
  use params
  use variables , only : x,y,qbar,sqbar
  use user_params , only : topoType,h0,hB,dhBdx
  implicit none
  real(dp) :: x0,x1,dxB0,hB1,hB2,dhB1dx,dhB2dx
  integer :: i,j

  if (.not.(allocated( hB)))   allocate( hB(0:nx+1,0:ny+1)) ;    hB=0.d0
  if (.not.(allocated(dhBdx))) allocate(dhBdx(0:nx+1,0:ny+1)) ; dhBdx=0.d0

  if ((topoType=='tianJFM').or.(topoType=='tianSeb')) then
     if (topoType=='tianJFM') dxB0=0.968d0 !JFM paper setup
     if (topoType=='tianSeb') dxB0=0.36d0  !JFM CORRECTED paper setup
     x0=xmax/4.d0 ; x1=3.d0*xmax/4.d0
     do i=0,nx+1
        hB1=h0*exp(-(x(i)-x0)**2/dxB0**2) !h0*exp(-(x(i)-x0)**2/2.d0/dxB0**2)
        hB2=h0*exp(-(x(i)-x1)**2/dxB0**2) !exp(-(x(i)-x1)**2/2.d0/dxB0**2)
        dhB1dx=-2.d0*(x(i)-x0)/dxB0**2*hB1
        dhB2dx=-2.d0*(x(i)-x1)/dxB0**2*hB2
         hB  (i,:)= hB1  + hB2
        dhBdx(i,:)=dhB1dx+dhB2dx
     end do
  endif
  if (topoType=='charney') then
     do j=0,ny+1
        do i=0,nx+1
           hB(i,j)=0.5*h0*(1.d0-cos(4.d0*pi/(xmax-xmin)*x(i)))*cos(pi*y(j)/(ymax-ymin))
        end do
     end do
  endif

  return
end subroutine computeTopography
!===============================================================================
!> Initialize seed for Random Number Generator
!===============================================================================
subroutine initRandom
  implicit none
  integer, dimension(8) :: vtime
  integer, dimension(:), allocatable :: seed
  integer :: size
  call date_and_time(values=vtime)
  call random_seed(size=size)
  allocate(seed(size))
  seed = (vtime(4)*(360000*vtime(5) + 6000*vtime(6) + 100*vtime(7) + vtime(8)))
  call random_seed(put=seed)
  return
end subroutine initRandom
