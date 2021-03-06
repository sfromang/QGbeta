!###########################################################
!###########################################################
!###########################################################
subroutine init
  use params
  use variables
  implicit none
  integer :: i,j

  ! Read basic simulation parameters
  call read_params

  ! Allocate arrays
  call allocate_workspace

  ! Initialize computational grid
  call init_grid

  ! Initialize FFT variables
  call initFFT

  ! Initialize problem
  call condinit

  ! Restart if required
  ndump=restart
  if (ndump>0) then
     call restart_run
     call user_restart
  endif

  return
end subroutine init
!###########################################################
!###########################################################
!###########################################################
subroutine read_params
  use params
  use variables
  implicit none
  namelist /start_params/restart,tlim,verbose
  namelist /scheme_params/nx,ny,nlayers,typeBc,dt,cfl
  namelist /model_params/xmin,xmax,ymin,ymax,beta,lambda,nuH
  namelist /output_params/dtdump,dthist

  call default_params

  open(unit=1,file='input',status='old')
  read(1,start_params)
  read(1,scheme_params)
  read(1,model_params)
  read(1,output_params)
  close(1)

  if (nlayers == 2) lambdaInvSq=1.d0/lambda**2

  return
end subroutine read_params
!###########################################################
!###########################################################
!###########################################################
subroutine default_params
  use params
  use variables
  implicit none

  ! Grid defaults resolution
  nx=200
  ny=100
  nlayers=1 
  
  ! Box default size
  xmin=0.d0
  xmax=16000.d3
  ymin=-4.d6
  ymax=4.d6
  
  ! Time default parameters
  restart=0
  tlim=1.d5
  dt=0.d0
  cfl=-1.d0

  ! Output default parameters
  dtdump=1.d3
  dthist=-1.d0

  ! Model physical parameters
  beta=0.d0    ! Perform simulations on the f-plane (Oruba: 4.8e-11)
  lambda=4.5d5 ! Default Rossby deformation radius (in m)
  nuH=0.d0     ! Hyperdiffusion coefficient (none by default)

  ! Boundary conditions
  typeBC=1    !=1 (periodic); =2 (P.Williams type)

  return
end subroutine default_params
!###########################################################
!###########################################################
!###########################################################
subroutine allocate_workspace
  use params
  use variables
  implicit none

  ! Allocate grid arrays
  allocate(x(0:nx+1))
  allocate(y(0:ny+1))

  ! Allocate arrays in physical space
  allocate(q     (0:nx+1,0:ny+1,1:nlayers)) ;    q   = 0.d0
  allocate(qm1   (0:nx+1,0:ny+1,1:nlayers)) ;  qm1   = 0.d0
  allocate(qp1   (0:nx+1,0:ny+1,1:nlayers)) ;  qp1   = 0.d0
  allocate(u     (0:nx+1,0:ny+1,1:nlayers)) ;    u   = 0.d0
  allocate(v     (0:nx+1,0:ny+1,1:nlayers)) ;    v   = 0.d0
  allocate(psi   (0:nx+1,0:ny+1,1:nlayers)) ;  psi   = 0.d0
  allocate(dqdt  (0:nx+1,0:ny+1,1:nlayers)) ; dqdt   = 0.d0

  allocate(psibar(0:nx+1,0:ny+1,1:nlayers)) ;  psibar = 0.d0
  allocate(  qbar(0:nx+1,0:ny+1,1:nlayers)) ;    qbar = 0.d0

  ! Allocate arrays in spectral space
  allocate(sq     (1:nx,1:ny,1:nlayers)) ;    sq   = 0.d0
  allocate(sqm1   (1:nx,1:ny,1:nlayers)) ;  sqm1   = 0.d0
  allocate(sqm2   (1:nx,1:ny,1:nlayers)) ;  sqm2   = 0.d0
  allocate(sqp1   (1:nx,1:ny,1:nlayers)) ;  sqp1   = 0.d0
  allocate(su     (1:nx,1:ny,1:nlayers)) ;    su   = 0.d0
  allocate(sv     (1:nx,1:ny,1:nlayers)) ;    sv   = 0.d0
  allocate(spsi   (1:nx,1:ny,1:nlayers)) ;  spsi   = 0.d0
  allocate(spsim1 (1:nx,1:ny,1:nlayers)) ;  spsim1 = 0.d0
  allocate(spsim2 (1:nx,1:ny,1:nlayers)) ;  spsim2 = 0.d0
  allocate(dsqdt  (1:nx,1:ny,1:nlayers)) ; dsqdt   = 0.d0
  allocate(dsqdtm1(1:nx,1:ny,1:nlayers)) ; dsqdtm1 = 0.d0
  allocate(dsqdtm2(1:nx,1:ny,1:nlayers)) ; dsqdtm2 = 0.d0

  allocate(spsibar(1:nx,1:ny,1:nlayers)) ;  spsibar = 0.d0
  allocate(  sqbar(1:nx,1:ny,1:nlayers)) ;    sqbar = 0.d0

  allocate(kx(1:nx,1:ny))
  allocate(ky(1:nx,1:ny))
  allocate(kx2(1:nx,1:ny))
  allocate(ky2(1:nx,1:ny))

  return
end subroutine allocate_workspace
!###########################################################
!###########################################################
!###########################################################
subroutine init_grid
  use params
  use variables
  implicit none
  integer :: i,j

  dx=(xmax-xmin)/nx
  dy=(ymax-ymin)/ny
  x(0) = xmin - half*dx
  y(0) = ymin - half*dy  
  do i=1,nx+1
     x(i) = x(0) + real(i)*dx
  end do
  do j=0,ny+1
     y(j)= y(0) + real(j)*dy
  end do   

  return
end subroutine init_grid
!===============================================================================
!> Initialize arrays for FFT depending on boundary conditions
!===============================================================================
subroutine initfft
  use params
  implicit none
  integer :: i,j
  real(dp) :: scalex,scaley

  scalex=twopi/(xmax-xmin)
  scaley=twopi/(ymax-ymin)

  if (typeBC==1) then
     
     ! Computing kx & ky
     do i=1,nx/2
        kx(i,:)=(i-1)
     end do
     do i=nx/2+1,nx
        kx(i,:)=(i-nx-1)
     end do
     kx = scalex * kx

     do j=1,ny/2
        ky(:,j)=(j-1)
     end do
     do j=ny/2+1,ny
        ky(:,j)=(j-ny-1)
     end do
     ky = scaley * ky

  endif

  if (typeBC==2) then
     
     ! Computing kx2 & ky2
     do i=1,nx/2
        kx2(i,:)=+(i-1)*(i-1)
     end do
     do i=nx/2+1,nx
        kx2(i,:)=+(i-nx-1)*(i-nx-1)
     end do
     kx2 = scalex**2 * kx2

     ky2(1:nx,1)=0.d0
     do j=2,ny
        ky2(1:nx,j)=+(j-1)*(j-1)/4.
     end do
     ky2 = scaley**2 * ky2

  endif

  return
end subroutine initfft
