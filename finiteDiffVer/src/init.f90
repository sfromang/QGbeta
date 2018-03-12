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
  namelist /scheme_params/nx,ny,nlayers,typeBc,dtvalue
  namelist /model_params/xmin,xmax,ymin,ymax,beta,lambda
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
  dt=2.d2
  dtvalue=0.d0

  ! Output default parameters
  dtdump=1.d3
  dthist=-1.d0

  ! Model physical parameters
  beta=0.d0    ! Perform simulations on the f-plane (Oruba: 4.8e-11)
  lambda=4.5d5 ! Default Rossby deformation radius (in m)

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

  allocate(x(0:nx+1))
  allocate(y(0:ny+1))
  allocate(q     (0:nx+1,0:ny+1,1:nlayers)) ;    q   =0.d0
  allocate(qm1   (0:nx+1,0:ny+1,1:nlayers)) ;  qm1   =0.d0
  allocate(qp1   (0:nx+1,0:ny+1,1:nlayers)) ;  qp1   =0.d0
  allocate(u     (0:nx+1,0:ny+1,1:nlayers)) ;    u   =0.d0
  allocate(v     (0:nx+1,0:ny+1,1:nlayers)) ;    v   =0.d0
  allocate(psi   (0:nx+1,0:ny+1,1:nlayers)) ;  psi   =0.d0
  allocate(dqdt  (0:nx+1,0:ny+1,1:nlayers)) ; dqdt   =0.d0

  allocate(psibar(0:nx+1,0:ny+1,1:nlayers)) ;  psibar=0.d0
  allocate(  qbar(0:nx+1,0:ny+1,1:nlayers)) ;    qbar=0.d0

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
