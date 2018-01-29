!===============================================================================
!> Precision module; define single and double precisions
!===============================================================================
module precision

  integer, parameter :: sp=kind(1.0E0) !< single precision
!#ifndef NPRE                            
!  integer, parameter :: dp=kind(1.0D0) !< default double precision
!#else                                   
!#if NPRE==4                            
!  integer, parameter :: dp=kind(1.0E0) !< double precision if NPRE==4
!#else                                   
  integer, parameter :: dp=kind(1.0D0) !< double precision is NPRE!=4
!#endif
!#endif

end module precision
!===============================================================================
!> Constant module; define constants of the code
!===============================================================================
module const
  use precision

  real(dp) :: bigreal = 1.0d+30               !< large real number
  real(dp) :: zero    = 0.0d0                 !< 0
  real(dp) :: one     = 1.0d0                 !< 1
  real(dp) :: two     = 2.0d0                 !< 2
  real(dp) :: three   = 3.0d0                 !< 3
  real(dp) :: four    = 4.0d0                 !< 4
  real(dp) :: two3rd  = 0.6666666666666667d0  !< 2/3
  real(dp) :: half    = 0.5d0                 !< 1/2
  real(dp) :: third   = 0.33333333333333333d0 !< 1/3
  real(dp) :: forth   = 0.25d0                !< 1/4
  real(dp) :: sixth   = 0.16666666666666667d0 !< 1/6
  real(dp) :: pi      = 2.d0*asin(1.d0)       !< \f$ \pi \f$
  real(dp) :: twopi   = 4.d0*asin(1.d0)       !< \f$ 2 \pi \f$

end module const
!===============================================================================
!> Parameters modules; parameters of the problem
!===============================================================================
module params
  use precision
  use const

  ! Physical parameters
  real(dp) :: beta=1.d0                !< beta parameter of the model
  real(dp) :: lambda=-1.d0             !< Rossby deformation radius (=-1 if one layer)
  real(dp) :: LambdaInvSq=+1.d0        !< One over Rossby deformation radius

  ! Boundary conditions parameters
  integer :: typeBC=1                  !< Boundary conditions type (see init.f90)

  ! Start parameters
  integer  :: restart         !< output index for restart
  real(dp) :: tlim=1000.      !< limit time
  logical  :: verbose=.false. !< verbosity
  logical  :: debug=.false.   !< turn-on debugging features

  ! Model parameters
  integer  :: nx=1           !< number of cells in a subdomain, x-direction
  integer  :: ny=1           !< number of cells in a subdomain, y-direction
  integer  :: nlayers=1      !< number of layers in z-direction
  integer  :: nxglob=1       !< total number of cells, x-direction
  integer  :: nyglob=1       !< total number of cells, y-direction
  real(dp) :: nu=0.d0        !< viscosity
  real(dp) :: xmin           !< minimum x coordinate
  real(dp) :: xmax           !< maximum x coordinate
  real(dp) :: ymin           !< minimum y coordinate
  real(dp) :: ymax           !< maximum y coordinate
  real(dp) :: dtvalue        !< set dt to a specific value

  ! Mesh and wavenumbers parameters
  real(dp) :: dx !< resolution in x-direction
  real(dp) :: dy !< resolution in y-direction
  real(dp), dimension(:,:), allocatable :: kx  !< wavenumbers in x-direction
  real(dp), dimension(:,:), allocatable :: ky  !< wavenumbers in y-direction
  real(dp), dimension(:,:), allocatable :: kx2 !< wavenumbers**2 in x-direction
  real(dp), dimension(:,:), allocatable :: ky2 !< wavenumbers**2 in y-direction

  ! Output parameters
  real(dp) :: dtdump           !< elapsed time between two outputs
  real(dp) :: dthist           !< elapsed time between two history
  character(LEN=10) :: io_type !< output format type

end module params
!===============================================================================
!> Variables module; variables of the problem
!===============================================================================
module variables
  use precision

  real(dp), dimension(:,:,:), allocatable   :: q,qm1,qp1 !< vorticity array
  real(dp), dimension(:,:,:), allocatable   :: psi       !< streamfunction array
  real(dp), dimension(:,:,:), allocatable   :: u,v       !< velocity array
  real(dp), dimension(:,:,:), allocatable   :: dqdt      !< time-derivative array
  real(dp), dimension(:    ), allocatable   :: x         !< x-coordinates array
  real(dp), dimension(:    ), allocatable   :: y         !< y-coordinates array

  complex(dp), dimension(:,:,:), allocatable   :: sq,sqm1,sqm2,sqp1  !< PV array in spectral space
  complex(dp), dimension(:,:,:), allocatable   :: spsi,spsim1,spsim2 !< streamfunction array in spectral space
  complex(dp), dimension(:,:,:), allocatable   :: su,sv              !< velocity array in spectral space
  complex(dp), dimension(:,:,:), allocatable   :: dsqdt              !< time-derivative array in spectr al space 
  complex(dp), dimension(:,:,:), allocatable   :: dsqdtm1,dsqdtm2    !< time-derivative array in spectral space

  !< analytical streamfunction and PV
  real(dp), dimension(:,:,:), allocatable   :: psibar, qbar
  complex(dp), dimension(:,:,:), allocatable   :: spsibar, sqbar

  !< Time variables
  real(dp) :: dt   !< timestep
  real(dp) :: time !< time

  !< Output variables
  integer :: ndump !< # of the current output

end module variables
!###########################################################
!###########################################################
!###########################################################
!###########################################################
module fftVar
  use precision

  integer :: inc,ierr
  integer :: lensav,lensavx,lensavyN,lensavyD
  integer :: lenwrk,lenwrkx,lenwrkyN,lenwrkyD

  real(dp), allocatable, dimension(:) :: wsave,wsavex,wsaveyN,wsaveyD
  real(dp), allocatable, dimension(:) :: work,workx,workyN,workyD

end module fftVar
