module user_params
  use params

  real(dp) :: Umax,sigma                          !< Jet velocity & jet width
  real(dp), dimension(:,:), allocatable :: psiJet !< Forced streamfunction

  !Dissipation variables
  real(dp) :: tau,tauM                            !< Forcing & frictionnal timescales
  real(dp) :: nuH                                 !< Hyperdiffusion coefficient

  !Forcing variables
  character(len=10) :: forcingType                !< Forcing type ('tian' or 'charney')
  logical :: varForcing                           !< Forcing is constant (==.false.) or variable
  real(dp) :: Amin,Amax,Ainc                      !< Forcing variations variables
  real(dp) :: A                                   !< Forcing amplitude
  real(dp) :: kappa                               !< Forcing dissipation "timescale"


  !Topographic variables
  character(len=10) :: topoType                   !< Topography type ('tian' or 'charney')
  real(dp) :: h0                                  !< Topography maximum height
  real(dp), dimension(:,:), allocatable :: hB     !< Topography at each cells
  

end module user_params

module noise_params
  use precision
  
  logical :: varNoise                             !< Variable noise if True, cst noise amplitude otherwise
  logical :: redTime,redSpace                     !< noise correlated in space (redSpace) and time (redTime)
  real(dp) :: noiseAmp                            !< Noise amplitude
  real(dp) :: noiseMinAmp,noiseMaxAmp             !< Noise minimum and maximum amplitude (in case it varies)
  real(dp) :: noiseInc                            !< noise relative incremental steps (in %)

  real(dp) :: w0,dw0                                       !< noise spectral caracteristics variables

  real(dp), dimension(:,:), allocatable :: noise,savNoise

end module noise_params
