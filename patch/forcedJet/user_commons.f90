module user_params
  use precision
  implicit none

  real(dp) :: Umax,sigma                          !< Jet velocity & jet width
  real(dp), dimension(:,:), allocatable :: psiJet !< Forced streamfunction

  real(dp) :: tau,tauM                            !< Forcing & frictionnal timescales
  real(dp) :: nuH                                 !< Hyperdiffusion coefficient

end module user_params
