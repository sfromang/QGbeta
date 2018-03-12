!###########################################################
!###########################################################
!###########################################################
subroutine user_restart
  use user_params
  use variables
  implicit none

  ! Compute bottom layer topography
  call computeTopography

  ! Compute vorticity forcing
  call getForcingAmp(time)
  call computeForcing

  return
end subroutine user_restart
