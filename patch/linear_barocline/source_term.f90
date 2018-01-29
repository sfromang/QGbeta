!###########################################################
!###########################################################
!###########################################################
subroutine source_term
  use user_params
  use params
  use variables
  implicit none
  real(dp), dimension(0:nx+1,0:ny+1) :: xi
  real(dp), dimension(0:nx+1,0:ny+1,nlayers) :: dtmp,dterm
  real(dp) :: psiR
  integer :: ilayer

  ! Get vorticity
  call getVorticity(psi(:,:,2),xi)

  ! 6th order hyperdiffusion
  call getDelSq(psi,dterm)
  call computeBC(dterm,nx,ny,nlayers)
  call getDelSq(dterm,dtmp)
  call computeBC(dtmp,nx,ny,nlayers)
  call getDelSq(dtmp,dterm)
  call computeBC(dterm,nx,ny,nlayers)
  dterm=nuH*dterm

  ! Small scale hyperdiffusion
  dqdt=dqdt-dterm

  return
end subroutine source_term
