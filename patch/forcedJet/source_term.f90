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
  real(dp) :: psiR,scalex,scaley
  integer :: ilayer
  complex(dp), dimension(nx,ny,nlayers) :: sdissip,sxi

  ! Jet forcing term
  do ilayer=1,nlayers
     dsqdt(:,:,ilayer)=dsqdt(:,:,ilayer)-\
                      (-1.d0)**ilayer/tau*(spsi(:,:,1)-spsi(:,:,2)-spsibar(:,:,1))*lambdaInvSq
  end do

  ! Get vorticity
  !sdissip=spsim1 !This is for leapfrog scheme
  sdissip=spsi
  call delSq(sdissip,nx,ny,nlayers)
  sxi=sdissip ; sxi(:,:,1)=0.d0

  ! Add bottom layer dissipation to PV tendency
  dsqdt = dsqdt - sxi/tauM

  return
end subroutine source_term
