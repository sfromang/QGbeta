!###########################################################
!###########################################################
!###########################################################
subroutine dissip(spsi,dsqdt)
  use params
  use variables , only : spsibar
  implicit none

  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(in ) :: spsi
  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(out) :: dsqdt
  complex(dp), dimension(1:nx,1:ny,1:nlayers) :: sdissip

  ! 6th order small scale hyperdiffusion
  sdissip=spsi-spsibar
  call delSq(sdissip,nx,ny,nlayers)
  call delSq(sdissip,nx,ny,nlayers)
  sdissip=nuH*sdissip

  dsqdt = dsqdt + sdissip

  return
end subroutine dissip
