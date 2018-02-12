!###########################################################
!###########################################################
!###########################################################
subroutine dissip(spsi,dsqdt)
  use params
  implicit none

  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(in ) :: spsi
  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(out) :: dsqdt
  complex(dp), dimension(1:nx,1:ny,1:nlayers) :: sdissip

  ! 6th order small scale hyperdiffusion
  sdissip=spsi
  call delSq(sdissip,nx,ny,nlayers)
  call delSq(sdissip,nx,ny,nlayers)
  call delSq(sdissip,nx,ny,nlayers)
  sdissip=nuH*sdissip

  dsqdt = dsqdt - sdissip

  return
end subroutine dissip
