!###########################################################
!> file poissonfft.f90
!! \brief
!! \b QGbeta:
!! This is poissonfft subroutines. It 
!! \details
!! Contains poissonfft(), initFFT()
!! \author
!! Sébastien Fromang
!! \date
!! \b created:          18-06-2015 
!! \b last \b modified: 18-07-2015
!===============================================================================
!> Solve Poisson equation
!> if "withsource" is False, solve:
!> del^2 psi = q(i,j)
!> if "withsource" is True, solve:
!> del^2 psi = q(i,j) + 2/lambda^2*psi(i,j) 
!===============================================================================
subroutine poissonfft(sq,spsi,withSource)
  use params
  implicit none

  logical , intent(in) :: withSource
  complex(dp), dimension(1:nx,1:ny), intent(in ) :: sq
  complex(dp), dimension(1:nx,1:ny), intent(out) :: spsi
  real(dp) :: factor
  integer :: i,j

  factor=0.d0
  if (withSource) factor=lambdaInvSq
  spsi(1,1)=0.d0
  do j=2,ny
     do i=2,nx
        spsi(i,j)=-sq(i,j)/(kx2(i,j)+ky2(i,j)+2.d0*factor)
     end do
  end do
  do i=2,nx
     spsi(i,1)=-sq(i,1)/(kx2(i,1)+ky2(i,1)+2.d0*factor)
  end do
  do j=2,ny
     spsi(1,j)=-sq(1,j)/(kx2(1,j)+ky2(1,j)+2.d0*factor)
  end do

  return
end subroutine poissonfft
