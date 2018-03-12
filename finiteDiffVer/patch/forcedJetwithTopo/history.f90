!===============================================================================
!> \file history.f90
!! \brief
!! \b QGbeta:
!! This is history subroutine.
!! \details
!! Contains history()
!! \author
!! Sébastien Fromang <sebastien.fromang@cea.fr>
!! \copyright
!! Copyrights 2015, CEA.
!! This file is distributed under the CeCILL-A & GNU/GPL licenses, see
!! <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html> and
!! <http://www.gnu.org/licenses/>
!! \date
!! \b created:          09-30-2015 
!! \b last \b modified: 09-30-2015
!<
!===============================================================================
!> \brief
!! This routine produces history file of the simulation.
!===============================================================================
subroutine history
  use params
  use variables
  implicit none
  character(LEN=80) :: filename="history.txt", hist_format
  logical :: fexist
  real(dp) :: xprobe,yprobe
  real(dp), dimension(nlayers) :: uprobe,umean
  integer :: iprobe,jprobe
  ! Mask variables
  integer :: cellCount,i,j
  real(dp) :: xmaskMin,xmaskMax,ymaskMin,ymaskMax
  real(dp), dimension(nlayers) :: umask

  if (verbose) print*, "> Writing history..."

  ! Compute zonal velocity at probe location and averaged over the channel
  ! xprobe=10.116d0 ; yprobe=0.d0
  !xprobe=3.71d0 ; yprobe=0.d0
  !xprobe=3.93d0 ; yprobe=0.d0
  xprobe=2.25d0 ; yprobe=0.d0 ! with MySetup (based on corrected JFM paper)
  xprobe=7.d0 ; yprobe=0.d0 ! with JFM paper setup
  iprobe=int((xprobe-xmin)/dx)+1
  jprobe=int((yprobe-ymin)/dy)+1
  uprobe=u(iprobe,jprobe,:)

  ! Compute mean zonal velocity inside mask
  xmaskMin=10.027753295 ; xmaskMax=10.192246705
  ymaskMin=-0.03079992784 ; ymaskMax=+0.03079992784
  umask=0.d0 ; cellCount=0
  do j=1,ny
     if ((y(j)>ymaskMin).and.(y(j)<ymaskMax)) then
        do i=1,nx
           if ((x(i)>xmaskMin).and.(x(i)<xmaskMax)) then
              cellCount=cellCount+1
              umask=umask+u(i,j,:)
           endif
        end do
     endif
  end do
  umask=umask/cellCount
  write (*,*) "   Number of cells included in mask: ",cellCount

  ! Compute mean zonal velocity 
  umean(1)=sum(u(1:nx,1:ny,1))/nx/ny
  if (nlayers==2) umean(2)=sum(u(1:nx,1:ny,2))/nx/ny

  inquire(file=filename, exist=fexist)
  if (fexist) then
     open(unit=2, file=filename, status="old", position="append")
  else
     open(unit=2, file=filename, status="unknown")
  endif
  if (nlayers==1) then
     hist_format = '(1X, 1PE12.5, 1X, 1PE12.5, 1X, 3(E14.5,1X))'
     write(2, hist_format) time, dt, uprobe(1), umask(1), umean(1)
  else
     hist_format = '(1X, 1PE12.5, 1X, 1PE12.5, 1X, 6(E14.5,1X))'
     write(2, hist_format) time, dt, uprobe(1), umask(1), umean(1), uprobe(2), umask(2), umean(2)
  endif
  close(2)

  return
end subroutine history
