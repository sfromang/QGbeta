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

  if (verbose) print*, "> Writing history..."

  inquire(file=filename, exist=fexist)
  if (fexist) then
     open(unit=2, file=filename, status="old", position="append")
  else
     open(unit=2, file=filename, status="unknown")
  endif
  hist_format = '(1X, 1PE12.5, 1X, 1PE12.5)'
  write(2, hist_format) time, dt
  close(2)

  return
end subroutine history
