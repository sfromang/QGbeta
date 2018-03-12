!###########################################################
!###########################################################
!###########################################################
subroutine output
  use noise_params
  use params
  use variables
  implicit none
  character(LEN=80) :: filename

  if (verbose) write (*,*) 'Entering output...'

  !if (allocated(noise)) u(1:nx,1:ny,1)=noise

  call get_filename(ndump,0,'qgout',filename)
  open(unit=10,file=filename,status='unknown',form='unformatted')
  write(10) time
  write(10) nx,ny,nlayers
  write(10) beta,lambda
  write(10) x
  write(10) y
  write(10) q
  write(10) psi
  write(10) u
  write(10) v
  write(10) qbar
  write(10) psibar
  write(10) qm1
  close(10)

  return
end subroutine output
!###########################################################
!###########################################################
!###########################################################
subroutine restart_run
  use params
  use variables
  implicit none
  character(LEN=80) :: filename

  if (verbose) write (*,*) 'Entering restart subroutine...'

  call get_filename(ndump,0,'qgout',filename)
  open(unit=10,file=filename,status='old',form='unformatted')
  read(10) time
  read(10) nx,ny,nlayers
  read(10) beta,lambda
  read(10) x
  read(10) y
  read(10) q
  read(10) psi
  read(10) u
  read(10) v
  read(10) qbar
  read(10) psibar
  read(10) qm1
  close(10)

  return
end subroutine restart_run
!###########################################################
!###########################################################
!###########################################################
subroutine get_filename(ndump,mype,prefix,filename)
  implicit none
!#if WITHMPI==1
!#include "mpif.h"
!#endif

  integer :: ndump,mype

  integer :: info
  character(LEN=6) :: snumfile,n_mype
  character(LEN=*) :: prefix
  character(LEN=80) :: filename,filedir,filecmd
!#if WITHMPI==1
!  integer :: ierr
!#endif

  call convtoasc(ndump,snumfile)
  call convtoasc(mype,n_mype  )
  filedir='output_'//trim(snumfile)//'/'
  filecmd='mkdir -p '//trim(filedir)
!#ifdef NOSYSTEM
!  if (mype==0) call PXFMKDIR(trim(filedir),len(trim(filedir)),O'755',info)
!#elif BLUEGENE==1
!  if (mype==0) call mkdir(trim(filedir)//'\0', %val(511)) ! dec(511) = oct(777)
!#else
  if (mype==0) call system(filecmd)
!#endif
!#if WITHMPI==1
!  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!#endif
  filename=trim(filedir)//trim(prefix)//'.'//trim(n_mype)
  
  return
end subroutine get_filename
!###########################################################
!###########################################################
!###########################################################
subroutine convtoasc(number,sstring)
!=======================================================================
! To convert an integer smaller than 999999 to a 6 characters string
!=======================================================================
  implicit none
  integer :: number, istring, num, nums10, i
  character(LEN=6) :: sstring
  character(LEN=10),parameter :: nstring="0123456789"
     
  num=1000000
  nums10=num/10
  do i=1,6
      istring=1+mod(number,num)/nums10
      sstring(i:i)=nstring(istring:istring)
      num=num/10
      nums10=nums10/10
  enddo
  
end subroutine convtoasc

