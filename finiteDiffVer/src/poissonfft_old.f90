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
!! \b last \b modified: 18-06-2015
!===============================================================================
!> Solve Poisson equation
!> if "withsource" is False, solve:
!> del^2 psi = q(i,j)
!> if "withsource" is True, solve:
!> del^2 psi = q(i,j) + 2/lambda^2*psi(i,j) 
!===============================================================================
subroutine poissonfft(q,psi,withSource,periodic)
  use params
  use fftVar
  implicit none
  logical , intent(in)                            :: withSource,periodic
  real(dp), dimension(0:nx+1,0:ny+1), intent(in ) :: q
  real(dp), dimension(0:nx+1,0:ny+1), intent(out) :: psi
  complex, dimension(1:nx,1:ny) :: fftPsi
  real, dimension(1:nx,1:ny) :: fftPsiR
  real(dp) :: scalex,scaley,kx,ky,ktotSq
  integer :: i,j

  ! Compute forward FFT transform for q
  if (periodic) then
     fftPsi=q(1:nx,1:ny)
     call cfft2f(nx,nx,ny,fftPsi,wsave,lensav,work,lenwrk,ierr)
  else
     fftPsiR=q(1:nx,1:ny)
     do i=1,nx
        !call cfft1f(ny,1,fftPsi(i,:),ny,wsavey,lensavy,worky,lenwrky,ierr)
        call sint1f(ny,1,fftPsiR(i,:),ny,wsavey,lensavy,worky,lenwrky,ierr)
     end do
     fftPsi=fftPsiR
     do j=1,ny
        call cfft1f(nx,1,fftPsi(:,j),nx,wsavex,lensavx,workx,lenwrkx,ierr)
     end do
  endif

  ! Compute scaling (beware kx=ky=0 special)
  scalex=2.d0*pi/(xmax-xmin)
  scaley=2.d0*pi/(ymax-ymin)

  ! Solve Poisson equation in Fourier space (include lambdaInvSq...).
  do j=1,ny
     do i=1,nx
        if (i.le.nx/2+1) then
           kx=(i-1)*scalex
        else
           kx=-(nx-i+1)*scalex
        endif
        if (j.le.ny/2+1) then
           ky=(j-1)*scaley
        else
           ky=-(ny-j+1)*scaley
        endif
        ktotSq=-kx*kx-ky*ky
        if (withSource) ktotSq=ktotSq-2.d0*lambdaInvSq
!        if ((i>1).or.(j>1)) fftPsi(i,j)=fftPsi(i,j)/ktotSq
!        if ((i>1).or.(j>1)) then
!           fftPsi(i,j)=half*fftPsi(i,j)*dx*dy/(cos(2.d0*pi*i/(nx+1))+cos(2.d0*pi*j/(ny+1))-2.d0)
!           fftPsi(i,j)=half*fftPsi(i,j)*dx*dy/(cos(2.d0*pi*i/(nx+1))+cos(pi*j/(ny+1))-2.d0)
!        endif
        if (periodic) then
           fftPsi(i,j)=half*fftPsi(i,j)*dx*dy/(cos(2.d0*pi*i/(nx+1))+cos(2.d0*pi*j/(ny+1))-2.d0)
        else
           fftPsi(i,j)=half*fftPsi(i,j)*dx*dy/(cos(2.d0*pi*i/(nx+1))+cos(pi*j/(ny+1))-2.d0)
        endif
     end do
  end do
!  fftPsi(1,1)=complex(0.d0,0.d0)

  ! Compute backward FFT transform to get psi
  if (periodic) then
     call cfft2b(nx,nx,ny,fftPsi,wsave,lensav,work,lenwrk,ierr)
     psi(1:nx,1:ny)=real(fftPsi)
  else
     do j=1,ny
        call cfft1b(nx,1,fftPsi(:,j),nx,wsavex,lensavx,workx,lenwrkx,ierr)
     end do
     fftPsiR=real(fftPsi)
     do i=1,nx
        !call cfft1b(ny,1,fftPsi(i,:),ny,wsavey,lensavy,worky,lenwrky,ierr)
        call sint1b(ny,1,fftPsiR(i,:),ny,wsavey,lensavy,worky,lenwrky,ierr)
     end do
     psi(1:nx,1:ny)=real(fftPsiR)
  endif

  return
end subroutine poissonfft
!===============================================================================
!> Initialize arrays for FFT depending on boundary conditions
!===============================================================================
subroutine initfft
  implicit none
  call init2dfft
  call init1dfftx
  call init1dsinty
  !call init1dffty
  return
end subroutine initfft
!===============================================================================
!> Initialize arrays for 2D FFT complex transform
!===============================================================================
subroutine init2dfft
  use params
  use fftVar
  implicit none

  if (verbose) write (*,*) 'Entering init2dfft...'

  lensav=2*(nx+ny)+16+int(log(real(nx)/log(2.)))+int(log(real(ny)/log(2.)))
  allocate(wsave(lensav)) ; wsave=0.D0
  call cfft2i(nx,ny,wsave,lensav,ierr)
  lenwrk=2*nx*ny
  allocate(work(lenwrk)) ; work=0.D0

  return
end subroutine init2dfft
!===============================================================================
!> Initialize arrays for 1D FFT complex transform (x-direction)
!===============================================================================
subroutine init1dfftx
  use params
  use fftVar
  implicit none

  if (verbose) write (*,*) 'Entering init1dfftx...'

  lensavx=2*nx+4+int(log(real(nx))/log(2.))
  allocate(wsavex(lensavx)) ; wsavex=0.
  call cfft1i(nx,wsavex,lensavx,ierr)
  lenwrkx=2*nx
  allocate(workx(lenwrkx)) ; workx=0.D0

  return
end subroutine init1dfftx
!===============================================================================
!> Initialize arrays for 1D FFT complex transform (y-direction)
!===============================================================================
subroutine init1dffty
  use params
  use fftVar
  implicit none

  if (verbose) write (*,*) 'Entering init1dffty...'

  lensavy=2*ny+4+int(log(real(ny))/log(2.))
  allocate(wsavey(lensavy)) ; wsavey=0.D0
  call cfft1i(ny,wsavey,lensavy,ierr)
  lenwrky=2*ny
  allocate(worky(lenwrky)) ; worky=0.D0

  return
end subroutine init1dffty
!===============================================================================
!> Initialize arrays for 1D sin transform (y-direction)
!===============================================================================
subroutine init1dsinty
  use params
  use fftVar
  implicit none

  if (verbose) write (*,*) 'Entering init1dsinty...'

  lensavy=ny/2+ny+4+int(log(real(ny))/log(2.))
  allocate(wsavey(lensavy)) ; wsavey=0.D0
  call sint1i(ny,wsavey,lensavy,ierr)
  lenwrky=2*ny+2
  allocate(worky(lenwrky)) ; worky=0.D0

  return
end subroutine init1dsinty
