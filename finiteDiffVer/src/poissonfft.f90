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
subroutine poissonfft(q,psi,withSource,periodic,neumann,dirichlet)
  use params
  use fftVar
  use user_params
  implicit none
  logical , intent(in) :: withSource
  logical , intent(in) :: periodic,dirichlet,neumann
  real(dp), dimension(0:nx+1,0:ny+1), intent(in ) :: q
  real(dp), dimension(0:nx+1,0:ny+1), intent(out) :: psi
  complex(dp), dimension(1:nx,1:ny) :: fftPsi
  real(dp), dimension(1:nx,1:ny) :: fftPsiR
  real(dp), dimension(1:ny) :: tmpYreal,tmpYimag
  real(dp) :: scalex,scaley,kx,ky,ktotSq,rsq,factor
  integer :: i,j

  rsq=dx/dy ; rsq=rsq*rsq

  ! Compute forward FFT transform for q
  if (periodic) then
     fftPsi=q(1:nx,1:ny)
     call cfft2f(nx,nx,ny,fftPsi,wsave,lensav,work,lenwrk,ierr)
  endif
  if (neumann) then
     fftPsiR=q(1:nx,1:ny)
     do i=1,nx
        call cost1f(ny,1,fftPsiR(i,:),ny,wsaveyN,lensavyN,workyN,lenwrkyN,ierr)
     end do
     fftPsi=fftPsiR
     do j=1,ny
        call cfft1f(nx,1,fftPsi(:,j),nx,wsavex,lensavx,workx,lenwrkx,ierr)
     end do
  endif
  if (dirichlet) then
     fftPsi=q(1:nx,1:ny)
     do j=1,ny
        call cfft1f(nx,1,fftPsi(:,j),nx,wsavex,lensavx,workx,lenwrkx,ierr)
     end do
     do i=1,nx
        tmpYreal=real(fftPsi(i,:))
        tmpYimag=aimag(fftPsi(i,:))
        if (i.gt.1) then
           call sint1f(ny,1,tmpYreal,ny,wsaveyD,lensavyD,workyD,lenwrkyD,ierr)
           call sint1f(ny,1,tmpYimag,ny,wsaveyD,lensavyD,workyD,lenwrkyD,ierr)
        else
           call cost1f(ny,1,tmpYreal,ny,wsaveyN,lensavyN,workyN,lenwrkyN,ierr)
           call cost1f(ny,1,tmpYimag,ny,wsaveyN,lensavyN,workyN,lenwrkyN,ierr)
        endif
        do j=1,ny
           fftPsi(i,j)=complex(tmpYreal(j),tmpYimag(j))
        end do
     end do
  endif

  ! Solve Poisson equation in Fourier space (include lambdaInvSq...).
  factor=1.d0/dx**2+1.d0/dy**2
  if (withSource) factor=factor+lambdaInvSq
  do j=1,ny
     do i=1,nx
        if (periodic) then
           if ((i==1).and.(j==1)) then
              fftPsi(i,j)=0.d0
           else
              fftPsi(i,j)=half*fftPsi(i,j)/(cos(twopi*(i-1)/nx)/dx**2+cos(twopi*(j-1)/ny)/dy**2-factor)
           endif
        endif
        if (neumann) then
           if (((i==1).and.(j==1))) then
              fftPsi(i,j)=0.d0
           else
              fftPsi(i,j)=half*fftPsi(i,j)/(cos(twopi*(i-1)/nx)/dx**2+cos(pi*(j-1)/(ny-1))/dy**2-factor)
           endif
        endif
        if (dirichlet) then
           if ((i==1).and.(j==1)) then
              fftPsi(i,j)=0.d0
           else if ((i==1).and.(j.gt.1)) then
               fftPsi(i,j)=half*fftPsi(i,j)/(1.d0/dx**2+cos(pi*(j-1)/(ny-1))/dy**2-factor)
           else if (i.gt.1) then
              fftPsi(i,j)=half*fftPsi(i,j)/(cos(twopi*(i-1)/nx)/dx**2+cos(pi*j/(ny+1))/dy**2-factor)
           endif
        endif
     end do
  end do

  ! Compute backward FFT transform to get psi
  if (periodic) then
     call cfft2b(nx,nx,ny,fftPsi,wsave,lensav,work,lenwrk,ierr)
     psi(1:nx,1:ny)=real(fftPsi)
  endif
  if (neumann) then
     do j=1,ny
        call cfft1b(nx,1,fftPsi(:,j),nx,wsavex,lensavx,workx,lenwrkx,ierr)
     end do
     fftPsiR=real(fftPsi)
     do i=1,nx
        call cost1b(ny,1,fftPsiR(i,:),ny,wsaveyN,lensavyN,workyN,lenwrkyN,ierr)
     end do
     psi(1:nx,1:ny)=fftPsiR
  endif
  if (dirichlet) then
     do i=1,nx
        tmpYreal=real(fftPsi(i,:))
        tmpYimag=aimag(fftPsi(i,:))
        if (i.gt.1) then
           call sint1b(ny,1,tmpYreal,ny,wsaveyD,lensavyD,workyD,lenwrkyD,ierr)
           call sint1b(ny,1,tmpYimag,ny,wsaveyD,lensavyD,workyD,lenwrkyD,ierr)
        else
           call cost1b(ny,1,tmpYreal,ny,wsaveyN,lensavyN,workyN,lenwrkyN,ierr)
           call cost1b(ny,1,tmpYimag,ny,wsaveyN,lensavyN,workyN,lenwrkyN,ierr)
        endif
        do j=1,ny
           fftPsi(i,j)=complex(tmpYreal(j),tmpYimag(j))
        end do
     end do
     do j=1,ny
        call cfft1b(nx,1,fftPsi(:,j),nx,wsavex,lensavx,workx,lenwrkx,ierr)
     end do
     fftPsiR=real(fftPsi)
     psi(1:nx,1:ny)=fftPsiR
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
  call init1dcosty
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
! subroutine init1dffty
!   use params
!   use fftVar
!   implicit none

!   if (verbose) write (*,*) 'Entering init1dffty...'

!   lensavy=2*ny+4+int(log(real(ny))/log(2.))
!   allocate(wsavey(lensavy)) ; wsavey=0.D0
!   call cfft1i(ny,wsavey,lensavy,ierr)
!   lenwrky=2*ny
!   allocate(worky(lenwrky)) ; worky=0.D0

!   return
! end subroutine init1dffty
!===============================================================================
!> Initialize arrays for 1D sin transform (y-direction)
!===============================================================================
subroutine init1dsinty
  use params
  use fftVar
  implicit none

  if (verbose) write (*,*) 'Entering init1dsinty...'

  lensavyD=ny/2+ny+4+int(log(real(ny))/log(2.))
  allocate(wsaveyD(lensavyD)) ; wsaveyD=0.D0
  call sint1i(ny,wsaveyD,lensavyD,ierr)
  lenwrkyD=2*ny+2
  allocate(workyD(lenwrkyD)) ; workyD=0.D0

  return
end subroutine init1dsinty
!===============================================================================
!> Initialize arrays for 1D cos transform (y-direction)
!===============================================================================
subroutine init1dcosty
  use params
  use fftVar
  implicit none

  if (verbose) write (*,*) 'Entering init1dcosty...'

  lensavyN=2*ny+ny+4+int(log(real(ny))/log(2.))
  allocate(wsaveyN(lensavyN)) ; wsaveyN=0.D0
  call cost1i(ny,wsaveyN,lensavyN,ierr)
  lenwrkyN=2*ny+2
  allocate(workyN(lenwrkyN)) ; workyN=0.D0

  return
end subroutine init1dcosty
