!###########################################################
!###########################################################
!###########################################################
subroutine source_term
  use user_params
  use noise_params
  use params
  use variables
  implicit none
  real(dp), dimension(0:nx+1,0:ny+1) :: xi
  real(dp), dimension(0:nx+1,0:ny+1,nlayers) :: dtmp,dterm
  real(dp), dimension(0:nx+1,0:ny+1,1:nlayers)  :: dqdt_topo
  complex(dp), dimension(1:nx,1:ny,1:nlayers)  :: dsqdt_topo
  real(dp) :: psiR
  integer :: ilayer,i,j
  !forcing parameters
  real(dp) :: f1,f2
  real(dp), dimension(0:nx+1,0:ny+1) :: qstar,dtermqstar
  !random number parameters
  real(dp) :: amp,sign

  if (verbose) write (*,*) 'Entering source_term...'

  ! Jet forcing term
  call getForcingAmp(time)
  if (varForcing) call computeForcing
  do ilayer=1,nlayers
     dsqdt(:,:,ilayer)=dsqdt(:,:,ilayer)-kappa*(sq(:,:,ilayer)-sqbar(:,:,ilayer))
  end do

  ! Bottom layer topography
  dqdt_topo=0.d0 ; dsqdt_topo=0.d0
  call su2u(su,u,nx,ny,nlayers)
  call svar2var(sv,v,nx,ny,nlayers)
  do j=1,ny
     do i=1,nx
        dqdt_topo(i,j,nlayers) = - (hB(i+1,j)*u(i+1,j,nlayers)-hB(i-1,j)*u(i-1,j,nlayers))/2.d0/dx \
                                 - (hB(i,j+1)*v(i,j+1,nlayers)-hB(i,j-1)*v(i,j-1,nlayers))/2.d0/dy
     end do
  end do
  call var2svar(dqdt_topo,dsqdt_topo,nx,ny,nlayers)
  dsqdt=dsqdt+dsqdt_topo

  ! TBD FOR FFT VERSION - SEB - 12/02/18
  ! ! Noise forcing (amp=1.d0 for first model that worked)
  ! !call getNoise(.true.,.true.,dt)
  ! call getNoiseAmp(time)            ! return noise amplitude
  ! call getNoise(dt) ! return noise with values between -0.5 and 0.5
  ! do j=1,ny
  !    do i=1,nx
  !       !call random_number(rvalue)
  !       do ilayer=1,nlayers
  !          sign=(-1)**(ilayer+1)
  !          dqdt(i,j,ilayer)=dqdt(i,j,ilayer)+2.d0*sign*noiseAmp*noise(i,j)
  !       end do
  !    end do
  ! end do

  return
end subroutine source_term
!===============================================================================
!> Return noise amplitude as a function of time
!===============================================================================
subroutine getNoiseAmp(time)
  use params , only : tlim
  use noise_params
  use precision
  implicit none
  real(dp), intent(in) :: time
  integer :: Nstep,tstep 

  if (varNoise) then
     Nstep=int((noiseMaxAmp-noiseMinAmp)/noiseInc)+1
     tstep=int(tlim/Nstep)
     noiseAmp=noiseMinAmp+noiseInc*int(time/tstep)
  endif

  write (*,*) "           Noise amplitude: ",noiseAmp
  
  return
end subroutine getNoiseAmp
!===============================================================================
!> Generate 2D noise map - with option to have white or red noise
!===============================================================================
subroutine getNoise(dt)
  use params
  use fftVar
  use noise_params
  implicit none
  ! Input parameters
  real(dp), intent(in) :: dt
  ! Local variables
  real(dp) :: waveno,rcoeff,tcorr
  integer :: i,j,di,dj
  complex(dp), dimension(1:nx,1:ny) :: fftNoise
  real(dp) :: rvalue

  if (.not.(allocated(noise   ))) allocate(noise   (1:nx,1:ny))
  if (.not.(allocated(savNoise))) allocate(savNoise(1:nx,1:ny))

  ! First generate white noise
  do j=1,ny
     do i=1,nx
        call random_number(rvalue)
        noise(i,j)=rvalue
     end do
  end do
  if ((.not.(redSpace)).and.(.not.(redTime))) noise=noise-5.d-1

  ! ! Next generate red noise (in space)
  ! ! the cut-off frequency w0 corresponds to a cut-off wavelength lambda ~ 1/w0 expressed in cells!!!
  ! ! so beware that the noise is resolution dependent!!!
  ! ! when doubling the resolution, w0 must be decreased by a factor of two to keep noise the same
  ! if (redSpace) then
  !    fftNoise=noise(1:nx,1:ny)
  !    call cfft2f(nx,nx,ny,fftNoise,wsave,lensav,work,lenwrk,ierr)
  !    fftNoise(1,1)=0.d0
  !    do j=1,ny
  !       do i=1,nx
  !          di = nx/2 + 1 - i
  !          dj = ny/2 + 1 - j
  !          waveno = sqrt( (0.5d0-abs(real(di))/real(nx))**2 + (0.5d0-abs(real(dj))/real(ny))**2 )
  !          fftNoise(i,j) = fftNoise(i,j) * exp(-((waveno-w0)/dw0)**2) ! Gaussian cut-off
  !          !fftNoise(i,j) = fftNoise(i,j) * exp(-((waveno-0.04d0)/0.01d0)**2) ! Gaussian cut-off
  !          !fftNoise(i,j) = fftNoise(i,j) * exp(-(abs(waveno/0.05d0))) ! Red noise
  !       end do
  !    end do
  !    call cfft2b(nx,nx,ny,fftNoise,wsave,lensav,work,lenwrk,ierr)
  !    noise(:,:) = real(fftNoise)
  !    noise = noise / max(maxval(noise),-minval(noise)) / 2.d0
  ! endif

  ! ! Finally generate red noise in time
  ! if (redTime) then
  !    ! Correlation coeff for temporally red noise
  !    !tcorr= 0.1d0 !1.d0 (ok transition), with wavenumber 0.05
  !    tcorr= 1.d0 !(ok transition), with wavenumber 0.05
  !    rcoeff=0.d0 !0.95d0
  !    !rcoeff=1.d0-dt/tcorr
  !    noise=savNoise*rcoeff+sqrt(1.d0-rcoeff**2)*noise
  ! endif

  ! Save copy of noise
  savNoise=noise

  return
end subroutine getNoise
