!###########################################################
!###########################################################
!###########################################################
subroutine condinit
  use params
  use variables
  implicit none
  integer :: i,j
  real(dp) :: xvortex,yvortex,Avortex,rvortex,dshift !Vortex position & vorticity amplitude
  real(dp) :: vup,vlow,yjet,djet
  real(dp), dimension(0:nx+1,0:ny+1) :: qvortex,psivortex
  namelist /init_params/vup,vlow,yjet,djet,xvortex,yvortex,dshift,rvortex,Avortex

  if (verbose) write (*,*) "Entering condinit.f90..."

  ! Default params for the jet
  yjet=0.d0 ; djet=1.4d6
  vup=50.d0 ; vlow=25.d0

  ! Deformation radius, default
  lambda=4.5e5
  lambdaInvSq=1.d0/lambda**2

  ! Defaults vortex parameters
  xvortex=5000.d3 ; yvortex=-1000.d3 
  dshift=1.d6
  rvortex=7.d5
  Avortex=8.d-5

  ! Read problem input data
  open(unit=1,file='input',status='old')
  read(1,init_params)
  close(1)

  ! Background state
  do j=0,ny+1
     do i=0,nx+1
        psibar(i,j,1)=-vup *sqrt(3.1415)/2.d0*djet*erf((y(j)-yjet)/djet)
        psibar(i,j,2)=-vlow*sqrt(3.1415)/2.d0*djet*erf((y(j)-yjet)/djet)
     end do
  end do
  call getQ(psibar,qbar) ; call computeBCzeroGrad(qbar,nx,ny,nlayers)

  ! Add perturbation to jet
  qvortex=0. ; psivortex=1.
  do j=1,ny
     do i=0,nx+1
        qvortex(i,j)=Avortex*exp(-((x(i)-xvortex-dshift)**2+(y(j)-yvortex)**2)/rvortex**2)
     end do
  end do
  call computeBCzeroGrad(qvortex,nx,ny)
  !call poissonfft(qvortex,psivortex,.false.)
  psi(:,:,2)=psibar(:,:,2)+psivortex
  qvortex=0. 
  do j=1,ny
     do i=0,nx+1
        qvortex(i,j)=Avortex*exp(-((x(i)-xvortex)**2+(y(j)-yvortex)**2)/rvortex**2)
     end do
  end do
  call computeBCzeroGrad(qvortex,nx,ny)
  !call poissonfft(qvortex,psivortex,.false.)
  !psi(:,:,1)=psibar(:,:,1)+psivortex 
  !call computeBC(psi,nx,ny,nlayers)

  ! Compute PV
  !call getQ(psi,q) ; call computeBC(q,nx,ny,nlayers)
  q(:,:,1)=qbar(:,:,1)+qvortex
  q(:,:,2)=qbar(:,:,2)+qvortex

  ! Compute FFT and store variable in spectral space
  !call var2svar(psibar,spsibar,nx,ny,nlayers)
  call var2svar(q,sq,nx,ny,nlayers)
  call sq2spsi(sq,spsi,nx,ny,nlayers)
  call svar2var(spsi,psi,nx,ny,nlayers)
  call var2svar(psibar,spsibar,nx,ny,nlayers)
  !call svar2var(spsi,psi,nx,ny,nlayers)
  !call spsi2sq(spsi,sq,nx,ny,nlayers,LambdaInvSq)
  !call svar2var(sq,q,nx,ny,nlayers)
  call computeBCzeroGrad(q,nx,ny,nlayers)
  call computeBCzeroGrad(psi,nx,ny,nlayers)

  ! Compute other variables (qm1,u,v) for saving purposes.
  sqm1=sq
  call spsi2su(spsi,su,nx,ny,nlayers,2.d0*pi/(ymax-ymin)) ! Compute u in spectral space
  call su2u(su,u,nx,ny,nlayers)                           ! Transform u back to real space
  call spsi2sv(spsi,sv,nx,ny,nlayers,2.d0*pi/(xmax-xmin)) ! Compute v in spectral space
  call svar2var(sv,v,nx,ny,nlayers)                       ! Transform v back to real space

  return
end subroutine condinit
