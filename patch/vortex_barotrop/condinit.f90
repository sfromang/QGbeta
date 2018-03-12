!###########################################################
!###########################################################
!###########################################################
subroutine condinit
  use params
  use variables
  implicit none
  integer :: i,j
  real(dp) :: x0,y0,dq !Vortex position & vorticity amplitude

  ! Background state
  do j=0,ny+1
     do i=0,nx+1
        !psi(i,j)=-37.5*sqrt(3.1415)/2./8.d-7*erf(8.d-7*(y(j)-6.4d5*sin(8.d-7*x(i)-3.1415/4.)))
        psi(i,j,1:nlayers)=-37.5*sqrt(3.1415)/2./8.d-7*erf(8.d-7*y(j))
     end do
  end do
  call getQ(psi,q) ; call computeBC(q,nx,ny,nlayers)

  ! Add perturbation to jet
  x0=2000.d3 ; y0=-800.d3
  do j=1,ny
     do i=1,nx
        dq=0.+1.5d-4*exp(-((x(i)-x0)**2+(y(j)-y0)**2)/2./460.e3**2)
        q(i,j,1:nlayers)=q(i,j,1:nlayers)+dq
     end do
  end do
  call computeBC(q,nx,ny,nlayers)

  ! Compute FFT and store variable in spectral space
  !call var2svar(psibar,spsibar,nx,ny,nlayers)
  call var2svar(q,sq,nx,ny,nlayers)
  call var2svar(psi,spsi,nx,ny,nlayers)
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
