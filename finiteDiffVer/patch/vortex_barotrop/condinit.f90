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

  ! Compute other variables (psi,u,v) for saving purposes.
  qm1=q
  call poisson(q,psi,.false.) ; call computeBC(psi,nx,ny,nlayers)
  call getVel(psi,u,v) ; call computeBC(u,nx,ny,nlayers) ; call computeBC(v,nx,ny,nlayers)

  return
end subroutine condinit
