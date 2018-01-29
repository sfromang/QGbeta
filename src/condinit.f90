!###########################################################
!###########################################################
!###########################################################
subroutine condinit
  use params
  use variables
  implicit none
  integer :: i,j
  namelist /init_params/

  ! Read problem input data
  open(unit=1,file='input' ,status='old')
  read(1,init_params)
  close(1)  

  ! Background state (here is dummy flow: u=v=0)
  do j=0,ny+1
     do i=0,nx+1
        psi(i,j,1:nlayers)=1.d0
     end do
  end do
  call getQ(psi,q) ; call computeBC(q,nx,ny,nlayers)

  ! Compute FFT and store variable in spectral space
  call var2spec(q,sq,nx,ny,nlayers)
  call var2spec(psi,spsi,nx,ny,nlayers)

  ! Compute other variables (psi,u,v) for saving purposes at initial time...
  qm1=q
  call poisson(q,psi,.false.) ; call computeBC(psi,nx,ny,nlayers)
  call getVel(psi,u,v) ; call computeBC(u,nx,ny,nlayers) ; call computeBC(v,nx,ny,nlayers)

  return
end subroutine condinit
