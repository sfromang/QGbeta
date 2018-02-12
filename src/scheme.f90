!###########################################################
!###########################################################
!###########################################################
subroutine getUpdate
  use params
  use variables
  use user_params
  implicit none
  real(dp), dimension(1:nx,1:ny) :: factor
  integer :: ilayer

  ! Solve Fourier equation in Fourier space (i.e. compute spsi from sq)
  call sq2spsi(sq,spsi,nx,ny,nlayers)

  ! Compute velocities from streamfunction & timestep
  call spsi2su(spsi,su,nx,ny,nlayers,2.d0*pi/(ymax-ymin)) ! Compute u in spectral space
  call spsi2sv(spsi,sv,nx,ny,nlayers,2.d0*pi/(xmax-xmin)) ! Compute v in spectral space

  ! Compute dt according to CFL condition (if needed)
  if (cfl>0.d0) then
     call su2u(su,u,nx,ny,nlayers)                           ! Transform u back to real space
     call svar2var(sv,v,nx,ny,nlayers)                       ! Transform v back to real space
     dt=cfl*min(dx/maxval(abs(u)),dy/maxval(abs(v)))
  endif

  ! Compute time derivative from velocities and vorticity
  call computeTimeDerivative(sq,spsi,su,sv,dsqdt)

  ! Add source term
  call source_term

  ! ! Compute q(t+dt) from q(t-dt) & dqdt using leapfrog scheme
  ! if (time>0.d0) then
  !    sqp1=sqm1+2.d0*dsqdt*dt
  ! else
  !    sqp1=sqm1+dsqdt*dt
  ! endif
  ! dsqdt=0.d0

  ! ! RAW filter
  ! call rawFilter

  !sqm1=sq ; sq=sqp1 ; spsim1=spsi

  ! Compute q(t+dt) using 3rd order Adam-Bashforth scheme
  if (time>0.d0) then
     sqp1=sq+dt/12.d0*(23.d0*dsqdt-16.d0*dsqdtm1+5.d0*dsqdtm2)
  else
     sqp1=sq+dsqdt*dt
  endif
  dsqdtm2=dsqdtm1 ; dsqdtm1=dsqdt ; dsqdt=0.d0
  sqm1=sq ; sq=sqp1 ; spsim1=spsi
  

end subroutine getUpdate
!###########################################################
!###########################################################
!###########################################################
subroutine computeTimeDerivative(sq,spsi,su,sv,dsqdt)
  use params
  implicit none

  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(in ) :: sq,spsi,su,sv
  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(out) :: dsqdt

  complex(dp), dimension(1:nx,1:ny,1:nlayers) :: sduqdx,sdvqdy,dsqdy,dsqdx
  real(dp), dimension(0:nx+1,0:ny+1,1:nlayers) :: u,v,q,dqdx,dqdy

  ! Dealisasing
  call dealiasing(sq,nx,ny,nlayers)
  call dealiasing(su,nx,ny,nlayers)
  call dealiasing(sv,nx,ny,nlayers)

  ! PV advection by velocities
  if (typeBC==1) then

     ! Compute d(uq)/dx in spectral space
     call spec2var(su,u,nx,ny,nlayers)
     call spec2var(sq,q,nx,ny,nlayers)
     call var2spec(u*q,sduqdx,nx,ny,nlayers)
     call derive(sduqdx,nx,ny,nlayers,1)

     ! Compute d(vq)/dy in spectral space
     call spec2var(sv,v,nx,ny,nlayers)
     call var2spec(v*q,sdvqdy,nx,ny,nlayers)
     call derive(sdvqdy,nx,ny,nlayers,2)

  else if (typeBC==2) then

     ! Compute u*dq/dx in spectral space
     call su2u(su,u,nx,ny,nlayers)
     dsqdx=sq
     call deriveTypeII(dsqdx,nx,ny,nlayers,2.d0*pi/(xmax-xmin),1)
     call svar2var(dsqdx,dqdx,nx,ny,nlayers)
     call var2svar(u*dqdx,sduqdx,nx,ny,nlayers)

     ! Compute v*dq/dy in spectral space
     call svar2var(sv,v,nx,ny,nlayers)
     dsqdy=sq
     call deriveTypeII(dsqdy,nx,ny,nlayers,2.d0*pi/(ymax-ymin),2)
     call su2u(dsqdy,dqdy,nx,ny,nlayers)
     call var2svar(v*dqdy,sdvqdy,nx,ny,nlayers)

  endif

  dsqdt=dsqdt-sduqdx-sdvqdy-beta*sv

  ! Compute small scale dissipation
  call dissip(spsi,dsqdt)

  return
end subroutine computeTimeDerivative
!###########################################################
!###########################################################
!###########################################################
subroutine dealiasing(svar,nx,ny,nlayers)
  use precision
  use params , only : typeBC
  implicit none
  integer, intent(in) :: nx,ny,nlayers
  complex(dp), dimension(nx,ny,nlayers) :: svar

  if (typeBC==2) then
     svar(nx/3:nx/2+nx/6,      :  ,:)=0.d0
     svar(    :         ,2*ny/3:ny,:)=0.d0
  endif

  return
end subroutine dealiasing
!###########################################################
!###########################################################
!###########################################################
subroutine dissip(spsi,dsqdt)
  use params
  implicit none

  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(in ) :: spsi
  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(out) :: dsqdt
  complex(dp), dimension(nx,ny,nlayers) :: sdissip

  ! 6th order small scale hyperdiffusion
  sdissip=spsi
  call delSq(sdissip,nx,ny,nlayers)
  call delSq(sdissip,nx,ny,nlayers)
  call delSq(sdissip,nx,ny,nlayers)
  sdissip=nuH*sdissip

  dsqdt = dsqdt - sdissip

  return
end subroutine dissip
!###########################################################
!###########################################################
!###########################################################
subroutine rawFilter
  use params
  use variables
  implicit none
  integer :: i,j,ilayer
  real(dp) :: d,alpha

  nu=0.2 ; alpha=0.53
  do ilayer=1,nlayers
     do j=1,ny
        do i=1,nx
           ! Compute filter displacement
           d = nu*(sqm1(i,j,ilayer)-2.d0*sq(i,j,ilayer)+sqp1(i,j,ilayer))/2.d0
           ! Apply filter
           sq  (i,j,ilayer) = sq  (i,j,ilayer) + d*alpha
           sqp1(i,j,ilayer) = sqp1(i,j,ilayer) + d*(alpha-1.d0)
        end do
     end do
  end do

  return
end subroutine rawFilter
!###########################################################
!###########################################################
!###########################################################
subroutine sq2spsi(sq,spsi,nx,ny,nlayers)
  use const
  implicit none
  integer :: nx,ny,nlayers
  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(in) :: sq
  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(out) :: spsi

  complex(dp), dimension(1:nx,1:ny) :: sq_s,sq_d,spsi_s,spsi_d
  integer :: ilayer

  ! Initialize arrays
  spsi_s=0.d0 ; spsi_d=0.d0

  ! Compute solution of Poisson equation
  sq_s=half*(sq(:,:,1)+sq(:,:,nlayers)) ; sq_d=half*(sq(:,:,1)-sq(:,:,nlayers))
  call poissonfft(sq_s,spsi_s,.false.)
  if (nlayers==2) call poissonfft(sq_d,spsi_d,.true.) 
  spsi(:,:,1)=spsi_s+spsi_d
  if (nlayers==2) spsi(:,:,2)=spsi_s-spsi_d

  return
end subroutine sq2spsi
!===============================================================================
!> Compute PV in spectral space from PSI in spectral space
!===============================================================================
subroutine spsi2sq(spsi,sq,nx,ny,nlayers,LambdaInvSq)
  use const
  implicit none
  integer :: nx,ny,nlayers
  complex(dp), dimension(1:nx,1:ny,1:nlayers) :: spsi,sq
  real(dp) :: LambdaInvSq

  integer :: ilayer

  sq = spsi
  call delSq(sq,nx,ny,nlayers)
  if (nlayers == 2) then
     sq(:,:,1)=sq(:,:,1)-(spsi(:,:,1)-spsi(:,:,2))*LambdaInvSq        
     sq(:,:,2)=sq(:,:,2)+(spsi(:,:,1)-spsi(:,:,2))*LambdaInvSq
  endif

  return
end subroutine spsi2sq
!###########################################################
!###########################################################
!###########################################################
! subroutine getVel(spsi,su,sv)
!   use params
!   implicit none
!   complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(in ) :: spsi
!   complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(out) :: su,sv
!   integer :: i,j,im1,ip1,jm1,jp1

!   do j=1,ny
!      su(:,j,:)=complex(0.d0,-1.d0)*ky(j)*spsi(:,j,:) 
!   end do
!   do i=1,nx
!      sv(i,:,:)=complex(0.d0,+1.d0)*kx(i)*spsi(i,:,:) 
!   end do

!   return
! end subroutine getVel
subroutine getVel(psi,u,v)
  use params
  implicit none
  real(dp), dimension(0:nx+1,0:ny+1,1:nlayers), intent(in ) :: psi
  real(dp), dimension(0:nx+1,0:ny+1,1:nlayers), intent(out) :: u,v
  integer :: i,j,im1,ip1,jm1,jp1

  do j=1,ny
     jm1=j-1 ; jp1=j+1
     do i=1,nx
        im1=i-1 ; ip1=i+1
        u(i,j,:)=-(psi(i  ,jp1,:)-psi(i  ,jm1,:))/2.d0/dy
        v(i,j,:)=+(psi(ip1,j  ,:)-psi(im1,j  ,:))/2.d0/dx
     end do
  end do
  call computeBCzeroGrad(u,nx,ny,nlayers)
  call computeBC(v,nx,ny,nlayers)
  v(:,0,:)=0.d0 ; v(:,ny+1,:)=0.d0

  return
end subroutine getVel
!###########################################################
!###########################################################
!###########################################################
subroutine getQ(psi,q)
  use params
  implicit none
  real(dp), dimension(0:nx+1,0:ny+1,1:nlayers), intent(in ) :: psi
  real(dp), dimension(0:nx+1,0:ny+1,1:nlayers), intent(out) :: q
  integer :: i,j,im1,ip1,jm1,jp1

  q=0.d0
  do j=1,ny
     jm1=j-1 ; jp1=j+1
     do i=1,nx
        im1=i-1 ; ip1=i+1
        q(i,j,:)=(psi(ip1,j,:)-2.d0*psi(i,j,:)+psi(im1,j,:))/dx**2+\
                 (psi(i,jp1,:)-2.d0*psi(i,j,:)+psi(i,jm1,:))/dy**2
     end do
  end do
  if (nlayers == 2) then
     q(:,:,1)=q(:,:,1)-(psi(:,:,1)-psi(:,:,2))*LambdaInvSq        
     q(:,:,2)=q(:,:,2)+(psi(:,:,1)-psi(:,:,2))*LambdaInvSq
  endif

  return
end subroutine getQ
!###########################################################
!###########################################################
!###########################################################
subroutine getVorticity(psi,xi)
  use params
  implicit none
  real(dp), dimension(0:nx+1,0:ny+1), intent(in ) :: psi
  real(dp), dimension(0:nx+1,0:ny+1), intent(out) :: xi
  integer :: i,j,im1,ip1,jm1,jp1

  xi=0.d0
  do j=1,ny
     jm1=j-1 ; jp1=j+1
     do i=1,nx
        im1=i-1 ; ip1=i+1
        xi(i,j)=(psi(ip1,j)-2.d0*psi(i,j)+psi(im1,j))/dx**2+\
                (psi(i,jp1)-2.d0*psi(i,j)+psi(i,jm1))/dy**2
     end do
  end do

  return
end subroutine getVorticity
!###########################################################
!###########################################################
!###########################################################
subroutine getDelSq(phi,delSqPhi)
  use params
  implicit none
  real(dp), dimension(0:nx+1,0:ny+1,nlayers), intent(in ) :: phi
  real(dp), dimension(0:nx+1,0:ny+1,nlayers), intent(out) :: delSqPhi
  integer :: i,j,im1,ip1,jm1,jp1

  delSqPhi=0.d0
  do j=1,ny
     jm1=j-1 ; jp1=j+1
     do i=1,nx
        im1=i-1 ; ip1=i+1
        delSqPhi(i,j,:)=(phi(ip1,j  ,:)-2.d0*phi(i,j,:)+phi(im1,j,:))/dx**2+\
                        (phi(i  ,jp1,:)-2.d0*phi(i,j,:)+phi(i,jm1,:))/dy**2
     end do
  end do

  return
end subroutine getDelSq
