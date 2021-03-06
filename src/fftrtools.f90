!=================================================================================
!> Compute svar in spectral space from var in real space (from psi-like variables)
!=================================================================================
subroutine var2svar(var,svar,nx,ny,nlayers)
  use, intrinsic :: iso_c_binding
  use precision
  implicit none
  include "fftw3.f03"

  integer,intent(in) :: nx,ny,nlayers
  real(dp), dimension(0:nx+1,0:ny+1,1:nlayers), intent(in) :: var
  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(out) :: svar

  integer :: ilayer,i,j
  complex(dp), dimension(1:nx) :: varIn,varOut
  real(dp), dimension(:),allocatable :: fIn,fOutReal,fOutImag
  type(C_PTR) :: plan
 
  ! Loop over layers
  do ilayer=1,nlayers

     ! Complex transform data in x using complex array
     do j=1,ny
        varIn=var(1:nx,j,ilayer) ! Convert input array to complex
        plan=fftw_plan_dft_1d(nx,varIn,varOut,FFTW_FORWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(plan,varIn,varOut)
        call fftw_destroy_plan(plan)
        svar(:,j,ilayer)=varOut
     end do

     ! Cosine transform for the i=1 component in y (beware real array needed in input)
     allocate(fIn(ny),fOutReal(ny),fOutImag(ny))
     fIn=real(svar(1,:,ilayer))
     plan=fftw_plan_r2r_1d(ny,fIn,fOutReal,FFTW_REDFT10,FFTW_ESTIMATE)
     call fftw_execute_r2r(plan,fIn,fOutReal)
     call fftw_destroy_plan(plan)
     fIn=aimag(svar(1,:,ilayer))
     plan=fftw_plan_r2r_1d(ny,fIn,fOutImag,FFTW_REDFT10,FFTW_ESTIMATE)
     call fftw_execute_r2r(plan,fIn,fOutImag)
     call fftw_destroy_plan(plan)
     do j=1,ny
        svar(1,j,ilayer)=complex(fOutReal(j),fOutImag(j))     
     end do

     ! Sine transform for the i>1 component in y (beware real array needed in input)
     do i=2,nx
        fIn=real(svar(i,:,ilayer))
        plan=fftw_plan_r2r_1d(ny,fIn,fOutReal,FFTW_RODFT10,FFTW_ESTIMATE)
        call fftw_execute_r2r(plan,fIn,fOutReal)
        call fftw_destroy_plan(plan)
        fIn=aimag(svar(i,:,ilayer))
        plan=fftw_plan_r2r_1d(ny,fIn,fOutImag,FFTW_RODFT10,FFTW_ESTIMATE)
        call fftw_execute_r2r(plan,fIn,fOutImag)
        call fftw_destroy_plan(plan)
        do j=1,ny
           svar(i,j,ilayer)=complex(fOutReal(j),fOutImag(j))     
        end do
     end do
     deallocate(fIn,fOutReal,fOutImag)

  end do

  return
end subroutine var2svar
!=================================================================================
!> Compute var in real space from svar in spactral space (from psi-like variables)
!=================================================================================
subroutine svar2var(svar,var,nx,ny,nlayers)
  use, intrinsic :: iso_c_binding
  use precision
  implicit none
  include "fftw3.f03"

  integer,intent(in) :: nx,ny,nlayers
  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(in) :: svar
  real(dp), dimension(0:nx+1,0:ny+1,1:nlayers), intent(out) :: var

  integer :: ilayer,i,j
  real(dp) :: scale
  complex(dp), dimension(1:nx) :: varIn,varOut
  complex(dp), dimension(1:nx,1:ny,1:nlayers) :: tmpsvar
  real(dp), dimension(:),allocatable :: fIn,fOutReal,fOutImag
  type(C_PTR) :: plan

 ! Loop over layers
  do ilayer=1,nlayers

     ! Cosine transform for the i=1 component in y (beware real array needed in input)
     allocate(fIn(ny),fOutReal(ny),fOutImag(ny))
     fIn=real(svar(1,:,ilayer))
     plan=fftw_plan_r2r_1d(ny,fIn,fOutReal,FFTW_REDFT01,FFTW_ESTIMATE)
     call fftw_execute_r2r(plan,fIn,fOutReal)
     call fftw_destroy_plan(plan)
     fIn=aimag(svar(1,:,ilayer))
     plan=fftw_plan_r2r_1d(ny,fIn,fOutImag,FFTW_REDFT01,FFTW_ESTIMATE)
     call fftw_execute_r2r(plan,fIn,fOutImag)
     call fftw_destroy_plan(plan)
     do j=1,ny
        tmpsvar(1,j,ilayer)=complex(fOutReal(j),fOutImag(j))     
     end do

     ! Sine transform for the i>1 component in y (beware real array needed in input)
     do i=2,nx
        fIn=real(svar(i,:,ilayer))
        plan=fftw_plan_r2r_1d(ny,fIn,fOutReal,FFTW_RODFT01,FFTW_ESTIMATE)
        call fftw_execute_r2r(plan,fIn,fOutReal)
        call fftw_destroy_plan(plan)
        fIn=aimag(svar(i,:,ilayer))
        plan=fftw_plan_r2r_1d(ny,fIn,fOutImag,FFTW_RODFT01,FFTW_ESTIMATE)
        call fftw_execute_r2r(plan,fIn,fOutImag)
        call fftw_destroy_plan(plan)
        do j=1,ny
           tmpsvar(i,j,ilayer)=complex(fOutReal(j),fOutImag(j))     
        end do
     end do
     deallocate(fIn,fOutReal,fOutImag)

     ! Complex transform data in x using complex array
     do j=1,ny
        varIn=tmpsvar(1:nx,j,ilayer) ! Convert input array to complex
        plan=fftw_plan_dft_1d(nx,varIn,varOut,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(plan,varIn,varOut)
        call fftw_destroy_plan(plan)
        var(1:nx,j,ilayer)=real(varOut)
     end do

  end do

  scale=2.d0*nx*ny
  var=var/scale

   return
end subroutine svar2var
!===============================================================================
!> Compute su from spsi
!===============================================================================
subroutine spsi2su(spsi,su,nx,ny,nlayers,scale)
  use, intrinsic :: iso_c_binding
  use precision
  implicit none
  include "fftw3.f03"

  integer,intent(in) :: nx,ny,nlayers
  real(dp),intent(in) :: scale
  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(in) :: spsi
  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(out) :: su

  integer :: j

  su = - spsi ! minus sign to account for -dpsidy

  !Manipulate Fourier coefficients to get derivatives for kx=1 mode
  do j=1,ny-1
     su(1,j,:)=-j/2.*su(1,j+1,:)
  end do
  su(1,ny,:)=0.d0

  !Manipulate Fourier coefficients to get derivatives for kx>1 modes
  do j=ny,2,-1
     su(2:nx,j,:)=+(j-1)/2.*su(2:nx,j-1,:)
  end do
  su(2:nx,1,:)=0.d0

  su = scale * su
  
  return
end subroutine spsi2su
!===============================================================================
!> Compute sv from spsi
!===============================================================================
subroutine spsi2sv(spsi,sv,nx,ny,nlayers,scale)
  use, intrinsic :: iso_c_binding
  use precision
  implicit none
  include "fftw3.f03"

  integer,intent(in) :: nx,ny,nlayers
  real(dp),intent(in) :: scale
  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(in) :: spsi
  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(out) :: sv

  integer :: i

  sv = + spsi ! plus sign to account for +dpsidx

  !Manipulate Fourier coefficients to get derivatives
  do i=1,nx/2
     sv(i,:,:)=+complex(0.d0,1.d0)*(i-1)*sv(i,:,:)
  end do
  do i=nx/2+1,nx
     sv(i,:,:)=+complex(0.d0,1.d0)*(i-nx-1)*sv(i,:,:)
  end do
  
  sv = scale * sv
  
  return
end subroutine spsi2sv
!===============================================================================
!> Compute u (zonal vel in real space) from su (zonal vel in spectral space)
!===============================================================================
subroutine su2u(su,u,nx,ny,nlayers)
  use, intrinsic :: iso_c_binding
  use precision
  implicit none
  include "fftw3.f03"

  integer,intent(in) :: nx,ny,nlayers
  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(in) :: su
  real(dp), dimension(0:nx+1,0:ny+1,1:nlayers), intent(out) :: u

  integer :: ilayer,i,j
  real(dp) :: scale
  complex(dp), dimension(1:nx) :: varIn,varOut
  complex(dp), dimension(1:nx,1:ny,1:nlayers) :: tmpsvar
  real(dp), dimension(:),allocatable :: fIn,fOutReal,fOutImag
  type(C_PTR) :: plan

 ! Loop over layers
  do ilayer=1,nlayers

     ! Cosine transform for the i=1 component in y (beware real array needed in input)
     allocate(fIn(ny),fOutReal(ny),fOutImag(ny))
     fIn=real(su(1,:,ilayer))
     plan=fftw_plan_r2r_1d(ny,fIn,fOutReal,FFTW_RODFT01,FFTW_ESTIMATE)
     call fftw_execute_r2r(plan,fIn,fOutReal)
     call fftw_destroy_plan(plan)
     fIn=aimag(su(1,:,ilayer))
     plan=fftw_plan_r2r_1d(ny,fIn,fOutImag,FFTW_RODFT01,FFTW_ESTIMATE)
     call fftw_execute_r2r(plan,fIn,fOutImag)
     call fftw_destroy_plan(plan)
     do j=1,ny
        tmpsvar(1,j,ilayer)=complex(fOutReal(j),fOutImag(j))     
     end do

     ! Sine transform for the i>1 component in y (beware real array needed in input)
     do i=2,nx
        fIn=real(su(i,:,ilayer))
        plan=fftw_plan_r2r_1d(ny,fIn,fOutReal,FFTW_REDFT01,FFTW_ESTIMATE)
        call fftw_execute_r2r(plan,fIn,fOutReal)
        call fftw_destroy_plan(plan)
        fIn=aimag(su(i,:,ilayer))
        plan=fftw_plan_r2r_1d(ny,fIn,fOutImag,FFTW_REDFT01,FFTW_ESTIMATE)
        call fftw_execute_r2r(plan,fIn,fOutImag)
        call fftw_destroy_plan(plan)
        do j=1,ny
           tmpsvar(i,j,ilayer)=complex(fOutReal(j),fOutImag(j))     
        end do
     end do
     deallocate(fIn,fOutReal,fOutImag)

     ! Complex transform data in x using complex array
     do j=1,ny
        varIn=tmpsvar(1:nx,j,ilayer) ! Convert input array to complex
        plan=fftw_plan_dft_1d(nx,varIn,varOut,FFTW_BACKWARD,FFTW_ESTIMATE)
        call fftw_execute_dft(plan,varIn,varOut)
        call fftw_destroy_plan(plan)
        u(1:nx,j,ilayer)=real(varOut)
     end do

  end do

  scale=2.d0*nx*ny
  u=u/scale

  return
end subroutine su2u
!===============================================================================
!> Compute derivative of f using FFTW based on real transform
!===============================================================================
subroutine derivetypeII(svar,nx,ny,nlayers,scale,dir)
  use precision
  implicit none

  integer :: nx,ny,nlayers,dir
  real(dp) :: scale
  complex(dp), dimension(nx,ny,nlayers) :: svar

  integer :: i,j,ilayer

  if (dir==1) then

     !Manipulate Fourier coefficients to get derivatives
     do i=1,nx/2
        svar(i,:,:)=+complex(0.d0,1.d0)*(i-1)*svar(i,:,:)
     end do
     do i=nx/2+1,nx
        svar(i,:,:)=+complex(0.d0,1.d0)*(i-nx-1)*svar(i,:,:)
     end do

  else if (dir==2) then
  
     !Manipulate Fourier coefficients to get derivatives for kx=1 mode
     do j=1,ny-1
        svar(1,j,:)=-j/2.*svar(1,j+1,:)
     end do
     svar(1,ny,:)=0.d0
     
     !Manipulate Fourier coefficients to get derivatives for kx>1 modes
     do j=ny,2,-1
        svar(2:nx,j,:)=+(j-1)/2.*svar(2:nx,j-1,:)
     end do
     svar(2:nx,1,:)=0.d0

  endif

  svar = scale * svar

  return
end subroutine derivetypeII
!===============================================================================
!> Compute Laplacian of svar in Fourier space
!===============================================================================
subroutine delSq(svar,nx,ny,nlayers)
  use params , only : kx2,ky2
  use precision
  implicit none

  integer :: nx,ny,nlayers
  complex(dp), dimension(nx,ny,nlayers) :: svar

  integer :: ilayer

  ! Compute delSq(svar)
  do ilayer=1,nlayers
     svar(:,:,ilayer) = - ( kx2 + ky2 ) * svar(:,:,ilayer)
  end do

  return
end subroutine delSq
