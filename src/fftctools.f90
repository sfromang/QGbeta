!===============================================================================
!> Compute FFT forward transform of array var into array svar
!===============================================================================
subroutine var2spec(var,svar,nx,ny,nlayers)
  use, intrinsic :: iso_c_binding
  use precision
  implicit none
  include "fftw3.f03"

  integer,intent(in) :: nx,ny,nlayers
  real(dp), dimension(0:nx+1,0:ny+1,1:nlayers), intent(in) :: var
  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(out) :: svar

  integer :: ilayer
  complex(dp), dimension(1:nx,1:ny) :: varIn,varOut
  type(C_PTR) :: plan

  ! Loop over layers
  do ilayer=1,nlayers

     ! Convert input array to complex
     varIn=var(1:nx,1:ny,ilayer)

     ! Make a plan for the FFT, and forward transform the data.
     plan=fftw_plan_dft_2d(ny,nx,varIn,varOut,FFTW_FORWARD,FFTW_ESTIMATE)
     call fftw_execute_dft(plan,varIn,varOut)
     call fftw_destroy_plan(plan)

     ! Store in final output array
     svar(:,:,ilayer)=varOut

  end do

  return
end subroutine var2spec
!===============================================================================
!> Compute FFT backward transform of array svar into array var
!===============================================================================
subroutine spec2var(svar,var,nx,ny,nlayers)
  use, intrinsic :: iso_c_binding
  use precision
  implicit none
  include "fftw3.f03"

  integer,intent(in) :: nx,ny,nlayers
  complex(dp), dimension(1:nx,1:ny,1:nlayers), intent(in) :: svar
  real(dp), dimension(0:nx+1,0:ny+1,1:nlayers), intent(out) :: var

  integer :: ilayer
  complex(dp), dimension(1:nx,1:ny) :: varIn,varOut
  type(C_PTR) :: plan

  ! Loop over layers
  do ilayer=1,nlayers

     ! Convert input array to complex
     varIn=svar(1:nx,1:ny,ilayer)

     ! Make a plan for the FFT, and forward transform the data.
     plan=fftw_plan_dft_2d(ny,nx,varIn,varOut,FFTW_BACKWARD,FFTW_ESTIMATE)
     call fftw_execute_dft(plan,varIn,varOut)
     call fftw_destroy_plan(plan)

     ! Store in final output array
     var(1:nx,1:ny,ilayer)=real(varOut(1:nx,1:ny))/nx/ny

  end do

  return
end subroutine spec2var
!===============================================================================
!> Compute derivative of f using FFTW based on complex transform
!===============================================================================
subroutine derive(svar,nx,ny,nlayers,dir)
  use precision
  use params , only : kx,ky
  use, intrinsic :: iso_c_binding 
  implicit none
  include "fftw3.f03"

  integer :: nx,ny,nlayers,dir
  complex(dp), dimension(nx,ny,nlayers) :: svar

  complex(dp), dimension(nx,ny) :: varIn,varOut
  type(C_PTR) :: plan

  integer :: i,j,ilayer

  do ilayer=1,nlayers

     ! Store variable in temporary array
     varIn=svar(:,:,ilayer)

     ! Make a plan for the FFT, and forward transform the data.
     plan=fftw_plan_dft_2d(ny,nx,varIn,varOut,FFTW_FORWARD,FFTW_ESTIMATE)
     call fftw_execute_dft(plan,varIn,varOut)
     call fftw_destroy_plan(plan)

     !Manipulate Fourier coefficients to get derivatives
     if (dir==1) then
        do i=1,nx
           varOut(i,:)=complex(0.d0,1.d0)*kx(i,:)*varOut(i,:)
        end do
     endif
     if (dir==2) then
        do j=1,ny
           varOut(:,j)=complex(0.d0,1.d0)*ky(:,j)*varOut(:,j)
        end do
     endif

     ! Overwrite initial array
     svar(:,:,ilayer)=varOut

   end do

  return
end subroutine derive
