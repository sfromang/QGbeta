program QGbeta
  use variables
  use params
  implicit none
  real(dp) :: tdump, thist

  ! Initialize flow
  time=0.d0
  call init
  tdump=time ; thist=time
  if (ndump==0) then
     call output
     call history
  endif
  ndump=ndump+1

  ! Loop on time & update variables
  do

     ! Calculate updated value for PV (q variables)
     call getUpdate

     ! Move to next time
     time=time+dt
     write (*,*) time,dt

     ! Write output
     if ( (time-dt .le. tdump+dtdump) .and. (time .gt. (tdump+dtdump)).and.(dtdump>0.d0)  ) then
        call output
        tdump = tdump + dtdump
        ndump = ndump + 1
     endif
     
     if ( (time-dt .le. thist+dthist) .and. (time .gt. (thist+dthist)).and.(dthist>0.d0)  ) then
        call history
        thist = thist + dthist
     endif

     ! Exit and terminate code
     if (time>tlim) exit

  end do

  write (*,*) 'Execution terminated.'

end program QGbeta

