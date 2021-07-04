module lastsaturation_tracer_mod
! ==================================================================================
use constants_mod, only : deg_to_rad
!==================================================================================
implicit none
private
!==================================================================================
! version information 
character(len=128) :: version='$Id: lastsaturation_tracer.f90 $'
character(len=128) :: tag='homemade'
!==================================================================================
! public interfaces
public :: lastsaturation_tracer 
!================================================================================== integer, allocatable, dimension(:) ::
integer :: nlat = 8
integer :: npres = 7
!q_index: 5:60
!==================================================================================
character(len=10), parameter :: mod_name = 'ls_tracer'
real :: missing_value = -999.
contains
! ==================================================================================    Atm(1)%q(isc:iec,jsc:jec,:,5:60) )
subroutine lastsaturation_tracer (is, ie, js, je, num_levels, &
        agrid, pfull, cond_dt_qg, q )

integer, intent(in)                 :: is, ie, js, je, num_levels
real, intent(in), dimension(:,:,:)  :: agrid
real, intent(in), dimension(:,:,:)  :: pfull,  cond_dt_qg
real, intent(inout), dimension(:,:,:,:)  :: q !(1:56)

real :: lat_tl
integer :: iLat, iPres, i, j, k
!master = (gid == masterproc)

do i=1,size(pfull,1)
   do j=1,size(pfull,2)
      do k=1,size(pfull,3)
         if ( (cond_dt_qg(i,j,k).lt.0) ) then
!reset tracers in other cells
            q(i,j,k, :) = 0.
!find the tidally-locked latitude
            lat_tl = acos(cos(agrid(i,j,2))*&
                cos(agrid(i,j,1)-270.*deg_to_rad ) ) / deg_to_rad !0:180
            iLat = int(lat_tl / 180.*8) +1 
!            if (iLat>8) iLat=8
!find the pressure
            if (pfull(i,j,k)>1e4 .and. pfull(i,j,k)<7e4) then
               iPres = int(pfull(i,j,k)/1e4)  
               q(i,j,k, (iLat-1)*7+iPres ) = 1.
            else  if (pfull(i,j,k) .ge. 7e4) then
               iPres = 7
               q(i,j,k, (iLat-1)*7+iPres ) = 1.
            endif

         endif
      enddo
   enddo
enddo   

return
end subroutine lastsaturation_tracer

! ==================================================================================

                                                                      
! ==================================================================================


end module lastsaturation_tracer_mod


