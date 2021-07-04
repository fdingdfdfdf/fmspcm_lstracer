program lstracer_auto_grid
  
  implicit none

  !  -------------------------------------------------------------
  !  Purpose: create a grid of LBL code runs, and submit them all
  !  Authour: Robin Wordsworth (2016)
  !  -------------------------------------------------------------

  integer i, iX, iY,iZ !X:ps0, Y:tau_dry, Z:omega

  integer, parameter :: nX = 1 
  integer, parameter :: nY = 1
  integer, parameter :: nZ = 2

  character(256) temp
  character(64) dname
  character(3) Xlab
  character(3) Ylab
  character(3) Zlab

  character(128) fms_exp  
  character(128) fms_tmp

  real(8) X_ar(1:nX)
  real(8) Y_ar(1:nY)
  real(8) Z_ar(1:nZ)

  !--------------------------------------------------------
  fms_exp='/n/home01/fdingdfdfdf/cubed_sphere_lbl/exp/'
  fms_tmp='/n/holyscratch01/wordsworth_lab/fding/fms_tmp/'
  ! define grid 
  ! surface pressure in Pa
  X_ar(1) = 1e5

  ! tau_dry     
  Y_ar(1) = 1.2
  !Y_ar(2) = 0.5
  !Y_ar(3) = 1.0
  !Y_ar(4) = 2.0
  !Y_ar(5) = 3.0

  !roration rate
  Z_ar(1) = 1.454e-6  
  Z_ar(2) = 7.272e-6

  !open(14,file='T_ar.out')
  !write(14,*) X_ar
  !close(14)

  !open(15,file='fCO2_ar.out')
  !write(15,*) Y_ar
  !close(15)

  !--------------------------------------------------------
  ! loop over variables

  print*,'Setting up the lstracer runs'

!  ps_loop : do iX=1,nX
  iX=1
  tau_loop :  do iY=1,nY
  omega_loop:    do iZ=1,nZ
  
        print*,'--- Input ---'
        print*,'surface pressure = ',X_ar(iX) /1e5,' bar'
        print*,'tau_dry          = ',Y_ar(iY)
        print*,'omega            = ',Z_ar(iZ),' rad/s'

        write(Xlab,'(i2)') iX
        write(Ylab,'(i2)') iY
        write(Zlab,'(i2)') iZ

        !dname = 'run_T'//trim(adjustl(Xlab))//'_C'//trim(adjustl(Ylab))
        dname = 'run_ps1bar_tau'//trim(adjustl(Ylab))//'_rot'//trim(adjustl(Zlab))

        call system('cp -rf '//trim(fms_exp)//'slow_ps5_tau3 '//trim(fms_exp)//trim(dname))

        ! update surface temperature 
        open(14,file=trim(fms_exp)//'slow_ps5_tau3/run/runscript')
        open(15,file=trim(fms_exp)//trim(dname)//'/run/runscript')
        do i=1,417
           read(14,'(a256)') temp
           if(i.eq.38)then
              write(15,'(a14,e14.7)') 'set omega   = ',Z_ar(iZ)
           elseif(i.eq.39)then
              write(15,'(a14,e14.7)') 'set ps0     = ',X_ar(iX)
           elseif(i.eq.40)then
              write(15,'(a14,e14.7)') 'set tau_dry = ',Y_ar(iY)
           else
              write(15,'(a)') trim(temp)
              !write(15,*) trim(adjustl(temp))
           endif
        enddo
        close(14)
        close(15)

!        call system('cp -rf '//trim(fms_tmp)//'slow_ps5_tau3 '//trim(fms_tmp)//trim(dname))

        print*,'Submitting run in directory: ',trim(dname)
        call system('cd '//trim(fms_exp)//trim(dname)//'/run && sbatch runscript')

  end do omega_loop
  end do tau_loop

end program lstracer_auto_grid

