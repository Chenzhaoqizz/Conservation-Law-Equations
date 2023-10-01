!======================================================================================
!        solve the Navier-Stokes equation in vorticity-stream function form
!--------------------------------------------------------------------------------------
program main
    use flow_var
    implicit none
    integer :: int_time1, int_time2

    nPrint = 0
    call read_input
    call alloc_data(Nx, Ny)
    call Initialize
    call system_clock(int_time1)

    ! Time advance
    do while ( Time <= End_Time )

        call Stream_Function
        call Vorticity_BC
        call Vorticity

        ! output result after iSave steps
        if ( mod( iStep, iSave ) == 0 ) then
            nPrint = nPrint + 1
            call output
        endif

        iStep = iStep + 1
        Time = Time + delta_t
        write(*,100) iStep, Time
        100 Format(" iStep = ", 1I6, ", Time = ", 1ES16.8)
    enddo
    
    call system_clock(int_time2)
    write(*,110) float( int_time2 - int_time1 ) / 1000.d0 
    110 Format('Finish run, execution time = ', 1X, 1ES15.8, 's. ') 

end program main
!--------------------------------------------------------------------------------------
