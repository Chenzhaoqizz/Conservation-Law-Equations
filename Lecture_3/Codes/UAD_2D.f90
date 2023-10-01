!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!   Solve the 2D unsteady advection diffusion equation with various time scheme
!==================================================================================
Program Unsteady_Advection_Diffusion_2D
    use Mod_Grid
    use Mod_UAD2d
    implicit none

    pi = 4.d0 * atan(1.d0);
    
    TestCase = 0;   ! 0 for vortex
                    ! 1 for Gauss pluse, some bugs not fixed

    call SetupDomain
    call SetupIC
    call SetupBC
    call SetupVelocity

    do while ( time<=EndTime .and. iTimeStep<=MaxStep )

        call Error_Func
        write(*,100) iTimeStep, time, dt, Error
        100 Format('iStep = ', 1I6, ', Time =', 1f12.8, ', dt =', 1f12.8, ', Error =', 1E15.8)

        if ( (time+dt) >= tprint ) then 
            call Output_Result
            tprint = tprint + dtprint
            nprint = nprint +1
        endif

        if ( TimeScheme == 0 ) then 
            call RK2_Scheme
        elseif ( TimeScheme == 1 ) then
            call RK3_Scheme
        elseif ( TimeScheme == 2 ) then
            call FTCS_Scheme
        elseif ( TimeScheme == 3 ) then
            ! call CN_Scheme
        else
            write(*,*) "No suitable time scheme, select in 0, 1, 2, 3. " 
            stop
        endif

        time = time + dt
        iTimeStep = iTimeStep + 1
        call SetupBC
        call SetupVelocity
    enddo

    stop
End Program Unsteady_Advection_Diffusion_2D
!==================================================================================