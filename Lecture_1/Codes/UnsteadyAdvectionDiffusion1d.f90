!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!   Solve the 1D unsteady advection diffusion equation with various time scheme
!==================================================================================
Program Unsteady_Advection_Diffusion_1d
    use Mod_UAD1d
    implicit none

    call SetupDomain
    call SetupIC
    call SetupBC

    do while (time<=EndTime .and. iTimeStep<MaxStep)

        if ( (time+dt) >= tprint ) then 
            call Output_Result
            tprint = tprint + dtprint
            nprint = nprint + 1
        endif

        if ( TimeScheme == 0 ) then 
            call FTCS_Scheme
        elseif ( TimeScheme == 1 ) then 
            call BTCS_Scheme
        elseif ( TimeScheme == 2 ) then 
            call CNCS_Scheme
        elseif ( TimeScheme == 3 ) then 
            call RK3_Scheme
        else
            write(*,*) "No sutiable timescheme, select in 0, 1, 2, 3. "
            stop
        endif 

        time = time+dt; 
        iTimeStep = iTimeStep+1; 
    enddo
    
    call Error_Func

    stop
End Program Unsteady_Advection_Diffusion_1d
!----------------------------------------------------------------------------------
