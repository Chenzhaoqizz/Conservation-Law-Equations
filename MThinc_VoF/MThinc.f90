!----------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------
!====================================================================================================
!  Ref: An interface capturing method with a continuous function: The THINC 
!       method with multi-dimensional reconstruction,2012,JCP.  
!
!  Note: 1. Real(Kind=8)的数据类型为双精度实型，占据的储存空间为8个字节，有效位数为15-16位
!        2. Fortran自由格式中每行不要超过132个字符 
!====================================================================================================
!====================================================================================================
Program MThincVoF2d_2021
    use Mod_grid 
    use Mod_flow
    implicit none
    integer :: subStep, nPrintStep
    real(8) :: phi_ini

    pi = 4.0*atan(1.d0)
    phi_ini = 0.0
    time = 0.d0
    iStep = 0

    print*, "MThinc VoF program!"

    ! step-0: read some input parameters
    call read_input
    nprint = 0
    nPrintStep = int( 0.5 / delta_t )
    write(*,*) "Output results every ", nPrintStep, " steps!"

    ! step-1: generate uniform mesh and initialize the problem
    call Initialize
    phi_ini = sum(phi) * hxi * hyi
    print*, "phi total = ", phi_ini

    call output("Ini")

    allocate( NormVectorX(Nx,Ny), NormVectorY(Nx,Ny), Lxx(Nx,Ny), Lxy(Nx,Ny), Lyy(Nx,Ny) )
    NormVectorX = 0.0; NormVectorY = 0.0; Lxx = 0.0; Lxy = 0.0; Lyy = 0.0;

    allocate( flux_x(Nx, Ny), flux_y(Nx, Ny) )
    flux_x = 0.0; flux_y = 0.0;

    ! step-2: start time advance by direction splitting scheme
    do while ( time <= End_Time )

        phio = phi
        iStep = iStep + 1
        time = time + delta_t
        subStep = mod( iStep, 2 )

        ! splitting direction by subStep
        if ( subStep == 1 ) then

            ! subStep = 0 

            ! first advection in X direction
            call Normal_Vectors
            call NumericalFlux_X
            call Advection_X

            ! then advection in Y direction
            call Normal_Vectors
            call NumericalFlux_Y
            call Advection_Y
        else if ( subStep == 0 ) then

            ! subStep = 1 

            ! first advection in Y direction
            call Normal_Vectors
            call NumericalFlux_Y
            call Advection_Y

            ! then advection in X direction
            call Normal_Vectors
            call NumericalFlux_X
            call Advection_X 
        endif

        write(*,100) iStep, Time
        100 format( "iStep = ", 1I5, ", Time = ", 1ES16.8 ) 
        ! if output result
        if ( Mod(iStep, nPrintStep) == 0 ) then 
            nprint = nprint + 1
            call output("grd")
        endif 
        
    enddo

    print*, "Finish run MThinc VoF program!" 
    Stop
End Program 

