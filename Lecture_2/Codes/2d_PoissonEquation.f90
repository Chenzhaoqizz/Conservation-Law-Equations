!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!  Solve the 2D Poisson equation with various iteration methods 
!  in uniform and ununiform grid 
!==================================================================================
Program Poisson_Equation_2d
    use Mod_Grid
    use Mod_Poisson2d
    implicit none
    integer :: i, j, k, int_time1, int_time2

    call Grid_Generation
    
    allocate( phi(Nx, Ny), phio(Nx, Ny) )
    allocate( phi1d(NIJ), phio1d(NIJ) )
    allocate( AW(NIJ), AE(NIJ), AS(NIJ), AN(NIJ), P(NIJ) )

    do j = js + 1, je - 1 
        do i = is + 1, ie - 1
            k = j + Ny * (i-1)
            AW(k) = 2.d0 / ( (x(i+1)-x(i-1)) * (x(i)-x(i-1)) )
            AE(k) = 2.d0 / ( (x(i+1)-x(i-1)) * (x(i+1)-x(i)) )
            AS(k) = 2.d0 / ( (y(j+1)-y(j-1)) * (y(j)-y(j-1)) )
            AN(k) = 2.d0 / ( (y(j+1)-y(j-1)) * (y(j+1)-y(j)) )
            P(k)  = - ( AW(k) + AE(k) + AS(k) + AN(k) )
        enddo
    enddo

    call Setup_Parameters
    call Setup_Boundary_Conditions

    call system_clock(int_time1)

    if ( IterScheme == 0 ) then
        call Jacobi_Iter
    elseif ( IterScheme == 1 ) then
        call GS_Iter
    elseif ( IterScheme == 2 ) then 
        call SOR_Iter
    elseif ( IterScheme == 3 ) then 
        ! call ADI_Iter
    elseif ( IterScheme == 4 ) then 
        call SIP_Iter
    elseif ( IterScheme == 5 ) then 
        ! Only for uniform grid 
        call MSD_Iter 
    elseif ( IterScheme == 6 ) then 
        ! Only for uniform grid 
        call CG_Iter
    elseif ( IterScheme == 7 ) then 
        call CGS_Iter
    elseif ( IterScheme == 8 ) then 
        call ICCG_Iter
    elseif ( IterScheme == 9 ) then 
        call CGStab_Iter
    elseif ( IterScheme == 10 ) then 
        call BiCGStab_Iter
    else
        write(*,*) "No sutiable Iteration method, select in 0 - 9." 
    endif 

    call system_clock(int_time2)
    write(*,100) float( int_time2 - int_time1 ) / 1000.d0 
    100 Format('Execution time = ', 3X, 1E15.8, 's. ') 

    call Output_Result

    stop
End Program Poisson_Equation_2d
!----------------------------------------------------------------------------------
