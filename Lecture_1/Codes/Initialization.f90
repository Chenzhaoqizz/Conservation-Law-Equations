!==================================================================================
! subroutine that set up the computational domain 
!----------------------------------------------------------------------------------
subroutine SetupDomain
    use Mod_UAD1d
    implicit none

    Nx = 811; 
    is = 1; ie = Nx; 
    x1 = 0.d0; x2 = 9.d0; 
    Length = x2 - x1;
    hx = Length / (Nx-1); 

    allocate( f(Nx), fo(Nx) ); 
    
    return
end subroutine SetupDomain

!==================================================================================
! subroutine that set up the initial condition 
!----------------------------------------------------------------------------------
subroutine SetupIC
    use Mod_UAD1d
    implicit none
    integer :: i
    real(8) :: x

    U = 1.d0; 
    D = 0.05d0; 
    CFL = 0.5d0; 

    EndTime = 2.5; time = 0.0; 
    
    ! dt = CFL*min( 0.5d0*hx*hx/D, 2.d0*D/(U*U) );
	dt = 0.00125

    tprint = 0.0; dtprint = 0.1; nprint = 0;

    MaxStep = 300000; iTimeStep = 0;

    TimeScheme = 3; ! 0 for FTCS, conditional stable
                    ! 1 for BTCS, unconditional stable
                    ! 2 for CN,   unconditional stable
                    ! 3 for RK3.  conditional stable
    out_path = 'out';

    do i = is, ie
        x = (i-1)*hx
        f(i) = exp( - (x-1.d0)**2 / (D*(4.d0*time+1.d0)) )
    enddo

    fo = f;

    return
end subroutine SetupIC

!==================================================================================
! subroutine that set up the boundary condition 
!----------------------------------------------------------------------------------
subroutine SetupBC
    use Mod_UAD1d
    implicit none
    real(8) :: c1, c2

    c1 = sqrt(4.d0*time + 1.d0)
    c2 = D*c1**2
    
    f(is) = exp( -(1.d0+U*time)**2/c2 ) / c1
    f(ie) = exp( -(8.d0-U*time)**2/c2 ) / c1

    return
end subroutine SetupBC
!----------------------------------------------------------------------------------

