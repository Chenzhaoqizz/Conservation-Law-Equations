!==================================================================================
! subroutine that generate grid in the computational domain 
!----------------------------------------------------------------------------------
subroutine Grid_Generation
    use Mod_Grid
    implicit none
    integer :: i, j
    real(8) :: pi, Re, x_ratio, y_ratio, beta

    pi = 4.d0*atan(1.d0);
    Re = 1000.d0;
    Is_Uniform_Grid = .TRUE. 

    Nx = 257; 
    Ny = 257; 
    NIJ = Nx * Ny;

    allocate( x(Nx), y(Ny), deltaX(Nx), deltaY(Ny) )

    is = 1; ie = Nx; 
    js = 1; je = Ny;

    xLeft =  0.d0;  
    xRight = 3.d0; 
    yBottom = - 1.d0; 
    yTop = 2.d0;

    xLength = abs(xRight - xLeft);
    yLength = abs(yTop - yBottom);

    ! Genarate grid in X direction 
    do i = is, ie
        if ( Is_Uniform_Grid ) then
            deltaX(i) = xLength * float(i-1) / float( ie - is )
        else if ( .not.(Is_Uniform_Grid) .and. Re < 200 ) then
            ! Gause-Lobatto distribution 
            x_ratio = 0.5d0 
            beta = 2.d0 * ( float( i - is ) / float( ie - is ) - 0.5 )
            deltaX(i) = xLength * ( ( x_ratio * sin( pi * beta / 2.d0 ) + ( 1.d0 - x_ratio ) * beta + 1.d0 ) / 2.d0 )
        else if ( .not.(Is_Uniform_Grid) .and. Re >= 200 ) then
            ! Hyperbolic tangent transformation 
            x_ratio = 1.0 * log( Re )
            beta = 2.d0 * ( float( i - is ) / float( ie - is ) - 0.5 )
            deltaX(i) = xLength * ( tanh( x_ratio * beta ) / tanh( x_ratio ) + 1.d0 ) / 2.d0 
        endif
        x(i) = xLeft + deltaX(i)
		! write(*,*) "Delta_X = ", deltaX(i)-1.d0
    enddo

    ! Genarate grid in Y direction 
    do j = js, je
        if ( Is_Uniform_Grid ) then
            deltaY(j) = yLength * float(j-1) / float( je - js )
        else if ( .not.(Is_Uniform_Grid) .and. Re < 200 ) then
            ! Gause-Lobatto distribution 
            y_ratio = 0.5d0 ! y_ratio in [-1, 1]
            beta = 2.d0 * ( float( j - js ) / float( je - js ) - 0.5 )
            deltaY(j) = yLength * ( ( y_ratio * sin( pi * beta / 2.d0 ) + ( 1.d0 - y_ratio ) * beta + 1.d0 ) / 2.d0 )
        else if( .not.(Is_Uniform_Grid) .and. Re >= 200 ) then
            ! Hyperbolic tangent transformation 
            y_ratio = 1.0 * log( Re )
            beta = 2.d0 * ( float( j - js ) / float( je - js ) - 0.5 )
            deltaY(j) = yLength * ( tanh( y_ratio * beta ) / tanh( y_ratio ) + 1.d0 ) / 2.d0 
        endif
        y(j) = yBottom + deltaY(j)
    enddo

end subroutine Grid_Generation

!==================================================================================
! subroutine that setup some parameters 
!----------------------------------------------------------------------------------
subroutine Setup_Parameters
    use Mod_Poisson2d
    implicit none
    
    Error = 1.d0
    MaxError = 1.0E-6

    IterScheme = 10
    ! 0 for Jacobi iteration
    ! 1 for Gauss-Seidler iteration
    ! 2 for SOR iteration
    ! 3 for ADI iteration 
    ! 4 for Stone's strong implicit method (SIP)
    ! 5 (Uniform_Grid only) for method of steepest descent (MSD), 
    ! 6 (Uniform_Grid only) for conjugate gradient method (CG)
    ! 7 for conjugate gradient square method (CGS)
	! 8 for incomplete cholesky conjugate gradient method (ICCG)
	! 9 for CGStab
	! 10 for BiCGStab

    out_path = 'out'

    nIter = 0
    MaxIter = 1E5

    return
end subroutine Setup_Parameters

!==================================================================================
! subroutine that setup the boundary values by analytical solution
!----------------------------------------------------------------------------------
subroutine Setup_Boundary_Conditions
    use Mod_Grid
    use Mod_Poisson2d
    implicit none
    integer :: i, j, k
    real(8) :: co

    phi = 0.d0; phio = 0.d0;
    phi1d = 0.d0; phio1d = 0.d0;

    do j = js, je

        k = j + Ny*(is-1)
        co = cos( x(is) - y(j) ) * exp( x(is) - y(j) )
        phi(is, j) = co
        phi1d(k) = co

        k = j + Ny*(ie-1)
        co = cos( x(ie) - y(j) ) * exp( x(ie) - y(j) )
        phi(ie, j) = co
        phi1d(k) = co

    enddo

    do i = is, ie

        k = js + Ny*(i-1)
        co = cos( x(i) - y(js) ) * exp( x(i) - y(js) )
        phi(i, js) = co
        phi1d(k) = co
        
        k = je + Ny*(i-1)
        co = cos( x(i) - y(je) ) * exp( x(i) - y(je) )
        phi(i, je) = co
        phi1d(k) = co

    enddo

    return
end subroutine Setup_Boundary_Conditions

!==================================================================================
! function that calculates the analytical solution at a point
!----------------------------------------------------------------------------------
real(8) function analytical(x0, y0)
    implicit none
    real(8), intent(in) :: x0, y0

    analytical = cos( x0 - y0 ) * exp( x0 - y0 )

end function analytical
