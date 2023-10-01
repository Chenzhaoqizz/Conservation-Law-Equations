!==================================================================================
!----------------------------------------------------------------------------------
subroutine SetupDomain
    use Mod_Grid
    use Mod_UAD2d
    implicit none
    integer :: i, j 

    Nx = 501; 
    Ny = 501;
    is = 1; ie = Nx; 
    js = 1; je = Ny;

    if ( TestCase == 0 ) then 
        x1 = -2.d0; x2 = 2.d0; 
        y1 = -2.d0; y2 = 2.d0;
    else if ( TestCase == 1 ) then 
        x1 = 0.d0; x2 = 9.d0; 
        y1 = 0.d0; y2 = 9.d0;
    endif

    xLength = x2 - x1;
    yLength = y2 - y1;

    hx = xLength / (ie - is);
    hy = yLength / (je - js);

    allocate(x(Nx), y(Ny));
    do i = is, ie 
        x(i) = hx * (i-1) + x1
    enddo

    do j = js, je
        y(j) = hy * (j-1) + y1
    enddo

    return
end subroutine SetupDomain

!==================================================================================
!----------------------------------------------------------------------------------
subroutine SetupIC
    use Mod_Grid
    use Mod_UAD2d
    implicit none

    integer :: i, j

    if ( TestCase == 0 ) then 
        D = 0.001d0; 
    else if ( TestCase == 1 ) then 
        D = 0.05d0; 
    endif

    CFL = 0.50d0; dt = 0.00025d0;
    EndTime = 2.5d0; time = 0.0d0; 

    MaxStep = 10000000; iTimeStep = 0;

    tprint = 0.d0; dtprint = 0.5d0; nprint = 0;

    TimeScheme = 0; ! 0 for RK-2, conditional stable
                    ! 1 for RK-3, conditional stable
                    ! 2 for FTCS. conditional stable
                    ! 3 for CN.   unconditional stable

    out_path = 'out';

    allocate( U(Nx, Ny), V(Nx, Ny), phi(Nx, Ny), phi_exat(Nx, Ny) );
    U = 0.d0; V = 0.d0; 
    
    ! Initialize value of phi by analytic solution here
    if ( TestCase == 0 ) then
        do j = js, je 
            do i = is, ie 
                phi(i,j) = cos( pi*x(i) ) * cos(pi*y(j))
            enddo 
        enddo
    else if ( TestCase == 1 ) then
        U = 1.0
        V = 1.0
        do j = js, je 
            do i = is, ie 
                phi(i,j) = exp( - ( (x(i)-1)**2 + (y(j)-1)**2 ) / D )
            enddo 
        enddo
    endif

    return 
end subroutine SetupIC

!==================================================================================
!     Setup boundary conditions for phi using analytic solution
!----------------------------------------------------------------------------------
subroutine SetupBC
    use Mod_Grid
    use Mod_UAD2d
    implicit none
    integer :: i, j
    real(8) :: t 

    if ( TestCase == 0 ) then
        t = -2.d0*D*pi*pi*time;
        do j = js, je 
            phi(is,j) = cos(pi*x(is)) * cos(pi*y(j)) * exp(t); ! Left BC
            phi(ie,j) = cos(pi*x(ie)) * cos(pi*y(j)) * exp(t); ! Right BC
        enddo 

        do i = is, ie 
            phi(i,js) = cos(pi*x(i)) * cos(pi*y(js)) * exp(t); ! Bottom BC
            phi(i,je) = cos(pi*x(i)) * cos(pi*y(je)) * exp(t); ! Top BC
        enddo
    else if ( TestCase == 1 ) then
        t = 4.d0 * time + 1.d0
        do j = js, je 
            phi(is,j) = (1.0/t)*exp( -( (x(is)-U(is,j)*time-1.0)**2 + (y(j)-V(is,j)*time-1.0)**2 ) / (D*t) ); ! Left BC
            phi(ie,j) = (1.0/t)*exp( -( (x(ie)-U(ie,j)*time-1.0)**2 + (y(j)-V(ie,j)*time-1.0)**2 ) / (D*t) ); ! Right BC
        enddo 

        do i = is, ie 
            phi(i,js) = (1.0/t)*exp( -( (x(i)-U(i,js)*time-1.0)**2 + (y(js)-V(i,js)*time-1.0)**2 ) / (D*t) ); ! Bottom BC
            phi(i,je) = (1.0/t)*exp( -( (x(i)-U(i,je)*time-1.0)**2 + (y(je)-V(i,je)*time-1.0)**2 ) / (D*t) ); ! Top BC
        enddo
    endif

    return 
end subroutine SetupBC

!==================================================================================
!     Setup velocity at time t 
!----------------------------------------------------------------------------------
subroutine SetupVelocity
    use Mod_Grid
    use Mod_UAD2d
    implicit none
    integer :: i, j 
    real(8) :: t, vel

    vel = 0.0
    if ( TestCase == 0 ) then
        t = -2.d0*D*pi*pi*time
        do i = is, ie 
            do j = js, je 
                U(i,j) = -pi * cos(pi*x(i)) * sin(pi*y(j)) * exp(t)
                V(i,j) =  pi * sin(pi*x(i)) * cos(pi*y(j)) * exp(t)
                vel = max(vel, sqrt( u(i,j)**2 + v(i,j)**2 ))
            enddo 
        enddo
    else if ( TestCase == 1 ) then
        U = 1.0
        V = 1.0
        vel = sqrt(2.d0)
    endif

    ! dt = max( CFL * max(hx, hy) / vel, dt ) 

    return
end subroutine SetupVelocity