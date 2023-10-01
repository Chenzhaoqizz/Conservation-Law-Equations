!======================================================================================
! sunroutine to initialize physical problems
!--------------------------------------------------------------------------------------
Subroutine Initialize
    use Mod_grid
    use Mod_flow
    implicit none

    allocate( u(Nx,Ny), v(Nx,Ny), phi(Nx,Ny), phio(Nx,Ny) )
    u = 0; v = 0; phi = 0
    is = 1; ie = Nx
    js = 1; je = Ny
    beta = 2.0d0

    select case(iproblem)
    case(1)
        call Ini_ShearFlow
    case(2)
        call Ini_ZaleskiRotation
    end select

    write(*,999)' delta_x = ',hxi,' delta_y = ',hyi,' dt = ',delta_t
    999 Format( A,1E16.9,2X,A,1E16.9,2X,A,1D16.9 )   

End Subroutine Initialize 
!======================================================================================
! initialize shear flow problem
!--------------------------------------------------------------------------------------
subroutine Ini_ShearFlow
    use Mod_grid
    use Mod_flow
    implicit none
    integer :: i, j
    real(8) :: xc, yc, rad, dist

    xLength = pi; yLength = pi;
    hxi = xLength / dble(ie-is)
    hyi = yLength / dble(je-js)
    beta = beta / min( hxi, hyi )

    xc = 0.2 * ( pi + 1.0 )
    yc = 0.5 * pi
    rad = 0.2 * pi

    ! Deformed interface in a shearing flow, Sec-4.2.2 
    do i = is, ie
        do j = js, je 
            dist = dSqrt( ( i * hxi - xc )**2 + ( j * hyi - yc )**2 )
            dist = rad - dist
            phi(i,j) = 0.5d0*( 1.d0 + dtanh( beta * dist ) ) 
            u(i,j) =  dsin(i*hxi)*dcos(j*hyi)
            v(i,j) = -dcos(i*hxi)*dsin(j*hyi) 
        enddo
    enddo
    phio = phi

end subroutine Ini_ShearFlow
!======================================================================================
! initialize Zaleski rotation problem
!--------------------------------------------------------------------------------------
Subroutine Ini_ZaleskiRotation
    use Mod_grid
    use Mod_flow
    implicit none
    integer :: i, j
    real(8) :: xc, yc, rad, dist

    xLength = 1.0; yLength = 1.0;
    hxi = xLength / dble(ie-is)
    hyi = yLength / dble(je-js)
    beta = beta / min( hxi, hyi ) 

    xc = 0.5
    yc = 0.75
    rad = 0.15

    do i=is,ie 
        do j=js,je 
            dist = dSqrt( (i*hxi - xc)**2 + (j*hyi - yc)**2 )
            dist = rad - dist
            if ( (i*hxi >= 0.525) .or. (i*hxi<=0.475) .or. (j*hyi>=0.85) ) then
                phi(i,j) = 0.5 * ( 1.0 + dtanh( beta * dist ) )
            endif
            u(i,j) = 0.5 - j * hyi 
            v(i,j) = i * hxi - 0.5 
        enddo 
    enddo 
    phio = phi

end Subroutine Ini_ZaleskiRotation
