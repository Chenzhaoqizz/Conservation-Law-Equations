!==================================================================================
! subroutine that uses SOR iteration 
! to solve the Poisson equation
!----------------------------------------------------------------------------------
subroutine SOR_Iter
    use Mod_Grid
    use Mod_Poisson2d
    implicit none 
    integer :: i, j, k
    real(8) :: f, beta, a, b, c, d, e, R2

    open(unit=16,file='Residual_SOR.dat',status='unknown')

    ! Begin iteration 
    do while ( nIter <= MaxIter .and. Error >= MaxError )

        phio1d = phi1d
        Error = 0.d0
        R2 = 0.d0
        beta = 1.75d0
    
        ! Iteration for inner grid point
        do i = is + 1, ie - 1
            do j = js + 1, je - 1 
                k = j + Ny * (i - 1)
                a = AW(k) * phi1d(k-Ny)
                b = AE(k) * phi1d(k+Ny)
                c = AS(k) * phi1d(k-1)
                d = AN(k) * phi1d(k+1)
                e = P(k)
                f = - 4.d0 * sin( x(i) - y(j) ) * exp( x(i) - y(j) )
                phi1d(k) = beta * ( f - a - b - c - d ) / e & 
                    & + ( 1.d0 - beta ) * phi1d(k)
                ! Error = Error + abs( phio1d(k) - phi1d(k) ) &
                !     & / abs( phi1d(k) )
            enddo
        enddo
    
        ! Calculate Residual using L2 Norm
        do i = is + 1, ie - 1
            do j = js + 1, je - 1 
                k = j + Ny * (i - 1)
                a = AW(k) * phi1d(k-Ny)
                b = AE(k) * phi1d(k+Ny)
                c = AS(k) * phi1d(k-1)
                d = AN(k) * phi1d(k+1)
                e = P(k) * phi1d(k)
                f = - 4.d0 * sin( x(i) - y(j) ) * exp( x(i) - y(j) )
                R2 = ( f - a - b - c - d - e )
                Error = Error + R2 * R2
            enddo
        enddo
    
        Error = sqrt(Error) / dble( NIJ )
        nIter = nIter + 1
        write(*,100) nIter, Error 
        write(16,110) nIter, Error
        100 Format('nIter = ', 1I6, 4X, 'Error = ', 1E15.8)
        110 Format(1I6, 2X, 1E15.8)
    enddo 

    close(unit=16) ! close file after writen


end subroutine SOR_Iter