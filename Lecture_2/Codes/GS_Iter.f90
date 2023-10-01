!==================================================================================
! subroutine that uses Gauss Seidile iteration 
! to solve the Poisson equation
!----------------------------------------------------------------------------------
subroutine GS_Iter
    use Mod_Grid
    use Mod_Poisson2d
    implicit none
    integer :: i, j, k
    real(8) :: f, a, b, c, d, e, R2

    open(unit=16,file='Residual_GS.dat',status='unknown')

    ! Begin iteration 
    do while ( nIter <= MaxIter .and. Error >= MaxError )

        Error = 0.d0
        phio1d = phi1d
    
        ! Iteration for inner grid point
        do i = is + 1, ie - 1
            do j = js + 1, je - 1 
                k = j + Ny * (i-1)
                a = AW(k) * phi1d(k-Ny)
                b = AE(k) * phi1d(k+Ny)
                c = AS(k) * phi1d(k-1)
                d = AN(k) * phi1d(k+1)
                e = P(k)
                f = - 4.d0 * sin( x(i) - y(j) ) * exp( x(i) - y(j) )
                
                phi1d(k) = ( f - a - b - c - d ) / e
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

        Error = Error / dble( NIJ )
        nIter = nIter + 1
        write(*,100) nIter, Error 
        write(16,110) nIter, Error
        100 Format('nIter = ', 1I6, 4X, 'Error = ', 1E15.8)
        110 Format(1I6, 2X, 1E15.8)
    enddo 

    close(unit=16) ! close file after writen


end subroutine GS_Iter
!----------------------------------------------------------------------------------