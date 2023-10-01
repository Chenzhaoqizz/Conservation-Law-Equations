!==================================================================================
! subroutine that uses conjugate gradient method (CG)
! to solve the Poisson equation, only suitable for 
! systematic positive define matrix 
!----------------------------------------------------------------------------------
subroutine CG_Iter
    use Mod_Grid
    use Mod_Poisson2d
    implicit none 
    integer :: i,j,k 
    real(8) :: f, a, b, c, d, e, R2, R20, alpha, beta 
    real(8), dimension(NIJ) :: Res, DRes 

    R2 = 0.d0
    R20 = 0.d0
    Res = 0.d0 
    DRes = 0.d0 
    beta = 0.d0 
    alpha = 0.d0 

    open(unit=16,file='Residual_CG.dat',status='unknown')

    ! Calculate initial residual vector 
    do i = is +1, ie-1
        do j = js+1, je-1
            k = j + Ny*(i-1) 
            a =  P(k) * phi1d(k)
            b = AN(k) * phi1d(k+1) 
            c = AS(k) * phi1d(k-1)
            d = AE(k) * phi1d(k+Ny)
            e = AW(k) * phi1d(k-Ny)
            f = -4.d0 * sin( x(i)-y(j) ) * exp( x(i) - y(j) ) 
            Res(k) = f - a - b - c - d - e 
            R2 = R2 + Res(k) * Res(k) 
        enddo
    enddo
    DRes = Res 
    R20 = R2 

    ! Begin iteration 
    do while ( nIter <= MaxIter .and. Error >= MaxError )

        ! Calculate the space step: alpha 
        R2 = 0
        do i = is +1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1) 
                R2 = R2 + Res(k) * Res(k) 
                Res(k) =  P(k)*DRes(k) & 
                    &  + AN(k)*DRes(k+1)  + AS(k)*DRes(k-1) & 
                    &  + AE(k)*DRes(k+Ny) + AW(k)*DRes(k-Ny) 
                alpha = alpha + DRes(k) * Res(k) 
            enddo
        enddo
        R20 = R2 
        alpha = R20 / (alpha + 1E-30) 

        ! Update the solution 
        phi1d(:) = phi1d(:) + alpha * DRes(:) 

        ! Calculate new residual vector using updated solution 
        R2 = 0.d0
        do i = is +1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1) 
                a =  P(k) * phi1d(k)
                b = AN(k) * phi1d(k+1) 
                c = AS(k) * phi1d(k-1)
                d = AE(k) * phi1d(k+Ny)
                e = AW(k) * phi1d(k-Ny)
                f = -4.d0 * sin( x(i)-y(j) ) * exp( x(i) - y(j) ) 
    
                Res(k) = f - a - b - c - d - e 
                R2 = R2 + Res(k) * Res(k) 
            enddo
        enddo

        beta = R2 / ( R20 + 1E-30 ) 

        ! Update the direction vector 
        DRes(:) = Res(:) + beta * DRes(:) 

        ! Check for convergence using L2 Norm 
        Error = Sqrt(R2) / dble(NIJ)
        nIter = nIter + 1
        write(*,100) nIter, Error 
        write(16,110) nIter, Error
        100 Format('nIter = ', 1I6, 4X, 'Error = ', 1E15.8)
        110 Format(1I6, 2X, 1E15.8)
    enddo 

    close(unit=16) ! close file after writen

end subroutine CG_Iter
