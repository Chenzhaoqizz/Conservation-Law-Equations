!==================================================================================
! subroutine that uses Stone's Strong implicit 
! method(SIP) to solve the Poisson equation
!----------------------------------------------------------------------------------
Subroutine SIP_Iter
    use Mod_Grid
    use Mod_Poisson2d
    implicit none
    integer :: i, j, k
    real(8) :: alpha, R2, a, b, c, d, e, f, p1, p2
    real(8), dimension(NIJ) :: LW, LS, UE, UN, LPR, Res

    open(unit=16,file='Residual_SIP.dat',status='unknown')

    ! Begin iteration 
    do while ( nIter <= MaxIter .and. Error >= MaxError )

        alpha = 0.975d0
        Res = 0
        R2 = 0
        LW = 0
        LS = 0
        UE = 0
        UN = 0
        LPR = 0
    
        ! Calculate elements of [L] and [U] matrices
        do i = is + 1, ie - 1
            do j = js + 1, je -1
                k = j + Ny * ( i-1 )
                LW(k) = AW(k) / ( 1.0 + alpha * UN(k-Ny) )
                LS(k) = AS(k) / ( 1.0 + alpha * UE(k-1) )
                p1 = alpha * LW(k) * UN(k-Ny)
                P2 = alpha * Ls(k) * UE(k-1)
                LPR(k) = 1.0 / ( P(k) + p1 + p2 - LW(k)*UE(k-Ny) - LS(k)*UN(k-1) )
                UN(k) = ( AN(k) - p1 ) * LPR(k)
                UE(k) = ( AE(k) - p2 ) * LPR(k) 
            enddo
        enddo
    
        ! Calculate residual vectors and make forward substitution
        do i = is + 1, ie - 1
            do j = js + 1, je -1
                k = j + Ny * ( i-1 )
                a =  P(k) * phi1d(k)
                b = AN(k) * phi1d(k+1)
                c = AS(k) * phi1d(k-1)
                d = AE(k) * phi1d(k+Ny)
                e = AW(k) * phi1d(k-Ny)
                f = - 4.d0 * sin( x(i) - y(j) ) * exp( x(i) - y(j) )
                Res(k) = f - a - b - c - d - e 
                R2 = R2 + Res(k) * Res(k) 
                Res(k) = LPR(k) * ( Res(k) - Ls(k)*Res(k-1) - Lw(k)*Res(k-Ny) )
            enddo
        enddo
    
        ! Make backward substitution and update solution
        do i = ie - 1, is + 1, -1
            do j = je -1, js + 1, -1
                k = j + Ny * ( i-1 )
                Res(k) = Res(k) - UN(j) * Res(k+1) - UE(k) * Res(k+Ny)
                phi1d(k) = phi1d(k) + Res(k) 
            enddo
        enddo
    
        ! Check convergence using L2 Norm
        Error = sqrt(R2) / dble(NIJ)
        nIter = nIter + 1
        write(*,100) nIter, Error 
        write(16,110) nIter, Error
        100 Format('nIter = ', 1I6, 4X, 'Error = ', 1E15.8)
        110 Format(1I6, 2X, 1E15.8)
    enddo 

    close(unit=16) ! close file after writen

End Subroutine SIP_Iter