!==================================================================================
! subroutine that uses Incomplete Cholesky Preconditioned  
! Conjugate Gradient method (ICCG) to solve the Poisson equation, 
! only suitable for systematic positive define matrix 
!----------------------------------------------------------------------------------
subroutine ICCG_Iter 
    use Mod_Grid
    use Mod_Poisson2d
    implicit none 
    integer :: i,j,k 
    real(8) :: f, a, b, c, d, e, R2, R20, alpha, beta 
    real(8) :: s0, sk
    real(8), dimension(NIJ) :: Res, Diag, Zk, Pk

    R2 = 0.d0
    R20 = 0.d0
    Res = 0.d0 
    Diag = 0.d0
    Zk = 0.d0
    Pk = 0.d0
    beta = 0.d0 
    alpha = 0.d0 
    s0 = 0.d0
    sk = 0.d0

    open(unit=16,file='Residual_ICCG.dat',status='unknown')

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
    R20 = R2 

    ! Preconditioning matrix diagonal
    do i = is + 1, ie -1 
        do j = js + 1, je -1
            k = j + Ny*(i-1)
            Diag(k) = 1.d0 / ( P(k) &
              - AW(k)**2 * Diag(k-Ny) &
              - AS(k)**2 * Diag(k-1 ) &
            )
        enddo
    enddo

    s0 = 1.0E20

    ! Begin iteration 
    do while ( nIter <= MaxIter .and. Error >= MaxError )
        ! Forward substitution
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                Zk(k) = (  Res(k) &
                  - AW(k) * Zk(k-Ny) &
                  - AS(k) * Zk(k-1)  &
                ) * Diag(k) 
            enddo
        enddo

        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                Zk(k) = Zk(k) / ( Diag(k) + 1.0E-32 )
            enddo
        enddo

        sk = 0.d0
        ! Backward substitution
        do i = ie-1, is+1, -1
            do j = je-1, js+1, -1
                k = j + Ny*(i-1)
                Zk(k) = Diag(k) * ( Zk(k) &
                  - AE(k) * Zk(k+Ny) &
                  - AN(k) * Zk(k+1)  &
                ) 
                sk = sk + Res(k) * Zk(k) 
            enddo
        enddo

        ! Calculate beta and new search vector
        beta = sk / (s0 + 1.0e-20)
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                Pk(k) = Zk(k) + beta * Pk(k)
            enddo
        enddo

        ! Calculate scalar product (Pk . A Pk) and alpha ( A Pk overwrites Zk ) 
        s0 = 0.d0
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                Zk(k) =  P(k) * Pk(k) &
                      + AE(k) * Pk(k+Ny) + AW(k) * Pk(k-Ny) &
                      + AN(k) * Pk(k+1)  + AS(k) * Pk(k-1) 
                s0 = s0 + Pk(k) * Zk(k)
            enddo
        enddo

        alpha = sk / (s0 + 1e-20)

        ! Calculate new residual amd update the solution
        R2 = 0.d0
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                phi1d(k) = phi1d(k) + alpha * Pk(k)
                Res(k) = Res(k) - alpha * Zk(k)
                R2 = R2 + Res(k) * Res(k) 
            enddo
        enddo
        s0 = sk

        ! Check for convergence using L2 Norm 
        Error = Sqrt(R2) / dble(NIJ)
        nIter = nIter + 1
        write(*,100) nIter, Error 
        write(16,110) nIter, Error
        100 Format('nIter = ', 1I6, 4X, 'Error = ', 1E15.8)
        110 Format(1I6, 2X, 1E15.8)

    enddo

    close(unit=16) ! close file after writen

end subroutine ICCG_Iter
