!==================================================================================
! subroutine that uses conjugate gradient squared method (CGStab) 
! to solve the Poisson equation, suitable for non-systematic  
! coefficient matrix 
!----------------------------------------------------------------------------------
subroutine CGStab_Iter
    use Mod_Grid
    use Mod_Poisson2d
    implicit none
    integer :: i,j,k 
    real(8) :: f, a, b, c, d, e, R2, R20, alpha, beta, gamma 
    real(8) :: s0, sk, beta0, omega
    real(8), dimension(NIJ) :: Res, Res0, Diag, Zk, Pk, Uk, Vk

    R2 = 0.d0
    R20 = 0.d0
    Res = 0.d0 
    Res0 = 0.d0
    Diag = 0.d0
    Zk = 0.d0
    Pk = 0.d0
    Uk = 0.d0
    Vk = 0.d0
    beta = 1.d0 
    beta0 = 1.d0
    alpha = 1.d0
    gamma = 1.d0 
    omega = 0.d0
    s0 = 0.d0
    sk = 0.d0

    open(unit=16,file='Residual_CGStab.dat',status='unknown')

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
    Res0 = Res

    ! Calculate elements of preconditioning matrix diagonal
    do i = is + 1, ie -1 
        do j = js + 1, je -1
            k = j + Ny*(i-1)
            Diag(k) = 1.d0 / ( P(k) &
              - AW(k) * Diag(k-Ny) * AE(k-Ny) &
              - AS(k) * Diag(k-1 ) * AN(k-1)  &
            )
        enddo
    enddo

    ! Begin iteration 
    do while ( nIter <= MaxIter .and. Error >= MaxError )

        ! Calculate beta and omega
        beta = 0.d0
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                beta = beta + Res(k) * Res0(k) 
            enddo
        enddo
        omega = beta*gamma / (alpha*beta0 + 1.0e-30)
        beta0 = beta

        ! Calculate Pk
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                Pk(k) = Res(k) + omega * ( Pk(k) - alpha * Uk(k) )
            enddo
        enddo

        ! Solve (M Zk = Pk) - forward substitution
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                Zk(k) = (   Pk(k) &
                          - AW(k) * Zk(k-Ny) &
                          - AS(k) * Zk(k-1)  &
                ) * Diag(k) 
            enddo
        enddo

        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                Zk(k) = Zk(k) / ( Diag(k) + 1.0e-30 )
            enddo
        enddo

        ! Backward substitution
        do i = ie-1, is+1, -1
            do j = je-1, js+1, -1
                k = j + Ny*(i-1)
                Zk(k) = (   Zk(k) &
                          - AE(k) * Zk(k+Ny) &
                          - AN(k) * Zk(k+1)  &
                ) * Diag(k) 
            enddo
        enddo

        ! Calculate Uk = A.Zk
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                Uk(k) =  P(k) * Zk(k) &
                      + AE(k) * Zk(k+Ny) + AW(k) * Zk(k-Ny) &
                      + AN(k) * Zk(k+1)  + AS(k) * Zk(k-1)  
            enddo
        enddo

        ! Calculate scalar product Uk.Reso and gamma
        s0 = 0.0
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                s0 = s0 + Uk(k) * Res0(k)
            enddo
        enddo
        gamma = beta / ( s0 + 1.0e-30 )

        ! Update the solution and calculate w (overwtite res - it is res-update)
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                phi1d(k) = phi1d(k) + gamma * Zk(k)
                Res(k) = Res(k) - gamma * Uk(k)
            enddo
        enddo

        ! Solve (M Y = W); Y overwrites Zk; forward substitution
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                Zk(k) = (   Res(k) &
                          - AW(k) * Zk(k-Ny) &
                          - AS(k) * Zk(k-1)  &
                ) * Diag(k) 
            enddo
        enddo

        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                Zk(k) = Zk(k) / (Diag(k) + 1.0e-30) 
            enddo
        enddo

        ! Backward substitution
        do i = ie-1, is+1, -1
            do j = je-1, js+1, -1
                k = j + Ny*(i-1)
                Zk(k) = (   Zk(k) &
                          - AE(k) * Zk(k+Ny) &
                          - AN(k) * Zk(k+1)  &
                ) * Diag(k) 
            enddo
        enddo

        ! Calculate V = A.Y (Vk = A.Zk)
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                Vk(k) =  P(k) * Zk(k) &
                      + AE(k) * Zk(k+Ny) + AW(k) * Zk(k-Ny) &
                      + AN(k) * Zk(k+1)  + AS(k) * Zk(k-1)  
            enddo
        enddo

        ! Calculate alpha
        s0 = 0.d0
        sk = 0.d0
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                s0 = s0 + Vk(k) * Res(k)
                sk = sk + Vk(k) * Vk(k)  
            enddo
        enddo

        alpha = s0 / ( sk + 1.0e-30 )

        ! Update the solution and residual vector
        R2 = 0.d0
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                phi1d(k) = phi1d(k) + alpha * Zk(k)
                Res(k) = Res(k) - alpha * Vk(k)
                R2 = R2 + Res(k) * Res(k) 
            enddo
        enddo

        ! Check for convergence using L2 Norm 
        Error = Sqrt(R2) / dble(NIJ)
        nIter = nIter + 1
        write(*,100) nIter, Error 
        write(16,110) nIter, Error
        100 Format('nIter = ', 1I6, 4X, 'Error = ', 1E15.8)
        110 Format(1I6, 2X, 1E15.8)

    enddo

    close(unit=16) ! close file after writen

end subroutine CGStab_Iter
