!==================================================================================
!   uses BiCGStab to solve the Poisson equation, 
!   suitable for non-systematic coefficient matrix 
!   preconditioned by Jacobi preconditioner
!----------------------------------------------------------------------------------
subroutine BiCGStab_Iter
    use Mod_Grid
    use Mod_Poisson2d
    implicit none

    integer :: i, j, k
    real(8) :: a, b, c, d, e, f, rho, rho_o 
    real(8) :: beta, alpha, omega
    real(8), dimension(NIJ) :: Res, Res0, Pk, Vk, Qk, Sk, Zk, Tk

    Pk = 0.d0
    Vk = 0.d0
    Qk = 0.d0
    Sk = 0.d0
    Zk = 0.d0
    Tk = 0.d0
    Res = 0.d0
    open(unit=16,file='Residual_BiCGStab.dat',status='unknown')

    ! Calculate initial residual vector
    do i = is+1, ie-1
        do j = js+1, je-1
            k = j + Ny*(i-1)
            a =  P(k) * phi1d(k)
            b = AN(k) * phi1d(k+1)
            c = AS(k) * phi1d(k-1)
            d = AE(k) * phi1d(k+Ny)
            e = AW(k) * phi1d(k-Ny)
            f = -4.d0 * sin( x(i) - y(j) ) * exp( x(i) - y(j) )
            Res(k) = f - a - b - c - d - e 
        enddo
    enddo

    Res0 = Res 
    rho = sum(Res*Res0)
    write(*,*) "rho = " , rho

    ! Begin iteration
    do while ( nIter <= MaxIter .and. Error >= MaxError )
    
        rho_o = rho
        rho = sum( Res*Res0 )

        ! Check for convergence using L2 Norm 
        Error = sqrt(rho) / dble(NIJ)

        if ( nIter .eq. 1 ) then 
            Pk = Res 
        else 
            beta = (alpha/omega)*(rho/rho_o)
            Pk = Res + beta*( Pk - omega*Vk )
        endif

        ! M.Q = pk, M(k) = -P(k) in Jacobi precondition 
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                Qk(k) = Pk(k) / (P(k))
            enddo
        enddo

        ! Calculate Vk = A.Qk
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                Vk(k) =  P(k) * Qk(k) &
                      + AE(k) * Qk(k+Ny) + AW(k) * Qk(k-Ny) &
                      + AN(k) * Qk(k+1)  + AS(k) * Qk(k-1)  
            enddo
        enddo

        alpha = sum( Res0 * Vk )
        alpha = rho / ( alpha + 1.0E-20 )
        Sk = Res - alpha * Vk

        ! M.Z = Sk, M(k) = -P(k) in Jacobi precondition 
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                Zk(k) = Sk(k) / (P(k))
            enddo
        enddo

        ! Calculate Tk = A.Zk
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                Tk(k) =  P(k) * Qk(k) &
                      + AE(k) * Qk(k+Ny) + AW(k) * Qk(k-Ny) &
                      + AN(k) * Qk(k+1)  + AS(k) * Qk(k-1)  
            enddo
        enddo

        omega = sum(Tk*Sk) / ( sum(Tk*Tk) + 1.0E-20 )

        ! Update the solution and residual vector
        do i = is+1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1)
                phi1d(k) = phi1d(k) + alpha * Qk(k) + omega * Zk(k)
                Res(k) = Sk(k) - omega * Tk(k)
            enddo
        enddo

        ! -----------
        nIter = nIter + 1
        write(*,100) nIter, Error 
        write(16,110) nIter, Error
        100 Format('nIter = ', 1I6, 4X, 'Error = ', 1E15.8)
        110 Format(1I6, 2X, 1E15.8)

    enddo

    close(unit=16) ! close file after writen

end subroutine BiCGStab_Iter


