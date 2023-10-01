!==================================================================================
! subroutine that uses method of steepest descent (MSD) 
! to solve the Poisson equation, only suitable for 
! systematic positive define matrix 
!----------------------------------------------------------------------------------
Subroutine MSD_Iter
    use Mod_Grid
    use Mod_Poisson2d
    implicit none 
    integer :: i, j , k
    real(8) :: f, a, b, c, d, e, R2, alpha, RAR
    real(8), dimension(NIJ) :: Res, ARes

    open(unit=16,file='Residual_MSD.dat',status='unknown')

    ! Begin iteration 
    do while ( nIter <= MaxIter .and. Error >= MaxError )

        R2 = 0.d0
        Res = 0.d0
        RAR = 0.d0 
        alpha = 0.d0
    
        ! Calculate the residual vectors 
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
    
        ! Calculate the space step: alpha 
        do i = is +1, ie-1
            do j = js+1, je-1
                k = j + Ny*(i-1) 
                ARes(k) =  P(k)*Res(k) & 
                    &   + AN(k)*Res(k+1)  + AS(k)*Res(k-1) & 
                    &   + AE(k)*Res(k+Ny) + AW(k)*Res(k-Ny) 
                
                RAR = RAR + Res(k) * ARes(k) 
            enddo
        enddo 
    
        alpha = R2 / (RAR + 1E-20) 
    
        ! Update the solution 
        phi1d(:) = phi1d(:) + alpha * Res(:) 
        ! do i = is +1, ie-1
        !     do j = js+1, je-1
        !         k = j + Ny*(i-1) 
        !         phi1d(k) = phi1d(k) + alpha * Res(k) 
        !     enddo
        ! enddo 

        ! Check for convergence using L2 Norm 
        Error = Sqrt( R2 ) / dble( NIJ ) 
        nIter = nIter + 1
        write(*,100) nIter, Error 
        write(16,110) nIter, Error
        100 Format('nIter = ', 1I6, 4X, 'Error = ', 1E15.8)
        110 Format(1I6, 2X, 1E15.8)
    enddo 

    close(unit=16) ! close file after writen

End Subroutine MSD_Iter
