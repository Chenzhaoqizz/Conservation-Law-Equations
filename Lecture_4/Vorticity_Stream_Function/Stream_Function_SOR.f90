!======================================================================================
!      solve the possion equation of stream function by SOR method
!--------------------------------------------------------------------------------------
subroutine Stream_Function_SOR
    use flow_var
    implicit none
    integer :: i, j, nIter
    real(8) :: beta, R2, Residual

    nIter = 0
    Residual = 100.d0
    R2 = 0.d0
    beta = 1.75d0

    write(*,*) "Solve possion equation by SOR "

    ! Begin iteration 
    do while ( nIter <= MaxIter .and. Residual >= Epsilon )

        ! Iteration for inner grid point
        do i = is+1, ie-1
            do j = js+1, je-1
                sf(i,j) = beta * ( - vort(i,j) & 
                    & - AE(i,j) * sf(i+1,j) - AW(i,j) * sf(i-1,j) &
                    & - AN(i,j) * sf(i,j+1) - AS(i,j) * sf(i,j-1) & 
                    & ) / AP(i,j) + ( 1.d0 - beta ) * sf(i,j) 
            enddo
        enddo

        Residual = 0.d0
        ! Calculate residual using L2 Norm of inner grid point
        do i = is+1, ie-1
            do j = js+1, je-1
                R2 =  - vort(i,j) - AE(i,j) * sf(i+1,j) - AW(i,j) * sf(i-1,j) &
                    & - AN(i,j) * sf(i,j+1) - AS(i,j) * sf(i,j-1) - AP(i,j) * sf(i,j) 
                Residual = Residual + R2**2
            enddo
        enddo

        Residual = sqrt( Residual ) /  dble( (ie-is) * (je-js) )
        nIter = nIter + 1
    enddo

    write(*,20) nIter, Residual 
    20 Format('nIter = ', 1I6, ', Residual = ', 1ES15.8)

end subroutine Stream_Function_SOR
