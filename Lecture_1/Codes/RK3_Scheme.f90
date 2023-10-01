!==================================================================================
!    RK3 scheme for the unsteady advection diffusion equation
!----------------------------------------------------------------------------------
subroutine RK3_Scheme
    use Mod_UAD1d
    implicit none
    integer :: j
    real(8) :: lamb, beta
    real(8), dimension(Nx) :: f1

    fo = f
    f1 = 0.0
    lamb = 0.5d0*U*dt/hx
    beta = D*dt/(hx*hx)

    ! RK3 Step 1
    do j = is+1, ie-1
        f(j) = fo(j) & 
            + ( beta + lamb ) * fo(j-1) &
            - 2.d0 * beta * fo(j) &
            + ( beta - lamb ) * fo(j+1) 
    enddo

    ! RK3 Step 2
    do j = is+1, ie-1
        f1(j) = 0.75d0*fo(j) + 0.25d0*f(j) + 0.25d0 * ( &
            ( beta + lamb ) * f(j-1) - 2.d0*beta*f(j) + &
            ( beta - lamb ) * f(j+1) &
        )
    enddo

    ! RK3 Step 3
    do j = is+1, ie-1
        f(j) = (1.d0/3.d0)*fo(j) + (2.d0/3.d0)*f1(j) + (2.d0/3.d0) * ( &
            ( beta + lamb ) * f1(j-1) - 2.d0*beta*f1(j) + &
            ( beta - lamb ) * f1(j+1) &
        )
    enddo

    return 
end subroutine RK3_Scheme

