!==================================================================================
! subroutine that uses CNCS scheme to solve 
! the unsteady advection diffusion equation
!----------------------------------------------------------------------------------
Subroutine CNCS_Scheme
    use Mod_UAD1d
    implicit none
    integer :: M, j
    real(8) :: beta, lamb
    real(8), dimension(Nx-2) :: Dia, Lower, Up, RM, x

    fo=f;
    M = Nx-2;

    lamb = 0.5d0*U*dt/hx;
    beta = D*dt/(hx*hx);

    ! Assemble the coefficient matrix
    Dia(1:M) = 2.d0 + 2.d0*beta;
    
    Lower(1) = 0.d0;
    Lower(2:M) = -(beta+lamb);

    Up(1:M-1) = lamb - beta;
    Up(M) = 0.d0;

    ! Assemble the Right hand side(RHS.)
    do j = 2, M-1
        RM(j) = (beta+lamb)*f(j) &
            + 2.d0*(1.d0-beta)*f(j+1) &
            + (beta-lamb)*f(j+2)
    enddo

    RM(1) = 2.d0*(beta+lamb)*f(1) &
        + 2.d0*(1.d0-beta)*f(2) &
        + (beta-lamb)*f(3)
    
    RM(M) = (beta+lamb)*f(M) &
        + 2.d0*(1.d0-beta)*f(M+1) &
        + 2.d0*(beta-lamb)*f(M+2)

    call TDMA(M,Dia,Lower,Up,RM,x)
    f(2:Nx-1) = x(1:M);
    
    call SetupBC

    return
End Subroutine CNCS_Scheme
!----------------------------------------------------------------------------------
