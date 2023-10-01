!==================================================================================
! subroutine that uses FTCS scheme to solve 
! the unsteady advection diffusion equation
!----------------------------------------------------------------------------------
Subroutine FTCS_Scheme
    use Mod_UAD1d
    implicit none 
    integer :: j
    real(8) :: lamb, beta

    fo = f; 
    lamb = 0.5d0*U*dt/hx;
    beta = D*dt/(hx*hx);

    do j = is+1, ie-1;
        f(j) = (beta+lamb)*fo(j-1) &
            + (1.d0-2.d0*beta)*fo(j) &
            + (beta-lamb)*fo(j+1) 
    enddo

    call SetupBC

    return
End Subroutine FTCS_Scheme
