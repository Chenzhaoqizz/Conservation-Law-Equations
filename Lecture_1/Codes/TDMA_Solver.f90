!==================================================================================
! subroutine that uses TDMA algorothm to solve matrix
!----------------------------------------------------------------------------------
Subroutine TDMA(M,Dia,Lower,Up,B,x)
    implicit none 
    integer :: M, j
    real(8), dimension(M) :: Dia, Lower, Up, B
    real(8), dimension(M) :: x, P, Q

    ! Forward eliminate
    P(1) = - Up(1)/Dia(1);
    Q(1) = B(1)/Dia(1);

    do j = 2, M
        P(j) = -Up(j) / ( Dia(j) + Lower(j)*P(j-1) );
        Q(j) = ( B(j) - Lower(j)*Q(j-1) ) &
            / ( Dia(j) + Lower(j)*P(j-1) );
    enddo

    ! Backward substitute
    x(M) = Q(M);

    do j = M-1, 1, -1
        x(j) = P(j)*x(j+1) + Q(j);
    enddo

    return
End Subroutine TDMA
!----------------------------------------------------------------------------------
