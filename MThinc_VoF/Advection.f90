!======================================================================================
! subroutine to advection flux in x direction
!--------------------------------------------------------------------------------------
Subroutine Advection_X
    use Mod_grid
    use Mod_flow
    implicit none
    integer :: i, j 

    do j = js+1 ,je-1 
        do i = js+1, ie-1
            phi(i,j) = phi(i,j) & 
                - ( flux_x(i,j) - flux_x(i-1,j) ) / hxi & 
                + phio(i,j)*( u(i,j) - u(i-1,j) )*delta_t / hxi 
        enddo
    enddo 
End Subroutine Advection_X
!======================================================================================
! subroutine to advection flux in y direction
!--------------------------------------------------------------------------------------
Subroutine Advection_Y
    use Mod_grid
    use Mod_flow
    implicit none
    integer :: i, j 

    do j = js+1 ,je-1 
        do i = js+1, ie-1
            phi(i,j) = phi(i,j) & 
                - ( flux_y(i,j) - flux_y(i,j-1) ) / hyi & 
                + phio(i,j)*( v(i,j) - v(i,j-1) )*delta_t / hyi 
        enddo
    enddo 
End Subroutine Advection_Y
!======================================================================================
