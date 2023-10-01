!======================================================================================
!                updata boundary values of vorticity
!--------------------------------------------------------------------------------------
subroutine Vorticity_BC
    use flow_var
    implicit none
    real(8) :: hy1, hy2, hx1, hx2

    hy1 = y(js+1) - y(js)
    hy2 = y(je) - y(je-1)
    hx1 = x(is+1) - x(is)
    hx2 = x(ie) - x(ie-1)

    ! Bottom boundary
    vort(:,js) = - 2.d0 * sf(:,js+1) / hy1**2

    ! Top boundary
    vort(:,je) = - 2.d0 * sf(:,je-1) / hy2**2 - 2.d0 * uTop / hy2

    ! Left boundary
    vort(is,:) = - 2.d0 * sf(is+1,:) / hx1**2

    ! Right boundary
    vort(ie,:) = - 2.d0 * sf(ie-1,:) / hx2**2

end subroutine Vorticity_BC
