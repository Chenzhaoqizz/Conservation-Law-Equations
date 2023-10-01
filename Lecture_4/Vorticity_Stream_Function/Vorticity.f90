!======================================================================================
!                solve vorticity equation by RK-3
!--------------------------------------------------------------------------------------
subroutine Vorticity
    use flow_var
    implicit none
    integer :: i, j
    real(8) :: hx, hxL, hxR, hy, hyL, hyR
    real(8) :: advection, diffusion
    real(8) :: advx, advy, visx, visy
    real(8), dimension(Nx, Ny) :: vorto, vort1

    vorto = vort

    ! RK-3, Step 1:
    vort1 = vort
    do i = is+1, ie-1
        hxL = x(i) - x(i-1)
        hxR = x(i+1) - x(i)
        hx = hxR + hxL
        do j = js+1, je-1
            hyL = y(j) - y(j-1)
            hyR = y(j+1) - y(j)
            hy = hyR + hyL

            advx = advection( sf(i,j+1), sf(i,j-1), hy, vort(i+1,j), vort(i-1,j), hx )
            advy = advection( sf(i+1,j), sf(i-1,j), hx, vort(i,j+1), vort(i,j-1), hy )
            visx = diffusion( vort1(i-1,j), vort1(i,j), vort1(i+1,j), hxL, hxR, hx )
            visy = diffusion( vort1(i,j-1), vort1(i,j), vort1(i,j+1), hyL, hyR, hy )

            vort(i,j) = vorto(i,j) + delta_t * ( - advx + advy + ( visx + visy ) / Re )
        enddo
    enddo

    ! RK-3, Step 2:
    vort1 = vort
    do i = is+1, ie-1
        hxL = x(i) - x(i-1)
        hxR = x(i+1) - x(i)
        hx = hxR + hxL
        do j = js+1, je-1
            hyL = y(j) - y(j-1)
            hyR = y(j+1) - y(j)
            hy = hyR + hyL

            advx = advection( sf(i,j+1), sf(i,j-1), hy, vort(i+1,j), vort(i-1,j), hx )
            advy = advection( sf(i+1,j), sf(i-1,j), hx, vort(i,j+1), vort(i,j-1), hy )
            visx = diffusion( vort1(i-1,j), vort1(i,j), vort1(i+1,j), hxL, hxR, hx )
            visy = diffusion( vort1(i,j-1), vort1(i,j), vort1(i,j+1), hyL, hyR, hy )

            vort(i,j) = 0.75d0 * vorto(i,j) + 0.25d0 * vort1(i,j) & 
                & + 0.25d0 * delta_t * ( - advx + advy + ( visx + visy ) / Re )
        enddo
    enddo

    ! RK-3, Step 3:
    vort1 = vort
    do i = is+1, ie-1
        hxL = x(i) - x(i-1)
        hxR = x(i+1) - x(i)
        hx = hxR + hxL
        do j = js+1, je-1
            hyL = y(j) - y(j-1)
            hyR = y(j+1) - y(j)
            hy = hyR + hyL

            advx = advection( sf(i,j+1), sf(i,j-1), hy, vort(i+1,j), vort(i-1,j), hx )
            advy = advection( sf(i+1,j), sf(i-1,j), hx, vort(i,j+1), vort(i,j-1), hy )
            visx = diffusion( vort1(i-1,j), vort1(i,j), vort1(i+1,j), hxL, hxR, hx )
            visy = diffusion( vort1(i,j-1), vort1(i,j), vort1(i,j+1), hyL, hyR, hy )

            vort(i,j) = (1.d0/3.d0) * vorto(i,j) + (2.d0/3.d0) * vort1(i,j) & 
                & + (2.d0/3.d0) * delta_t * ( - advx + advy + ( visx + visy ) / Re )
        enddo
    enddo
    
end subroutine Vorticity
!======================================================================================
!--------------------------------------------------------------------------------------
real(8) function advection( strfL, strfR, dsf, vortL, vortR, dvort )
    implicit none
    real(8), intent(in) :: dsf, dvort
    real(8), intent(in) :: strfL, strfR, vortL, vortR

    advection = ( ( strfR - strfL ) / dsf ) * ( ( vortR - vortL ) / dvort )
end function advection
!======================================================================================
!--------------------------------------------------------------------------------------
real(8) function diffusion( vtL, vt, vtR, dL, dR, d )
    implicit none
    real(8), intent(in) :: vtL, vt, vtR, dL, dR, d

    diffusion = ( ( vtR - vt ) / dR - ( vtL - vt ) / dL ) / ( 0.5d0 * d )
end function diffusion
