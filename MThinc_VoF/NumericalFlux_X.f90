!======================================================================================
! subroutine to calculate numerical flux in x direction
!--------------------------------------------------------------------------------------
Subroutine NumericalFlux_X 
    use Mod_grid
    use Mod_flow
    implicit none
    integer :: i, j, iup
    real(8) :: alpha, x1, x2, y1, y2
    real(8) :: rp1, rp2, gs1, gs2, A1, A2, A3, A4
    real(8) :: A, B1, B2, Q, AA, BB, CC, DD1, DD2, D

    do j = js+1, je-1  
        do i = is+1, ie-1 
            if ( u(i,j) >= 0.0 ) then 
                iup = i
                alpha = 1.0
                ! 局部坐标系下的积分区间 [x1, x2], [y1, y2] 
                x1 = 1.0 - delta_t * u(i,j) / hxi
                x2 = 1.0 
                y1 = 0.0
                y2 = 1.0 
            else if ( u(i,j) < 0.0 ) then 
                iup = i+1
                alpha = -1.0
                x1 = 0.0
                x2 = - delta_t * u(i,j) / hxi
                y1 = 0.0
                y2 = 1.0
            endif

            if ( (phi(iup,j) < eplison) .or. (phi(iup,j) > 1.0 - eplison) ) then 
                ! 远离界面的网格面上的数值通量, simple upwind method is used. 
                flux_x(i,j) = delta_t * u(i,j) * phi(iup,j) 
            else if ( (phi(iup,j) >= eplison) .and. (phi(iup,j) <= 1.0 - eplison) ) then 
                ! 界面附近的网格面上的数值通量
                ! 主法线方向求解析解，副法线方向求高斯积分解
                if ( abs(NormVectorX(iup,j)) >= abs(NormVectorY(iup,j)) ) then 
                    ! X方向为主法线方向，Y方向为副法线方向
                    rp1 = ( 1.0 + 1.0 / sqrt(3.0) ) / 2.0 
                    rp2 = ( 1.0 - 1.0 / sqrt(3.0) ) / 2.0 
                    A   = exp( 2.0*beta*NormVectorX(iup,j) )
                    B1  = exp( 2.0*beta*( Lyy(iup,j)*(rp1**2)/2.0 + ( NormVectorY(iup,j) - Lyy(iup,j)/2.0 )*rp1 ) ) 
                    B2  = exp( 2.0*beta*( Lyy(iup,j)*(rp2**2)/2.0 + ( NormVectorY(iup,j) - Lyy(iup,j)/2.0 )*rp2 ) ) 
                    Q   = exp( 2.0*beta*NormVectorX(iup,j)*(2.0*phi(iup,j) - 1.0 ) ) 
                    AA  = A * B1 * B2 * ( A - Q ) 
                    BB  = A * ( B1 + B2 ) * ( 1.0 - Q ) 
                    CC  = 1.0 - A * Q 
                    DD1 = ( -BB + sqrt( BB**2 - 4.0*AA*CC ) )/( 2.0*AA ) 
                    DD2 = ( -BB - sqrt( BB**2 - 4.0*AA*CC ) )/( 2.0*AA ) 
                    D   = log( max(DD1, max(0.0, DD2)) )/(2.0*beta)
                    
                    gs1 = ( 1.0 + 1.0 / sqrt(3.0) ) / 2.0 
                    gs2 = ( 1.0 - 1.0 / sqrt(3.0) ) / 2.0 
                    A1 = cosh( beta*( NormVectorX(iup,j)*x2 + Lyy(iup,j)*(gs1**2)/2.d0 + ( NormVectorY(iup,j) - Lyy(iup,j)/2.d0 )*gs1 + D ) ) 
                    A2 = cosh( beta*( NormVectorX(iup,j)*x2 + Lyy(iup,j)*(gs2**2)/2.d0 + ( NormVectorY(iup,j) - Lyy(iup,j)/2.d0 )*gs2 + D ) ) 
                    A3 = cosh( beta*( NormVectorX(iup,j)*x1 + Lyy(iup,j)*(gs1**2)/2.d0 + ( NormVectorY(iup,j) - Lyy(iup,j)/2.d0 )*gs1 + D ) ) 
                    A4 = cosh( beta*( NormVectorX(iup,j)*x1 + Lyy(iup,j)*(gs2**2)/2.d0 + ( NormVectorY(iup,j) - Lyy(iup,j)/2.d0 )*gs2 + D ) ) 
                    flux_x(i,j) = delta_t*u(i,j)/2.d0 + ( alpha*hxi/( 4.d0*beta*NormVectorX(iup,j) ) )*Log( (A1*A2)/(A3*A4) ) 
                else if ( abs(NormVectorX(iup,j)) < abs(NormVectorY(iup,j)) ) then 
                    ! Y方向为主法线方向，X方向为副法线方向
                    rp1 = ( 1.0 + 1.0 / sqrt(3.0) ) / 2.0 
                    rp2 = ( 1.0 - 1.0 / sqrt(3.0) ) / 2.0 
                    A   = exp( 2.0*beta*NormVectorY(iup,j) ) 
                    B1  = exp( 2.0*beta*( Lxx(iup,j)*(rp1**2)/2.0 + ( NormVectorX(iup,j) - Lxx(iup,j)/2.0 )*rp1 ) ) 
                    B2  = exp( 2.0*beta*( Lxx(iup,j)*(rp2**2)/2.0 + ( NormVectorX(iup,j) - Lxx(iup,j)/2.0 )*rp2 ) ) 
                    Q   = exp( 2.0*beta*NormVectorY(iup,j)*(2.0*phi(iup,j) - 1.0 ) ) 
                    AA  = A * B1 * B2 * ( A - Q ) 
                    BB  = A * ( B1 + B2 ) * ( 1.0 - Q ) 
                    CC  = 1.0 - A * Q 
                    DD1 = ( -BB + sqrt( BB**2 - 4.0*AA*CC ) )/( 2.0*AA ) 
                    DD2 = ( -BB - sqrt( BB**2 - 4.0*AA*CC ) )/( 2.0*AA ) 
                    D   = log( max(DD1, max(0.0, DD2)) )/(2.0*beta) 
                    
                    if ( u(i,j) >= 0.0 ) then 
                        gs1 = (  1.0/sqrt(3.0)-1.0+2.0*hxi/( delta_t*u(i,j)+1.0e-16 ) )*( delta_t*u(i,j)/(2.0*hxi) ) 
                        gs2 = ( -1.0/sqrt(3.0)-1.0+2.0*hxi/( delta_t*u(i,j)+1.0e-16 ) )*( delta_t*u(i,j)/(2.0*hxi) ) 
                    else if ( u(i,j) < 0.0 ) then 
                        gs1 = (  1.0/sqrt(3.0) + 1.0 )*( -delta_t*u(i,j)/(2.0*hxi) ) 
                        gs2 = ( -1.0/sqrt(3.0) + 1.0 )*( -delta_t*u(i,j)/(2.0*hxi) ) 
                    endif
                    
                    A1 = cosh( beta*( Lxx(iup,j)*(gs1**2)/2.d0 + ( NormVectorX(iup,j) - Lxx(iup,j)/2.d0 ) * gs1 + NormVectorY(iup,j)*y2 + D ) ) 
                    A2 = cosh( beta*( Lxx(iup,j)*(gs2**2)/2.d0 + ( NormVectorX(iup,j) - Lxx(iup,j)/2.d0 ) * gs2 + NormVectorY(iup,j)*y2 + D ) ) 
                    A3 = cosh( beta*( Lxx(iup,j)*(gs1**2)/2.d0 + ( NormVectorX(iup,j) - Lxx(iup,j)/2.d0 ) * gs1 + NormVectorY(iup,j)*y1 + D ) ) 
                    A4 = cosh( beta*( Lxx(iup,j)*(gs2**2)/2.d0 + ( NormVectorX(iup,j) - Lxx(iup,j)/2.d0 ) * gs2 + NormVectorY(iup,j)*y1 + D ) ) 
                    flux_x(i,j) = delta_t*u(i,j)/2.d0 + ( delta_t*u(i,j)/( 4.d0*beta*NormVectorY(iup,j) ) )*Log( (A1*A2)/(A3*A4) ) 
                endif
            endif
        enddo
    enddo

End Subroutine NumericalFlux_X
!======================================================================================
