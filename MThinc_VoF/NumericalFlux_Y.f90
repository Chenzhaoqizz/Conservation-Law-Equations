!======================================================================================
! subroutine to calculate numerical flux in x direction
!--------------------------------------------------------------------------------------
Subroutine NumericalFlux_Y
    use Mod_grid
    use Mod_flow
    implicit none
    integer :: i, j, jup
    real(8) :: alpha, x1, x2, y1, y2
    real(8) :: rp1, rp2, gs1, gs2, A1, A2, A3, A4
    real(8) :: A, B1, B2, Q, AA, BB, CC, DD1, DD2, D

    do i = is+1, ie-1
        do j = js+1, je-1 
            if ( v(i,j) >= 0.0 ) then 
                jup = j
                alpha = 1.0
                x1 = 0.0
                x2 = 1.0
                y1 = 1.0 - delta_t*v(i,j)/hyi
                y2 = 1.0 
            else if ( v(i,j) < 0.0 ) then 
                jup = j+1
                alpha = -1.0 
                x1 = 0.0
                x2 = 1.0
                y1 = 0.0
                y2 = -delta_t*v(i,j)/hyi
            endif

            if ( (phi(i,jup) < eplison) .or. (phi(i,jup) > 1.0 - eplison) ) then 
                ! 远离界面的网格面上的数值通量, simple upwind method is used. 
                flux_y(i,j) = delta_t * v(i,j) * phi(i,jup) 
            else if ( (phi(i,jup) >= eplison) .and. (phi(i,jup) <= 1.0 - eplison) ) then 
                ! 界面附近的网格面上的数值通量
                ! 主法线方向求解析解，副法线方向求高斯积分解
                if ( abs(NormVectorX(i,jup)) >= abs(NormVectorY(i,jup)) ) then 
                    ! X方向为主法线方向，Y方向为副法线方向
                    rp1 = ( 1.0 + 1.0 / sqrt(3.0) ) / 2.0 
                    rp2 = ( 1.0 - 1.0 / sqrt(3.0) ) / 2.0 
                    A   = exp( 2.0*beta*NormVectorX(i,jup) )
                    B1  = exp( 2.0*beta*( Lyy(i,jup)*(rp1**2)/2.0 + ( NormVectorY(i,jup) - Lyy(i,jup)/2.0 )*rp1 ) ) 
                    B2  = exp( 2.0*beta*( Lyy(i,jup)*(rp2**2)/2.0 + ( NormVectorY(i,jup) - Lyy(i,jup)/2.0 )*rp2 ) ) 
                    Q   = exp( 2.0*beta*NormVectorX(i,jup)*(2.0*phi(i,jup) - 1.0 ) ) 
                    AA  = A * B1 * B2 * ( A - Q ) 
                    BB  = A * ( B1 + B2 ) * ( 1.0 - Q ) 
                    CC  = 1.0 - A * Q 
                    DD1 = ( -BB + sqrt( BB**2 - 4.0*AA*CC ) )/( 2.0*AA ) 
                    DD2 = ( -BB - sqrt( BB**2 - 4.0*AA*CC ) )/( 2.0*AA ) 
                    D   = log( max(DD1, max(0.0, DD2)) )/(2.0*beta)
                    
                    if ( v(i,j) >= 0.0 ) then 
                        gs1 = (  1.0/sqrt(3.0)-1.0+2.0*hyi/( delta_t*v(i,j)+1.0e-16 ) )*( delta_t*v(i,j)/(2.0*hyi) ) 
                        gs2 = ( -1.0/sqrt(3.0)-1.0+2.0*hyi/( delta_t*v(i,j)+1.0e-16 ) )*( delta_t*v(i,j)/(2.0*hyi) ) 
                    else if ( v(i,j) < 0.0 ) then 
                        gs1 = (  1.0/sqrt(3.0) + 1.0 )*( -delta_t*v(i,j)/(2.0*hyi) ) 
                        gs2 = ( -1.0/sqrt(3.0) + 1.0 )*( -delta_t*v(i,j)/(2.0*hyi) ) 
                    endif

                    A1 = cosh( beta*( NormVectorX(i,jup)*x2 + Lyy(i,jup)*(gs1**2)/2.d0 + ( NormVectorY(i,jup) - Lyy(i,jup)/2.d0 )*gs1 + D ) ) 
                    A2 = cosh( beta*( NormVectorX(i,jup)*x2 + Lyy(i,jup)*(gs2**2)/2.d0 + ( NormVectorY(i,jup) - Lyy(i,jup)/2.d0 )*gs2 + D ) ) 
                    A3 = cosh( beta*( NormVectorX(i,jup)*x1 + Lyy(i,jup)*(gs1**2)/2.d0 + ( NormVectorY(i,jup) - Lyy(i,jup)/2.d0 )*gs1 + D ) ) 
                    A4 = cosh( beta*( NormVectorX(i,jup)*x1 + Lyy(i,jup)*(gs2**2)/2.d0 + ( NormVectorY(i,jup) - Lyy(i,jup)/2.d0 )*gs2 + D ) ) 
                    flux_y(i,j) = delta_t*v(i,j)/2.d0 + ( delta_t*v(i,j)/( 4.d0*beta*NormVectorX(i,jup) ) )*Log( (A1*A2)/(A3*A4) ) 
                else if ( abs(NormVectorX(i,jup)) < abs(NormVectorY(i,jup)) ) then 
                    ! Y方向为主法线方向，X方向为副法线方向
                    rp1 = ( 1.0 + 1.0 / sqrt(3.0) ) / 2.0 
                    rp2 = ( 1.0 - 1.0 / sqrt(3.0) ) / 2.0 
                    A   = exp( 2.0*beta*NormVectorY(i,jup) ) 
                    B1  = exp( 2.0*beta*( Lxx(i,jup)*(rp1**2)/2.0 + ( NormVectorX(i,jup) - Lxx(i,jup)/2.0 )*rp1 ) ) 
                    B2  = exp( 2.0*beta*( Lxx(i,jup)*(rp2**2)/2.0 + ( NormVectorX(i,jup) - Lxx(i,jup)/2.0 )*rp2 ) ) 
                    Q   = exp( 2.0*beta*NormVectorY(i,jup)*(2.0*phi(i,jup) - 1.0 ) ) 
                    AA  = A * B1 * B2 * ( A - Q ) 
                    BB  = A * ( B1 + B2 ) * ( 1.0 - Q ) 
                    CC  = 1.0 - A * Q 
                    DD1 = ( -BB + sqrt( BB**2 - 4.0*AA*CC ) )/( 2.0*AA ) 
                    DD2 = ( -BB - sqrt( BB**2 - 4.0*AA*CC ) )/( 2.0*AA ) 
                    D   = log( max(DD1, max(0.0, DD2)) )/(2.0*beta) 
                    
                    gs1 = ( 1.0 + 1.0 / sqrt(3.0) ) / 2.0 
                    gs2 = ( 1.0 - 1.0 / sqrt(3.0) ) / 2.0 
                    
                    A1 = cosh( beta*( Lxx(i,jup)*(gs1**2)/2.d0 + ( NormVectorX(i,jup) - Lxx(i,jup)/2.d0 ) * gs1 + NormVectorY(i,jup)*y2 + D ) ) 
                    A2 = cosh( beta*( Lxx(i,jup)*(gs2**2)/2.d0 + ( NormVectorX(i,jup) - Lxx(i,jup)/2.d0 ) * gs2 + NormVectorY(i,jup)*y2 + D ) ) 
                    A3 = cosh( beta*( Lxx(i,jup)*(gs1**2)/2.d0 + ( NormVectorX(i,jup) - Lxx(i,jup)/2.d0 ) * gs1 + NormVectorY(i,jup)*y1 + D ) ) 
                    A4 = cosh( beta*( Lxx(i,jup)*(gs2**2)/2.d0 + ( NormVectorX(i,jup) - Lxx(i,jup)/2.d0 ) * gs2 + NormVectorY(i,jup)*y1 + D ) ) 
                    flux_y(i,j) = delta_t*v(i,j)/2.d0 + ( alpha*hyi/( 4.d0*beta*NormVectorY(i,jup) ) )*Log( (A1*A2)/(A3*A4) ) 
                endif
            endif 

        enddo
    enddo

End Subroutine NumericalFlux_Y
!======================================================================================
