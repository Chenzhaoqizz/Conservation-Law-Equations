!======================================================================================
! subroutine to calculate normal vector and curvature by volume fraction
!--------------------------------------------------------------------------------------
subroutine Normal_Vectors
    use Mod_grid
    use Mod_flow
    implicit none
    integer :: i,j
    real(8), dimension(Nx, Ny) :: MCorX, MCorY, mx, my, NxCor, NyCor

    MCorX = 0.0; MCorY = 0.0; mx = 0.0; my = 0.0;
    NxCor = 0.0; NyCor = 0.0;

    ! 计算体积分数φ的梯度
    do i = is+1, ie-1
        do j = js+1, je-1
            MCorX(i,j) = ( phi(i+1,j) + phi(i+1,j+1) - phi(i,j) - phi(i,j+1) ) / ( hxi + hxi ) 
            MCorY(i,j) = ( phi(i,j+1) + phi(i+1,j+1) - phi(i,j) - phi(i+1,j) ) / ( hyi + hyi ) 
        enddo
    enddo

    ! 将周围几个网格的体积分数φ的梯度在界面附近进行简单代数平均 
    do i = is+1, ie-1
        do j = js+1, je-1
            if ( (phi(i,j) >= eplison) .and. (phi(i,j) <= 1.0 - eplison) ) then 
                mx(i,j) = 0.25*( MCorX(i,j) + MCorX(i+1,j) + MCorX(i,j+1) + MCorX(i+1,j+1) ) 
                my(i,j) = 0.25*( MCorY(i,j) + MCorY(i+1,j) + MCorY(i,j+1) + MCorY(i+1,j+1) ) 
            endif
        enddo
    enddo

    ! 基于体积分数φ的梯度,计算单位法向量,只在界面附近进行
    do i = is+1, ie-1
        do j = js+1, je-1
            if ( (phi(i,j) > eplison) .and. (phi(i,j) < 1.0 - eplison) ) then 
                NormVectorX(i,j) = hxi*mx(i,j) / dSqrt( mx(i,j)**2 + my(i,j)**2 + 1.0e-16 ) 
                NormVectorY(i,j) = hyi*my(i,j) / dSqrt( mx(i,j)**2 + my(i,j)**2 + 1.0e-16 ) 
                NxCor(i,j) = MCorX(i,j) / dSqrt( MCorX(i,j)**2 + MCorY(i,j)**2 + 1.0e-16 ) 
                NyCor(i,j) = MCorY(i,j) / dSqrt( MCorX(i,j)**2 + MCorY(i,j)**2 + 1.0e-16 ) 
            endif
        enddo
    enddo

    ! 基于单位法向量,计算空间曲率
    do i = is+1, ie-1
        do j = js+1, je-1
            if ( (phi(i,j) > eplison) .and. (phi(i,j) < 1.0 - eplison) ) then 
                Lxx(i,j) = 0.50 * hxi * ( NxCor(i,j-1) + NxCor(i,j) - NxCor(i-1,j-1) - NxCor(i-1,j) ) 
                Lyy(i,j) = 0.50 * hyi * ( NyCor(i-1,j) + NyCor(i,j) - NyCor(i-1,j-1) - NyCor(i,j-1) ) 
                Lxy(i,j) = 0.25 * hxi * ( NxCor(i,j-1) + NxCor(i,j) - NxCor(i-1,j-1) - NxCor(i-1,j) ) &
                         + 0.25 * hyi * ( NyCor(i-1,j) + NyCor(i,j) - NyCor(i-1,j-1) - NyCor(i,j-1) ) 
            endif
        enddo
    enddo

end subroutine Normal_Vectors 
