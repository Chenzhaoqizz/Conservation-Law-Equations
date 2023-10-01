!***********************************************************************
!    solves div (1/r grad ps) = S on a regular staggered grid
!    default boundary conditions for ps: Neumann 
!-----------------------------------------------------------------------
subroutine ppe_bicgstab2D(hxi,hyi,nx,ny,r,ps,S,resmax,maxit)
    implicit double precision (a-h,o-z)
    dimension ps(nx+2,ny+2),S(nx+2,ny+2),r(nx+2,ny+2)
    dimension res(nx+2,ny+2)
    dimension reso(nx+2,ny+2)
    dimension pk(nx+2,ny+2),tk(nx+2,ny+2)
    dimension zk(nx+2,ny+2),vk(nx+2,ny+2)
    dimension sk(nx+2,ny+2),qk(nx+2,ny+2)
    dimension aR(nx+2,ny+2),aL(nx+2,ny+2),aT(nx+2,ny+2),aB(nx+2,ny+2)
    dimension aC(nx+2,ny+2)

    res=0.0d0
    reso=0.0d0
    pk=0.0d0
    tk=0.0d0
    zk=0.0d0
    vk=0.0d0
    sk=0.0d0
    qk=0.0d0
    two=2.0d0

    nxp1=nx+1
    nyp1=ny+1
    nxp2=nx+2
    nyp2=ny+2

    aR(2:nxp1,2:nyp1) = two * hxi * hxi &
        &    / ( r(3:nxp2,2:nyp1) + r(2:nxp1,2:nyp1) )
    aL(2:nxp1,2:nyp1) = two * hxi * hxi &
        &    / ( r(2:nxp1,2:nyp1) + r(1:nx,2:nyp1)   )
    aT(2:nxp1,2:nyp1) = two * hyi * hyi &
        &    / ( r(2:nxp1,3:nyp2) + r(2:nxp1,2:nyp1) )
    aB(2:nxp1,2:nyp1) = two * hyi * hyi & 
        &    / ( r(2:nxp1,2:nyp1) + r(2:nxp1,1:ny)   )
    aC=aR+aL+aT+aB

    res(2:nxp1,2:nyp1) = S(2:nxp1,2:nyp1) & 
        & - ( aR(2:nxp1,2:nyp1) * ps(3:nxp2,2:nyp1) & 
        & +   aL(2:nxp1,2:nyp1) * ps(1:nx,2:nyp1)   & 
        & +   aT(2:nxp1,2:nyp1) * ps(2:nxp1,3:nyp2) & 
        & +   aB(2:nxp1,2:nyp1) * ps(2:nxp1,1:ny)   & 
        & -   aC(2:nxp1,2:nyp1) * ps(2:nxp1,2:nyp1) ) 

    reso = res 

    do l=1,maxit

        rho_o = rho;
        rho = sum(res*reso);

        if(l.eq.1)then
            pk=res
        else
            beta=(rho/rho_o)*(alpha/omega)
            pk=res+beta*(pk-omega*vk)
        endif
        
        ! M P = Pi, Jacobi precondition, here, M = aC 
        zk=pk/aC

        ! Neumann BC's
        zk(1:nxp2,1)=zk(1:nxp2,2)
        zk(1:nxp2,nyp2)=zk(1:nx+2,nyp1)
        zk(1,1:nyp2)=zk(2,1:nyp2)
        zk(nxp2,1:nyp2)=zk(nxp1,1:nyp2) 
        
        vk(2:nxp1,2:nyp1) = aR(2:nxp1,2:nyp1)*zk(3:nxp2,2:nyp1) & 
            & + aL(2:nxp1,2:nyp1)*zk(1:nx,2:nyp1)   & 
            & + aT(2:nxp1,2:nyp1)*zk(2:nxp1,3:nyp2) & 
            & + aB(2:nxp1,2:nyp1)*zk(2:nxp1,1:ny)   & 
            & - aC(2:nxp1,2:nyp1)*zk(2:nxp1,2:nyp1) 
        
        rtvk = sum( reso * vk )
        alpha = rho / ( rtvk + 1.e-30 )
        sk = res - alpha * vk

        ! M s = si, Jacobi precondition, here, M = aC
        qk = sk / aC

        ! Neumann BC's
        qk(1:nxp2,1)=qk(1:nxp2,2)
        qk(1:nxp2,nyp2)=qk(1:nxp2,nyp1)
        qk(1,1:nyp2)=qk(2,1:nyp2)
        qk(nxp2,1:nyp2)=qk(nxp1,1:nyp2)

        tk(2:nxp1,2:nyp1) = aR(2:nxp1,2:nyp1)*qk(3:nxp2,2:nyp1) &
            & + aL(2:nxp1,2:nyp1)*qk(1:nx,2:nyp1)    & 
            & + aT(2:nxp1,2:nyp1)*qk(2:nxp1,3:nyp2)  & 
            & + aB(2:nxp1,2:nyp1)*qk(2:nxp1,1:ny)    & 
            & - aC(2:nxp1,2:nyp1)*qk(2:nxp1,2:nyp1) 
        
        tksk = sum(tk*sk)
        tktk = sum(tk*tk)
        omega = tksk/(tktk+1.e-30)
        ps = ps + alpha*zk + omega*qk
        res = sk - omega*tk

        ! check norm
        rnorm = maxval( dabs(res) )  
        write(*,100) l, rnorm 

        if( rnorm .lt. resmax )then
            ps = ps - minval(ps) !If Neumann BC's
            ! write(*,100) l, rnorm 
            return
        endif

        100 Format('nIter = ', 1I6, 4X, 'Error = ', 1E15.8)

    end do

return
end subroutine ppe_bicgstab2D
