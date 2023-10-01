!======================================================================================
!             solve the possion equation of stream function here
!--------------------------------------------------------------------------------------
subroutine Stream_Function
    use flow_var, only : iPossion_Equation
    implicit none

    ! solve possion equation 
    select case(iPossion_Equation)
    case(1)
        call Stream_Function_SOR
    case(2)
        ! call Stream_Function_SIP
    case(3)
        ! call Stream_Function_BiCGStab
    end select

end subroutine Stream_Function
