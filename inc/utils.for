      function var_r(theta)
        implicit none
        intent (in)  :: theta
        real*8       :: theta, var_r
        real*8, save :: r_0, alpha, b_z, p_dyn

        if (r_0.eq.0) then ! initialization
           b_z = 0.
           p_dyn = 2.
           r_0 = (10.22+1.29*tanh(0.184*(B_z+8.14)))*p_dyn**(-0.1515152)
           alpha = (0.58-0.007*B_z)*(1+0.024*log(p_dyn))
        end if
        var_r = r_0*(2/(1+cos(theta)))**alpha
      end function

      module utils
      contains
      function dipole(m, r)
        implicit none
        intent (in)  :: m, r
        real*8       :: m(3), r(3), dipole(3), s_m, mod_2, mod_5
        integer      :: i, j

        mod_2 = 0. ! r**2
        do i = 1, 3
           mod_2 = mod_2 + r(i)**2
        end do
!        mod = sqrt(mod_2)
        mod_5 = mod_2**2*sqrt(mod_2) ! r**5

        do i = 1, 3
           s_m = 0.
           do j = 1, 3
              s_m = s_m + m(j)*r(j)
           end do
           dipole(i) = (3*s_m*r(i)-m(i)*mod_2 )/mod_5
       end do
      end function
      end module
