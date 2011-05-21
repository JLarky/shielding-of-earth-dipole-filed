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
