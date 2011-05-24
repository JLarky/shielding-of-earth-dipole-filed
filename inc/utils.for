      module utils
        type point
          real*8 :: r(3), n(3)
        end type
      contains
      subroutine var_r_and_n(theta, r, nx, nz)
        implicit none
        intent (in)  :: theta
        intent (out) :: r, nx, nz
        real*8       :: theta, r, nx, nz
        real*8       :: er_x, er_z, et_x, et_z, dr_dt, norm
        real*8, save :: r_0, alpha, b_z, p_dyn

        if (r_0.eq.0) then ! initialization
           b_z = 0.
           p_dyn = 2.
           r_0 = (10.22+1.29*tanh(0.184*(B_z+8.14)))*p_dyn**(-0.1515152)
           alpha = (0.58-0.007*B_z)*(1+0.024*log(p_dyn))
        end if
        r = r_0*(2/(1+cos(theta)))**alpha

        er_x = cos(theta)
        er_z = sin(theta)
        et_x = -sin(theta)
        et_z = cos(theta)
        dr_dt = r*alpha*sin(theta)/(1.+cos(theta))
        nx = er_x - (dr_dt/r)*et_x
        nz = er_z - (dr_dt/r)*et_z
        norm = sqrt(1. + (dr_dt/r)**2)
        nx = nx/norm
        nz = nz/norm
      end subroutine

!     \vec B = \frac{3(m.r)r-m r^2} {r^5}
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
