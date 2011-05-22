      program test
      implicit none

      print *, 'var_r'
      call test_var_r
      print *, 'dipole'
      call test_dipole

      contains
        subroutine assert(a,b)
           real*8  :: a,b
           if (1d-7 .lt. (abs(a-b)/(abs(a)+abs(b)+1d-7))) then
              print *, a,'=/=',b
              stop 1
           end if
        end subroutine
        subroutine test_var_r
          real*8 :: theta, r, r_0, var_r, tested_r(21)
          data tested_r /
     _     10.281916970713286 , 10.59600833616862 , 10.937382777323126 ,
     _     11.310124655377043, 11.719217999982913 , 12.170816490619728 ,
     _     12.67261902950815 , 13.234403264701184 , 13.868801764862072 ,
     _     14.592462511112778, 15.427839922922313 , 16.406063957908475 ,
     _     17.57174504791317 , 18.99146717248451 , 20.76984091800765 ,
     _     23.082563576099566, 26.25274109038196 , 30.957868786141684 ,
     _     38.949834584790224, 57.02884687466857 , 234.50355590680456 /
          integer :: i

          r_0 = 10.1
          i = 0
          do r = 10, -10, -1.
             i = i+1
             theta = acos(r/r_0)
!             print *, theta, var_r(theta)
             call assert(var_r(theta),tested_r(i))
          end do
        end subroutine
        subroutine test_dipole
          use utils
          real*8 :: m(3), r(3), d(3)
          real*8 :: XGSW,YGSW,ZGSW,BXGSW,BYGSW,BZGSW
          m(1:3) = 0.; m(3) = -30000.
          r(1:3) = 0.; r(1) = -10.
          d = dipole(m, r)
          call assert(d(1), 0.d0)
          call assert(d(2), 0.d0)
          call assert(d(3), 30.d0)
!          print *, d
          xGSW = r(1);yGSW = r(2);zGSW = r(3);
          call recalc_08 (2008,106,4,15,0,-400.d0,0.d0,0.d0)
          call DIP_08 (XGSW,YGSW,ZGSW,BXGSW,BYGSW,BZGSW)
!          print *, BXGSW,BYGSW,BZGSW
        end subroutine
      end program
