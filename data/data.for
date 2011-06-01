      program data
        use utils
        implicit none
        integer :: fd1, fd2, fd3
        logical :: debug
        parameter(fd1 = 10)
        parameter(fd2 = 11)
        parameter(fd3 = 12)
!        parameter(DEBUG = .True.)
        parameter(DEBUG = .False.)
        integer :: N, i
        real*8  :: theta, r, dr, old_dr
        type(point) :: r_0, old_r
        type(point), allocatable :: subset(:)

!     subsolar point (theta = 0) 
        allocate(subset(1))
        subset = get_subset(0d0,1)
        r_0 = subset(1)
        deallocate(subset)
        print *, 'r_0', r_0

        open(fd1, file='data.dat')
        open(fd2, file='points.dat')
        open(fd3, file='XZ.dat')
        N = 60                  ! number of points with same r (distance from Earth)
        if (debug) then
        write (fd1, *) 2, 'debug'
        else
        write (fd1, *) N, 'number of points in subset'
        end if

!     initial subset with X \approx -100 Re
        theta = 2.85
        allocate(subset(N))
        subset = get_subset(theta,N)
        old_r = subset(1)
        dr = dist(subset(1), subset(2))
        call save_subset(subset, theta)

 1      do ! this loop will drop subsets which are too close with each other
           ! too close means r1 and r2 closer that 1Re
           ! or closer that dr, which is distanse between closest points in subset

           theta = theta -0.0005
           if (theta.lt.0) then
              goto 2 ! exit.
           end if

           subset = get_subset(theta,N)
           dr = dist(subset(1), subset(2))
!           print *, 'dr', dr

           if (dist(subset(1), old_r) .le. 1*max(dr, 1.)) then
              goto 1
           end if

           call save_subset(subset, theta)
           old_r = subset(1)

!           r=sqrt(subset(1)%r(1)**2+subset(1)%r(2)**2+subset(1)%r(3)**2) ! dist from (0,0,0)
!           print *, r
           r = dist(r_0, subset(1)) ! dist from r_0
           print *, r
        end do

 2      continue

        close(fd1)
        close(fd2)
        close(fd3)

      contains
        function dist(r1,r2)
          type(point) :: r1,r2
          real*8 :: dist
          integer :: i
          dist = 0.
          do i = 1, 3
             dist = dist + (r1%r(i)-r2%r(i))**2
          end do
          dist = sqrt(dist)
        end function
        function get_subset(theta, N)
          intent(in) :: theta, N
          integer :: N, i
          type(point) :: get_subset(N)
          real*8  :: theta, r_n, var_r, r, alpha
          real*8  :: nx, nz
          ! r = R(\theta). nx, nz are normals
          call var_r_and_n(theta, r, nx, nz)

          get_subset(1:N)%r(1)=r*cos(theta) ! r_x
          get_subset(1:N)%n(1)=nx ! n_x
          r_n  = r*sin(theta)
          do i = 1, N
             alpha = (i-1.0)/n*2*3.1415928
             get_subset(i)%r(2)=r_n*sin(alpha)
             get_subset(i)%r(3)=r_n*cos(alpha)
             get_subset(i)%n(2)=nz*sin(alpha)
             get_subset(i)%n(3)=nz*cos(alpha)
          end do
        end function
        subroutine save_subset(subset, theta)
          type(point) :: subset(:)
          integer :: i
          real*8 :: theta
          if (debug) then
          write(fd1,*) subset(1),subset(16)
          else
          write(fd1,*) subset
          end if
          do i = 1, size(subset)
             write(fd2,*) subset(i)%r, subset(i)%n, theta
             if ((i.eq.1).or.(2*(i-1).eq.N)) then
                write(fd3,*) subset(i)%r, subset(i)%n, theta
             end if
          end do
        end subroutine
      end program
