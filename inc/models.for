      module jl_models
       type jl_model
        integer :: NLIN,NNON,NTOT,INDEPVAR,NDEPVAR
        real*8,allocatable :: A(:)
       end type
       type (jl_model) :: modpar
       contains
        function get_a(NTOT, A_file)
          INTEGER,   intent(in) :: NTOT
          CHARACTER, intent(in) :: A_file*(*)
          real*8                :: get_a(NTOT)
          open (unit=1,file=A_file)
          read (1,'(G15.6)') get_A
          close(1)
        end function
      end module



cOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
c
      SUBROUTINE MODEL_DIP_SH
     _(ID, A, XI, F, DER, IA, NTOT, NLIN, NNON, INDEPVAR, NDEPVAR, ICO)
C
C        ***  N.A. Tsyganenko ***  28.10.1997  ***
C
C  Calculates dependent model variables and their derivatives with respect to
c  both linear and nonlinear parameters, to be used in the RMS minimization procedure.
C
C      Description of parameters:
C
C  ID  - number of a current point from the data set (initial assignments can
c        be made for ID=1 only, saving thus CPU time)
C  A   - input vector containing model parameters;
C  XI  - input vector containing independent variables;
C  F   - output double precision vector containing
C        calculated values of dependent variables;
C  DER   - output double precision vector containing
C        calculated values for derivatives of dependent
C        variables with respect to model parameters;
C  IA -  integer vector; IA(L)=0 means that the Lth parameter should be kept
c        fixed in the search, and hence there is no need in calculating
c        derivatives of field components with regard to this parameter;
c        IA(L)=1 means that the derivatives are needed.
c
c  NTOT - total number of model parameters (NTOT=NLIN+NNON)
c
c
C  ICO - input parameter coded as follows:
c
c       ICO = 0 - derivatives with respect to nonlinear parameters must be calculated;
c       ICO = 1 - only derivatives with respect to linear parameters must be calculated
C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C
      use utils
	 implicit none
C
         INTENT(IN)    :: ID, A, XI, IA, NTOT, NLIN
         INTENT(IN)    :: NNON, INDEPVAR, NDEPVAR, ICO
         INTENT(OUT)   :: F, DER

         integer :: i, id, ico, NTOT, NLIN, NNON, INDEPVAR, NDEPVAR
         real*8 :: x,y,z,nx,ny,nz,df_dn
         real*8 :: F(NDEPVAR),DER(NDEPVAR,NTOT),A(NTOT),XI(INDEPVAR),
     _     IA(NTOT)
         real*8 :: b(NNON), fi(NNON), U(NNON,3)
c

c     10,5,5,6,1
      if ((NTOT.ne.10).or.(NLIN.ne.5).or.(NNON.ne.5).or.
     _            (INDEPVAR.ne.6).or.(NDEPVAR.ne.1)) then
         print *, 'function MODELVEC was called with impropriet args'
         stop 1
      end if

c     INPUT

!      print *, 'xi', xi

         X=(XI(1))
         Y=XI(2)
         Z=XI(3)
         nx=XI(4)
         ny=XI(5)
         nz=XI(6)

         b(1:NNON) = A(NLIN+1:NNON)+.1 ! to del
!      print *, 'anlin', a(NLIN+1:NTOT),'b', b
!         print *, 'b', b

         U(1:NLIN,1:3) = 0.
         do i = 1, NLIN
            U(i,1)=U(i,1)+exp(sqrt(2.)*b(i)*x)
     _           *cos(b(i)*y)*sin(b(i)*z)*sqrt(2.)*b(i)
            U(i,2)=U(i,2)+exp(sqrt(2.)*b(i)*x)
     _           *(-sin(b(i)*y))*sin(b(i)*z)*b(i)
            U(i,3)=U(i,3)+exp(sqrt(2.)*b(i)*x)
     _           *cos(b(i)*y)*cos(b(i)*z)*b(i)
         end do
!         print *, 'fi', fi, '<'

!         print *, 'a', a
!         print *, 'f', f1, f2
!         print *, 'd', df_dn

        F(1)=0.d0
        saved_f(1:3) = 0.
        DER(1, 1:NLIN) = 0.

        do I=1,NLIN
!          print *, 'a', i, a(i)
           ! - \grad U
           saved_f = saved_f -A(I)*U(I,1:3)

           if (abs(nx).gt.0.) then
              DER(1, i) = DER(1, i)+U(I,1)/nx
           end if

           if (abs(ny).gt.0.) then
              DER(1, i) = DER(1, i)+U(I,2)/ny
           end if

           if (abs(nz).gt.0.) then
              DER(1, i) = DER(1, i)+U(I,3)/nz
           end if

           F(1)=F(1)+A(I)*DER(1,I)
!           print *
!           print *, nx,ny,nz
!           print *, saved_f(1:3)
!           print *, U(i,1:3)
!           print *, DER(i,1)
!           print *, F(1)
!           if (i.gt.1) then; stop; end if
        end do

!        print *, 'f', x, y, z, f(1)
!        print *, 'f', x, y, z, der(1,1)

!        print *, df_dn, b1, f1, df_dn*b1*f1
!        print *, df_dn, b2, f2, df_dn*b2*f2

      RETURN
      END
C


