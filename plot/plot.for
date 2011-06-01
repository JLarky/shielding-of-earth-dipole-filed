      module shared
       implicit none
        INTEGER NTOT,NLIN,NNON,INDEPVAR,NDEPVAR ! call modfun
        REAL*8,allocatable :: A(:) !(NTOT)

      end module

      module t_model
        type Tmodel
          INTEGER NTOT,NLIN,NNON,INDEPVAR,NDEPVAR
        end type Tmodel
        contains
          function MODFUN(A)
            MODFUN=A
          end function
      end module t_model

      module p_models
       contains
        subroutine md_T02_CUSP(PS,X,Y,Z,BX,BY,BZ) !(A,NTOT,NLIN,NNON,INDEPVAR,NDEPVAR,modfun)
          IMPLICIT REAL*8 (A-H,O-Z)
          print *, PS,X,Y,Z,BX,BY,BZ
        end subroutine
        subroutine md_DIP_08(x_notused_PS,X,Y,Z,BX,BY,BZ)
          IMPLICIT REAL*8 (A-H,O-Z)
          intent(in)  :: x_notused_PS,X,Y,Z;
          intent(out) :: BX,BY,BZ
          COMMON /GEOPACK1/ AB(10),SPS,CPS,ABB(3),PS,CD(18) ! call recalc_08
          call recalc_08 (2008,47,04,45,0,-400.d0,0.d0,0.d0) !  HERE JUST TO SPECIFIY THE EARTH'S DIPOLE MOMENT
          !PS=PSI*0.01745329
          PS=x_notused_PS
          SPS=DSIN(PS);CPS=DCOS(PS)
          call DIP_08 (X,Y,Z,BX,BY,BZ)
        end subroutine
      end module


      program plot
      use jl_models
      use cusp
      IMPLICIT NONE
      EXTERNAL model_dip_sh_fp

!      call plot_trace(20,5,15,6,1,model_dip_sh_fp,
!      call plot_trace(40,10,30,6,1,model_dip_sh_fp,
!      call plot_trace(32,8,24,6,1,model_dip_sh_fp,
!      call plot_trace(10,5,5,6,1,model_dip_sh_fp,
      call plot_trace(16,8,8,6,1,model_dip_sh_fp,
     _ '../fitting/model_par.dat')

c      call plot_trace(20,16,4,10,3,TR_T02_CUSP,
c     _ '../out/model_par_m10-10.dat')

c      call p_profile(20,16,4,10,3,TR_T02_CUSP,
c     _ '../out/model_par_m10-10.dat', 'profile.dat')

c      call plot_trace(15,13,2,8,3,TR_T02,
c     _ '../out/model_par_back_m10-10.dat')

c      call p_profile(15,13,2,8,3,TR_T02,
c     _ '../out/model_par_back_m10-10.dat', 'profile_b.dat')

      end program

      subroutine p_profile(NTOT,NLIN,NNON,INDEPVAR,NDEPVAR,modfun,
     _ A_file, O_file)
      use shared, ntot_=>ntot,nlin_=>nlin,nnon_=>nnon,
     _ indepvar_=>indepvar,ndepvar_=>ndepvar
      use p_models
      use jl_models

      implicit real*8 (a-h,o-z)
      CHARACTER   A_file*(*),O_file*(*)

      ntot_=ntot;nlin_=nlin;nnon_=nnon;
      indepvar_=indepvar;ndepvar_=ndepvar;

      if (allocated(a)) deallocate(a)
      allocate(a(ntot))
      a=get_a(ntot,a_file)
      R=6.0
      do I=-90,90
	  Phi=3.14159*I/180.0
	  X=R*DCOS(Phi)
	  Y=0.0
	  Z=R*DSIN(Phi)

	  PS=.0
	  call modfun (IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
          call md_DIP_08(PS,X,Y,Z,BXGSW,BYGSW,BZGSW)
c	  call TR_T02 (IOPT,PARMOD,PS,X,Y,Z,BXGSW,BYGSW,BZGSW)

          open (unit=10,file=O_file)
94        FORMAT(3F9.2,A,3F9.2,A,3F9.2,A,3F9.2)
	  DD=SQRT(BXGSW**2+BYGSW**2+BZGSW**2)
	  WD=SQRT((BXGSW+BX)**2+(BYGSW+BY)**2+(BZGSW+BZ)**2)
	  write (10,94) X,Y,Z,'    MF= ', BX,BY,BZ, '    DF= ', 
     _	  BXGSW,BYGSW,BZGSW, '    ABS= ', DD,WD,DD-WD
      end do
      print *, 'profile done'
      close(10)
      end subroutine p_profile


      subroutine plot_trace(NTOT_,NLIN_,NNON_,INDEPVAR_,NDEPVAR_,modfun,
     _ A_file)
      use shared
      use jl_models
      use utils
      IMPLICIT REAL * 8 (A - H, O - Z)
      CHARACTER         A_file*(*)
      EXTERNAL modfun

      CHARACTER*12 LIN_OUT,NUM_OUT
      CHARACTER*12 C,D

      DIMENSION XC(5000),YC(5000),ZC(5000),STOREX(10000),STOREY(10000),
     + STOREZ(10000),NUMPNT(200)
c

      COMMON /GEOPACK1/ AB(10),SPS,CPS,ABB(3),PS,CD(18) ! call recalc_08


      EXTERNAL DIP_08, null_field

      NTOT=NTOT_;NLIN=NLIN_;NNON=NNON_;
      INDEPVAR=INDEPVAR_;NDEPVAR=NDEPVAR_;

      if (allocated(a)) deallocate(a)
      allocate (A(NTOT))
      print *, 1

      A=get_a(NTOT, A_file)
      print *, 1, a

      PRINT *, A

      call recalc_08 (2008,47,04,45,0,-400.d0,0.d0,0.d0) !  HERE JUST TO SPECIFIY THE EARTH'S DIPOLE MOMENT


      LIN_OUT= 'lines.dat'
      NUM_OUT='numpnt.dat'

       PRINT *, '    ENTER TILT ANGLE (DEGREES)'
c       read *, PSI
	PSI=0.0
       PS=PSI*0.01745329
       SPS=DSIN(PS)
       CPS=DCOS(PS)

	goto 2333
c	goto 2334
345     print *, 'Put x,y,z'
c	read *, x, y, z
	x=1.0;y=2.0;z=3.0;
	call modfun (IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
	print *, BX,BY,BZ
	stop
	goto 345

 2333   continue
c	goto 999

2334   DIST=90.
      Y=0.

      INDD=0
C
      DO 1 I=1,31
      T=59.+I*1.
      TE=T*0.01745329
      COST=COS(TE)
      DO 1 K=1,2
      IF (I.EQ.31.AND.K.EQ.2) GOTO 2
      F=180.*(K-1)
      FI=F*0.01745329
      X=COST*COS(FI)
      Z=SIN(TE)
      XI=X*CPS+Z*SPS
      ZI=Z*CPS-X*SPS
      print *, 1
c
      CALL TRACE_08(XI,Y,ZI,1.d0,1.0d0,1.d-5,DIST,1.d0,IOPT,PARMOD,
!     * modfun,null_field,XF,YF,ZF,XC,YC,ZC,LP,5000)
     * modfun,DIP_08,XF,YF,ZF,XC,YC,ZC,LP,5000)
c
        PRINT *, '  I=', I, ' (OF 31)'
c
        ILINE=ILINE+1
        NUMPNT(ILINE)=LP
          DO 22 N=1,LP
           IPOINT=IPOINT+1
           STOREX(IPOINT)=XC(N)
           STOREY(IPOINT)=YC(N)
   22      STOREZ(IPOINT)=ZC(N)
C
      IF(ABS(XF)+ABS(ZF).GT.5.) GOTO 1
      XMI=1.
C
      IF (K.EQ.1.AND.INDD.EQ.0) THEN
         DO 1380 ILM=1,LP
         IF (XC(ILM).LT.XMI) XMI=XC(ILM)
 1380    CONTINUE
         IF (XMI.GT.-5.) GOTO 110
            INDD=1
            GOTO 113
      ENDIF
      IF (INDD.EQ.1.AND.K.EQ.1) GOTO 1
C
C   In some cases the first open line for K=1 (a bit poleward from the
C   cusp) appears to be closed. In such a case XF1 could be assigned
C   to the nightside field line, rather than to the dayside one (as it
C   must be). To avoid this, we control the minimal value of the X coor-
C   dinate (XMI).  XMI>-10 and K=1 means that the lines starting at noon
c   meridian still belong entirely to the dayside sector.
C
      IF (INDD.EQ.1) GOTO 1
 113  XF2=XF
      ZF2=ZF
      GO TO 1
 110  XF1=XF
      ZF1=ZF
C
C   The above quantities (XF1,ZF1 and XF2,ZF2) specify the coordinates of
C   the landing points in the southern hemisphere corresponding, respectively,
C   to the dayside and nightside last closed field lines.
C
   1  CONTINUE
   2  CONTINUE

      print *, '  xf1,xf2,zf1,zf2,sps,cps=',xf1,xf2,zf1,zf2,sps,cps

      XL1=ATAN2((XF1*SPS+ZF1*CPS),(XF1*CPS-ZF1*SPS))
      XL2=ATAN2((XF2*SPS+ZF2*CPS),(XF2*CPS-ZF2*SPS))
      DL=XL1-XL2
      SDL=1.d0*DL/.01745329 !1.d0 means 1 degree interval in footpoint latitude
      LL=SDL-0.2
       IF (LL.LT.1) GOTO 13
      DTET=DL/(LL+1)
C
      DO 11 I=1,LL
      TE=XL1-I*DTET
      X=COS(TE)
      Z=SIN(TE)
      XI=X*CPS+Z*SPS
      ZI=Z*CPS-X*SPS

      CALL TRACE_08 (XI,Y,ZI,-1.d0,1.0d0,1.d-5,DIST,1.d0,IOPT,PARMOD,
     * modfun, DIP_08,XF,YF,ZF,XC,YC,ZC,LP,5000)

       PRINT *, '  I=', I, ' (OF ',LL,')'
C
        ILINE=ILINE+1
        NUMPNT(ILINE)=LP
          DO 33 N=1,LP
           IPOINT=IPOINT+1
           STOREX(IPOINT)=XC(N)
           STOREY(IPOINT)=YC(N)
   33      STOREZ(IPOINT)=ZC(N)
C
  11  CONTINUE
C
  13   CONTINUE
C
            OPEN(UNIT=1,FILE=LIN_OUT)
            WRITE(1,888) (STOREX(N),STOREY(N),STOREZ(N),N=1,IPOINT)
  888       FORMAT(3(1X,F7.2))
            CLOSE(UNIT=1)
C
           DO I=1,ILINE
             I1=ILINE+2-I
             I2=I1-1
             NUMPNT(I1)=NUMPNT(I2)
           END DO
C
             NUMPNT(1)=ILINE
C
             OPEN(UNIT=1,FILE=NUM_OUT)
             WRITE(1,777) NUMPNT
             CLOSE(UNIT=1)
 777         FORMAT(I3)
C
c
999   END subroutine plot_trace
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c
c

c     null field
      subroutine null_field (X,Y,Z,BX,BY,BZ)
      bx = 0
      by = 0
      bz = 0
      end subroutine

c     новые каспы
      subroutine model_dip_sh_fp (IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
        use cusp
        use shared
        use utils
        implicit none
        intent(in)  :: IOPT,PARMOD,PS,X,Y,Z
        intent(out) :: BX,BY,BZ
        REAL*8  :: F(NDEPVAR),DER(NDEPVAR,NTOT),XI(INDEPVAR)
        integer :: IA(NTOT)
        real*8  :: parmod, ps, x, y, z, bx, by, bz
        integer :: iopt
          
       XI(1)=X !X=XI(1)
       XI(2)=Y !Y=XI(2)
       XI(3)=Z !Z=XI(3)
       XI(4)=PS !PS=XI(4)
       XI(5)=PS !PS=XI(4)
       XI(6)=PS !PS=XI(4)

       saved_f(1) = 1.
!       print *, saved_f
!        print *, 'a', a, size(a)
        call model_dip_sh(0,A,XI,F,DER,IA,
     _                        NTOT,NLIN,NNON,INDEPVAR,NDEPVAR,0)
!        print *, xi
!        print *, saved_f

!       f(1:3)=1*f(1:3)
        f = saved_f
      BX=f(1)
      BY=f(2)
      BZ=f(3)
      end subroutine
