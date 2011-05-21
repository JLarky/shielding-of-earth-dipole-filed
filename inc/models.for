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
      SUBROUTINE MODEL_CF_RC_T3 
     _(ID, A, XI, F, DER, IA, NTOT, NLIN, NNON,INDEPVAR, NDEPVAR, ICO)
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
	 IMPLICIT  REAL * 8  (A - H, O - Z)
C
         INTENT(IN)    :: ID, A, XI, IA, NTOT, NLIN
         INTENT(IN)    :: NNON, INDEPVAR, NDEPVAR, ICO
         INTENT(OUT)   :: F, DER

         DIMENSION  F(NDEPVAR),DER(NDEPVAR,NTOT),A(NTOT),XI(INDEPVAR),
     _     IA(NTOT)
c
      COMMON /FPAR/ XAPPAD

c     15,13,2,8,3
      if ((NTOT.ne.15).or.(NLIN.ne.13).or.(NNON.ne.2).or.
     _            (INDEPVAR.ne.8).or.(NDEPVAR.ne.3)) then
         print *, 'function MODELVEC was called with impropriet args'
         stop 1
      end if
         

c     INPUT

         X=XI(1)
         Y=XI(2)
         Z=XI(3)
         PS=XI(4)
         PD=XI(5)
         W1=XI(6)
         W2=XI(7)
         W3=XI(8)
         SPS=DSIN(PS)
         
	XAPPAD=A(14)
        RCSCALE=A(15)

       CALL FIELDS (X,Y,Z,PS,PD,RCSCALE,CFX,CFY,CFZ,RCX,RCY,RCZ,
     * T1X,T1Y,T1Z,T2X,T2Y,T2Z,T3X,T3Y,T3Z)
     
        CPD=sqrt(sqrt(PD/2.0))
c       BX=a(1)*CFX+(a(2)+a(12)*W2+a(13)*W3)*RCX+(a(3)+a(6)*CPD+ a(9)*W1)*T1X
c                                               +(a(4)+a(7)*CPD+a(10)*W1)*T2X
c                                               +(a(5)+a(8)*CPD+a(11)*W1)*T3X

        DER(1, 1)=CFX
        DER(2, 1)=CFY
        DER(3, 1)=CFZ
        DER(1, 2)=RCX
        DER(2, 2)=RCY
        DER(3, 2)=RCZ
        DER(1, 3)=T1X
        DER(2, 3)=T1Y
        DER(3, 3)=T1Z
        DER(1, 4)=T2X
        DER(2, 4)=T2Y
        DER(3, 4)=T2Z
        DER(1, 5)=T3X
        DER(2, 5)=T3Y
        DER(3, 5)=T3Z
        DER(1, 6)=T1X*CPD
        DER(2, 6)=T1Y*CPD
        DER(3, 6)=T1Z*CPD
        DER(1, 7)=T2X*CPD
        DER(2, 7)=T2Y*CPD
        DER(3, 7)=T2Z*CPD
        DER(1, 8)=T3X*CPD
        DER(2, 8)=T3Y*CPD
        DER(3, 8)=T3Z*CPD
        DER(1, 9)=T1X*W1
        DER(2, 9)=T1Y*W1
        DER(3, 9)=T1Z*W1
        DER(1,10)=T2X*W1
        DER(2,10)=T2Y*W1
        DER(3,10)=T2Z*W1
        DER(1,11)=T3X*W1
        DER(2,11)=T3Y*W1
        DER(3,11)=T3Z*W1
        DER(1,12)=RCX*W2
        DER(2,12)=RCY*W2
        DER(3,12)=RCZ*W2
        DER(1,13)=RCX*W3
        DER(2,13)=RCY*W3
        DER(3,13)=RCZ*W3
        


        F(1)=0.d0               !     COMPONENTS OF THE TOTAL EXTERNAL FIELD
        F(2)=0.d0
        F(3)=0.d0

        DO I=1,NLIN
        F(1)=F(1)+A(I)*DER(1,I)
        F(2)=F(2)+A(I)*DER(2,I)
        F(3)=F(3)+A(I)*DER(3,I)
        END DO

      RETURN
      END
C



cOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
c
      SUBROUTINE MODEL_CF_RC_T3_CUSP 
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
	 IMPLICIT  REAL * 8  (A - H, O - Z)
C
         INTENT(IN)    :: ID, A, XI, IA, NTOT, NLIN
         INTENT(IN)    :: NNON, INDEPVAR, NDEPVAR, ICO
         INTENT(OUT)   :: F, DER

         DIMENSION  F(NDEPVAR),DER(NDEPVAR,NTOT),A(NTOT),XI(INDEPVAR),
     _     IA(NTOT)
c
      COMMON /par/
     *  ampl0,ampl1,dp0,c,Dtheta,PhiC0,alpha1,alpha2,YC0,powf,powth
      COMMON /FPAR/ XAPPAD

c     20,16,4,10,3
      if ((NTOT.ne.20).or.(NLIN.ne.16).or.(NNON.ne.4).or.
     _            (INDEPVAR.ne.10).or.(NDEPVAR.ne.3)) then
         print *, 'function MODELVEC was called with impropriet args'
         stop 1
      end if


c     INPUT

         X=XI(1)
         Y=XI(2)
         Z=XI(3)
         PS=XI(4)
         PD=XI(5)
         BZimf=XI(6)         
         W1=XI(7)
         W2=XI(8)
         W3=XI(9)
         W5=XI(10)
         SPS=DSIN(PS)


      ampl0=-1.5
      ampl1=-0.5
      dp0=0.04
      c=0.5
      Dtheta=0.6
c      Phic0=0.24
c      Phic0=A(19)
      Phic0=A(19)+A(20)*W5
      alpha1=0.13
      alpha2=0.03
      YC0=5
      powf=5
      powth=5

	XAPPAD=A(17)
        RCSCALE=A(18)


       CALL FIELDS (X,Y,Z,PS,PD,RCSCALE,CFX,CFY,CFZ,RCX,RCY,RCZ,
     * T1X,T1Y,T1Z,T2X,T2Y,T2Z,T3X,T3Y,T3Z)
       CALL CUSP_DEPRESSION (PS,X,Y,Z,C_BX,C_BY,C_BZ)
     
	SPD=sqrt(PD/2.0)
        CPD=sqrt(SPD)
        SBZ=BZimf/5.0
        
c       BX=a(1)*CFX+(a(2)+a(12)*W2+a(13)*W3)*RCX+(a(3)+a(6)*CPD+ a(9)*W1)*T1X
c                                               +(a(4)+a(7)*CPD+a(10)*W1)*T2X
c                                               +(a(5)+a(8)*CPD+a(11)*W1)*T3X
c						+(a(14)+a(15)SPD+a(16)*SBZ)*CUSP_X

        DER(1, 1)=CFX
        DER(2, 1)=CFY
        DER(3, 1)=CFZ
        DER(1, 2)=RCX
        DER(2, 2)=RCY
        DER(3, 2)=RCZ
        DER(1, 3)=T1X
        DER(2, 3)=T1Y
        DER(3, 3)=T1Z
        DER(1, 4)=T2X
        DER(2, 4)=T2Y
        DER(3, 4)=T2Z
        DER(1, 5)=T3X
        DER(2, 5)=T3Y
        DER(3, 5)=T3Z
        DER(1, 6)=T1X*CPD
        DER(2, 6)=T1Y*CPD
        DER(3, 6)=T1Z*CPD
        DER(1, 7)=T2X*CPD
        DER(2, 7)=T2Y*CPD
        DER(3, 7)=T2Z*CPD
        DER(1, 8)=T3X*CPD
        DER(2, 8)=T3Y*CPD
        DER(3, 8)=T3Z*CPD
        DER(1, 9)=T1X*W1
        DER(2, 9)=T1Y*W1
        DER(3, 9)=T1Z*W1
        DER(1,10)=T2X*W1
        DER(2,10)=T2Y*W1
        DER(3,10)=T2Z*W1
        DER(1,11)=T3X*W1
        DER(2,11)=T3Y*W1
        DER(3,11)=T3Z*W1
        DER(1,12)=RCX*W2
        DER(2,12)=RCY*W2
        DER(3,12)=RCZ*W2
        DER(1,13)=RCX*W3
        DER(2,13)=RCY*W3
        DER(3,13)=RCZ*W3
        DER(1,14)=C_BX
        DER(2,14)=C_BY
        DER(3,14)=C_BZ
        DER(1,15)=C_BX*SPD
        DER(2,15)=C_BY*SPD
        DER(3,15)=C_BZ*SPD
        DER(1,16)=C_BX*SBZ
        DER(2,16)=C_BY*SBZ
        DER(3,16)=C_BZ*SBZ
        


        F(1)=0.d0               !     COMPONENTS OF THE TOTAL EXTERNAL FIELD
        F(2)=0.d0
        F(3)=0.d0

        DO 1 I=1,NLIN
        F(1)=F(1)+A(I)*DER(1,I)
        F(2)=F(2)+A(I)*DER(2,I)
   1    F(3)=F(3)+A(I)*DER(3,I)

      RETURN
      END
C


