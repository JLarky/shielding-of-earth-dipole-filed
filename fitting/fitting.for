c
      P r o g r a m   RUNMODEL
C
c----------------------------------------------------------------------
c
      use utils
      IMPLICIT  REAL * 8  (A - H, O - Z)

      CHARACTER*100 DATANAME
      CHARACTER*100 OUTFNAME
      CHARACTER*13 COMMENT
      CHARACTER*7 SCNAME

      CHARACTER*80 NAMEPARAM,NAMESIMPLEX, NAMEPAR, NAMESIM

c----- new string and integer variables for reading sequential subsets of data -----
      character*4 seqstring
      integer*4 seqnum
c-----------------------------------------------------------------------------------


      PARAMETER (NLIN=5,NNON=5,NNON1=NNON+1,NTOT=NLIN+NNON)
      PARAMETER (NUMVTOT=7)       ! NUMVTOT=INDEPVAR+NDEPVAR (adjust in each individual case)
      PARAMETER (NPNT_UPP=1000000)  ! UPPER LIMIT ON THE NUMBER OF DATA POINTS
      PARAMETER (N_ELEM_B=NUMVTOT*NPNT_UPP) ! =NUMVTOT*NPNT_UPP  -
      PARAMETER (NDEPVAR=1)
c       1 normal field component
C
c  ATENTION:  Some of above parameters presents in common block /siplex_param/
c                in order to use in simplex.for's functions SRMS,FUNK,SIMPLEX,AMOEBA,AMOTRY
C
      PARAMETER (INDEPVAR=6)             !   ADJUSTED TO X,Y,Z,nx,ny,nz
C      independent variables (arguments): in this case, we have three coordinates,
c        tilt angle
C

      PARAMETER  (DATANAME='../data/data.dat')   !   ADJUST, WHEN USING A NEW DATA SET
      PARAMETER (OUTFNAME="model_par.dat")


      real*8  :: m_dip(3), dip_b(3)
      integer :: subset_len, subset_i, s_i
      type(point) :: cp
      type(point), allocatable :: subset(:)



C    IN THIS PARTICULAR CODE, THE ORDER OF DATA IN A RECORD IS AS FOLLOWS:
C                      XGSM,YGSM,ZGSM,PS,BxGSM,ByGSM,BzGSM

      DIMENSION B(NUMVTOT,NPNT_UPP),BBX(NPNT_UPP),BBY(NPNT_UPP),
     * BBZ(NPNT_UPP), BBXM(NPNT_UPP),BBYM(NPNT_UPP),BBZM(NPNT_UPP),
     * WEIGHT(NPNT_UPP)

      DIMENSION P(NNON1,NNON),XI(INDEPVAR),F(NDEPVAR),DER(NDEPVAR,NTOT)
      DIMENSION  A(NTOT),XXXXX(NNON)
      DIMENSION SCORE(200),BAV(200)

      common /q/ Q,Q_old,B_rms
      COMMON /TOLERANCE/ TOL
      COMMON /NAMES/ NAMEPARAM,NAMESIMPLEX

      COMMON /GEOPACK1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,
     * SPS,CPS,DS3,CGST,SGST,PSI,A11,A21,A31,A12,A22,A32,A13,A23,A33,
     * E11,E21,E31,E12,E22,E32,E13,E23,E33
c
      EXTERNAL MODEL_DIP_SH ! dipole shield
c     You must change modelname few times in this file
C
C            The array B contains experimental data in the following order:
c     each (B(i,*),i=1,NUMVTOT) contains:   x,y,z,ps,pd,bx,by,bz,..
c       (see right after the   33 READ ....  statement below)
c
c   A is the vector of independent model parameters (linear+nonlinear ones)
c
c   The array B  contains experimental data
c
      DIMENSION IA(NTOT),NUM(NTOT)
c
c  the integer array IA contains information on which of the parameters
c  should be treated as variable: IA(L)=1 means that the L-th parameter
c  is free and will be fitted to data; IA(L)=0 means that the parameter
c   is fixed
c  the array NUM stores the numbers of fixed parameters which are entered
c   at prompt from the keyboard
C
      DATA  IA /NTOT*1/
c

      COMMON /simplex_param/ NLINc,NNONc,NTOTc,NDEPVARc
c      is using in simplex, see comments above

c     due to fortran restrictions common values can't be parameters, so instead of
c     COMMON /simplex_param/ NLIN,NNON,NTOT,NDEPVAR
c     we must do that:
      NLINc=NLIN
      NNONc=NNON
      NTOTc=NTOT
      NDEPVARc=NDEPVAR
      

            NAMEPARAM="nameparam.txt"
            NAMESIMPLEX="namesimplex.txt"

            NF=0
c--------------------------------------------------
      IF (NF.NE.0) THEN
        print *, '  Please enter their numbers: '
        read *,(NUM(I),I=1,NF)

        DO I=1,NF
          IA(NUM(I))=0
        ENDDO
C
      ENDIF

c--------------------------------------------------
c INITIAL VALUES FOR MODEL PARAMETERS
	A(1:NTOT)=1.
        !TODO
        do i = nlin+1, ntot
           A(i)=1.d0*(i-5.)
        end do
        print *, a
C--------------------------------------------------
C
         NVNP=0  !  NVNP IS A COUNTER OF VARIABLE NONLINEAR PARAMETERS

         KKIK=1            !   KKIK IS A COUNTER OF FIXED LINEAR PARAMETERS
         DO 1616 I=1,NTOT
          IF (I.GT.NLIN) THEN
            IF (IA(I).EQ.1) NVNP=NVNP+1
            IF (IA(I).EQ.0) THEN
             COMMENT='   (FIXED)'
             ELSE
             COMMENT='   (VARIABLE)'
             ENDIF
          ELSE
            IF (I.EQ.NUM(KKIK)) KKIK=KKIK+1
          ENDIF
 1616    CONTINUE

c
          NVNP1=NVNP+1  !   THIS IS FOR SIMPLEX
C
c   initial values have been specified
C****** !!! NUMBER OF ITERATIONS IS NOW FROZEN !!! ***************
c      print *,'   How many ITERATIONS do you want to make ?'
c      read *, NITER
c*****************************************************************
C
C      NITER=70
       NITER=1

c      PRINT *,
c     *'  ENTER THE TOLERANCE FACTOR FOR SVD (BETWEEN 1.E-6 AND 1.E-15)'
c      READ *, TO

      TOL=1.d-6
C
      ICOMPAR=1
C
C     NOW, READ THE EXPERIMENTAL DATA SET.   IN THE FUTURE, THIS PART OF THE CODE (BETWEEN ASTERISKS) SHOULD BE
C     FURTHER DEVELOPED, SO THAT THE ARRAY B, INSTEAD OF (PDYN,DST,BYIMF,BZIMF)  WOULD CONTAIN
C     MORE ADVANCED PARAMETERS, E.G.,  (1) LINEAR FILTERS OF V-IMF-DEPENDENT REGRESSORS AND/OR DST-INDEX
C                                      (2) A COMBINATION OF DST AND ITS TIME DERIVATIVE
C                                      (3) ASY-INDEX
C
C
      OPEN (UNIT=1,FILE=DATANAME,STATUS='OLD')
C-------------------------------------------------------- ******************

      read(1, *) subset_len
      print *, subset_len

      allocate(subset(subset_len))

      call recalc_08 (2008,47,04,45,0,-400.d0,0.d0,0.d0)

      L=1
   33 READ (1, *,END=34,ERR=27) subset
      do subset_i = 1, subset_len

 777  FORMAT(I5,I4,2I3,3F8.2,3F9.1,F6.1,F6.1,f8.2,6F7.2)

C    CONVERT X,Y,Z,BX,BY,BZ TO SM !!!

          VX=-400.D0
          VY=0.D0
          VZ=0.D0
!          CALL RECALC_08(IYEAR,IDAY,IHOUR,MIN,0,VX,VY,VZ)
          PS=PSI

          cp = subset(subset_i)

          K=1
C
          B(K:K+2,L)=cp%r
          K=K+3
          B(K:K+2,L)=cp%n
          K=K+3

          r = sqrt(cp%r(1)**2+cp%r(2)**2+cp%r(3)**2)

!          m_dip = (/1000., 0., 0./)

          !TODO:
!          dip_b = 10.*dipole(m_dip, cp%r)
          call DIP_08 (cp%r(1),cp%r(2),cp%r(3),BXGSW,BYGSW,BZGSW)
          dip_b = (/BXGSW,BYGSW,BZGSW/)

          B(K, L) = 0.
          ! scalar mult
          do s_i = 1, 3
             B(K, L) = B(K, L) + cp%n(s_i)*dip_b(s_i)
          end do
!          print *, 'b_d', B(K, L)

          K=K+1

C--------------------------------------------------------*******************
          L=L+1
       end do

        GOTO 33

 27     print *,'  Error: L=',L
c        PAUSE
   34   L=L-1
        CLOSE(1)

        NPOINTS =L
C
      PRINT *,' '
      PRINT *,'                  NUMBER OF POINTS IN THE DATASET: L=',L
      PRINT *,' '
c
C    NOW CALCULATE THE AVERAGE (UNWEIGHTED) FIELD IN THE DATASET
C
      SUM=0.
      DO 35 I=1,L
  35  SUM=SUM+B(INDEPVAR+1,I)**2+B(INDEPVAR+2,I)**2+B(INDEPVAR+3,I)**2
C
      SUM=SQRT(SUM/NPOINTS)
C
      PRINT *,' '
      PRINT *,
     * '  AVERAGE UNWEIGHTED EXTERNAL B MAGNITUDE OVER THE DATASET:',SUM
      PRINT *,' '
c
C
       q_old=SUM
       B_rms=SUM
C

       IF (NNON.EQ.0) THEN
          Q=FUNK(XXXXX,NPOINTS,NTOT,NLIN,A,IA,INDEPVAR,NDEPVAR,B,
     *    WEIGHT,MODEL_DIP_SH)
          GOTO 12345
       ENDIF

       print *, 'a>', a

       CALL SIMPLEX (NPOINTS,NTOT,NLIN,A,IA,INDEPVAR,NDEPVAR,B,WEIGHT,
     _ NITER,MODEL_DIP_SH,P,NVNP1,NVNP)

C
C---------------------------------------------------

12345  CONTINUE
       print *, 'a>', q, a

      IF (ICOMPAR.EQ.1) THEN

       OPEN(UNIT=2,FILE='compar.dat')

        DO I=1,NPOINTS

          DO IN=1,INDEPVAR
          XI(IN)=B(IN,I)
          ENDDO

          CALL MODEL_DIP_SH (ID,A,XI,F,DER,IA,NTOT,NLIN,NNON,INDEPVAR,
     *      NDEPVAR,1)
          BBX(I)=B(INDEPVAR+1,I)
          BBY(I)=B(INDEPVAR+2,I)
          BBZ(I)=B(INDEPVAR+3,I)
          BBXM(I)=F(1)
!          BBYM(I)=F(2)
!          BBZM(I)=F(3)
!          print *, B(INDEPVAR+1:INDEPVAR+3,I)
!          print *, f(1)
C
          WRITE (2,202) (XI(M),M=1,3),
     *  F(1),B(INDEPVAR+1,I)!,F(2),B(INDEPVAR+2,I),F(3),B(INDEPVAR+3,I)
       end do

       CLOSE(2)
      ENDIF    !   END OF CASE ICOMPAR=1
 202  FORMAT(9F10.5)
C---------------------------------------------------
C
c  Now, evaluate correlation coefficients between observed and model
c     components of the fitted vector
c
      AVX=0.
      AVY=0.
      AVZ=0.
      AVXM=0.
      AVYM=0.
      AVZM=0.

      DO I=1,NPOINTS
      AVX=AVX+BBX(I)
!      AVY=AVY+BBY(I)
!      AVZ=AVZ+BBZ(I)
      AVXM=AVXM+BBXM(I)
!      AVYM=AVYM+BBYM(I)
! 1917 AVZM=AVZM+BBZM(I)
      end do
      AVX=AVX/NPOINTS
!      AVY=AVY/NPOINTS
!      AVZ=AVZ/NPOINTS
      AVXM=AVXM/NPOINTS
!      AVYM=AVYM/NPOINTS
!      AVZM=AVZM/NPOINTS

      RX=0.
      RY=0.
      RZ=0.
      SX=0.
      SY=0.
      SZ=0.
      SXM=0.
      SYM=0.
      SZM=0.
      DO I=1,NPOINTS
      RX=RX+(BBX(I)-AVX)*(BBXM(I)-AVXM)
!      RY=RY+(BBY(I)-AVY)*(BBYM(I)-AVYM)
!      RZ=RZ+(BBZ(I)-AVZ)*(BBZM(I)-AVZM)
      SX=SX+(BBX(I)-AVX)**2
!      SY=SY+(BBY(I)-AVY)**2
!      SZ=SZ+(BBZ(I)-AVZ)**2
      SXM=SXM+(BBXM(I)-AVXM)**2
!      SYM=SYM+(BBYM(I)-AVYM)**2
! 1918 SZM=SZM+(BBZM(I)-AVZM)**2
      end do
      RX=RX/SQRT(SX*SXM)
!      RY=RY/SQRT(SY*SYM)
!      RZ=RZ/SQRT(SZ*SZM)
      PRINT *, '  CORRELATIONS:'
      PRINT 1920, RX!,RY,RZ
      OPEN(UNIT=1,FILE='correlat_pos.par')
      WRITE (1,1920) RX!,RY,RZ
 1920 FORMAT('    RX=',F7.3,'   RY=',F7.3,'   RZ=',F7.3)
      CLOSE(1)

C----------------------------------------------------------------------
C----- Re-writing 'Bestfar_so_far.par' with the appropriate output name

      NAMEPAR="bla.txt"
c      print 111, 1.0d0
c            open (unit=1,file=NAMEPAR)
c            read (1,111) A
c            read (1,112) Q
c            read (1,113) B_rms
111         format(G15.6)
112         format('     Q=',G15.6)
113         format(' B_rms=',G15.6)
c            close(1)
C      A=1.0d0
C      Q=1.0d0
C      B_rms=1.0d0

c------ Writing new output file ------------------------

c      OUTFNAME="model_par.dat"
            open (unit=1,file=OUTFNAME,status='UNKNOWN')
            write (1,111) A
            write (1,112) Q
            write (1,113) B_rms
            close(1)

C----------------------------------------------------------------------
      STOP
      END

