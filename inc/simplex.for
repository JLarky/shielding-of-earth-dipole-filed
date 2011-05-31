c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C(((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((((
C
      module Msimplex
       interface
         double precision function srms
     _     (N, A, XI, MX, MY, B, WEIGHT, F, D, MODEL, IA, M)
           IMPLICIT REAL*8 (A-H,O-Z)
           COMMON /simplex_param/ NLIN,NNON,NTOT,NDEPVAR
           integer, intent(in) :: IA(M)
           real*8 :: A(NTOT),XI(MX),B(MX+MY,*),WEIGHT(*),
     _     F(*), D(NDEPVAR,*)
           external model
         end function
       end interface
      end module Msimplex


         double precision FUNCTION SRMS
     *    (N, A, XI, MX, MY, B, WEIGHT, F, D, MODEL, IA, M)      !    ADDED XI IN THE LIST OF PARAMETERS,
C                                                                        TO AUTOMATICALLY CHANGE ITS LENGTH WHEN NEEDED
C      Calculates standard RMS deviation of observed
C  dependent variables from the computed ones.
C
      COMMON /simplex_param/ NLIN,NNON,NTOT,NDEPVAR
      integer, intent(in) :: IA(M)
      real*8 :: F(*),D(NDEPVAR,*),A(*),XI(MX),
     * B(MX+MY,*),WEIGHT(*)
C
      SF = 0.0D0
      DO K=1,N
         XI(1:MX)=B(1:MX,K)
         CALL MODEL (K, A, XI, F, D, IA, M, NLIN, NNON, MX, MY, 1)
         DO L=1,MY
            SF=SF+(B(MX+L,K)-F(L))**2
         END DO
      END DO
      SRMS = SF/N
      END FUNCTION SRMS
C
C================================================================================================================
C
      DOUBLE PRECISION FUNCTION FUNK (X,N,M,ML,A,IA,MX,MY,B,WEIGHT,
     *  MODEL)  !   ADDED:  NNON, WEIGHT
C
C    X - VECTOR OF VARIABLE NONLINEAR PARAMETERS
C    N - # OF POINTS IN THE DATA SET
C    M - TOTAL NUMBER OF PARAMETERS (LINEAR+NONLINEAR, BOTH VARIABLE AND FIXED)
C    ML - NUMBER OF LINEAR PARAMETERS (BOTH VARIABLE AND FIXED)
C    NNON - NUMBER OF NONLINEAR PARAMETERS (VARIABLE ONLY)
C    A - VECTOR OF ALL MODEL PARAMETERS (LINEAR+NONLINEAR, BOTH VARIABLE AND FIXED)
C    IA - FLAG VECTOR, IDENTIFYING WHICH PARAMETERS ARE VARIABLE (IA(i)=1) AND WHICH ARE FIXED (IA(k)=1)
C    MX - # OF INDEPENDENT MODEL INPUT VARIABLES (IN THIS CASE, MX=4: X,Y,Z,PS)
C    MY - # OF DEPENDENT MODEL VARIABLES (MY=3, CORRESPONDING TO 3 COMPONENTS OF B VECTOR)
C    B - ARRAY CONTAINING N EXPERIMENTAL RECORDS, EACH RECORD INCLUDES MEASURED
C                          INPUT VARIABLES (X,Y,Z,PS) AND 3 B COMPONENTS
C    WEIGHT - ARRAY CONTAINING PRE-DETERMINED WEIGHTS FOR EACH DATA POINT
C    MODEL - NAME OF A MODEL FIELD SUBROUTINE, TO BE USED FOR CALCULATING THE MODEL
C    FIELD COMPONENTS FOR EACH (X,Y,Z,PS) IN THE EXPERIMENTAL DATASET

      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL MODEL

      COMMON /simplex_param/ NLIN,NNON,NTOT,NDEPVAR
c      PARAMETER (NLIN=13,NNON=2,NNON1=3,NTOT=15)
c      PARAMETER (NTOT=15)
c      PARAMETER (NNON=2)
c      PARAMETER (NUMVTOT=11)       ! NUMVTOT=INDEPVAR+NDEPVAR (adjust in each individual case)
c      PARAMETER (NPNT_UPP=1000000)  ! UPPER LIMIT ON THE NUMBER OF DATA POINTS
c      PARAMETER (N_ELEM_B=11000000) ! =NUMVTOT*NPNT_UPP  -
c      PARAMETER (NDEPVAR=3)
C
      DIMENSION B(*),WEIGHT(*)
      DIMENSION A(M),IA(M),X(*),F(NDEPVAR),D(NDEPVAR,NTOT),RH(150),
     * FM(150,150),      ! THE SIZE OF ARRAYS ASSUMED FOR THE
     * BB(150),XI(MX),U(150,150),V(150,150),W(150),XX(150)              !  MAXIMUM OF 150 LINEAR PARAMETERS; NO NEED TO CHANGE
C                                                               ANYTHING, AS LONG AS THE ACTUAL NUMBER IS LOWER THAN THAT

!      print *, 't', 4
      NVP=0
      DO 111 I=1,NNON
      IF (IA(ML+I).EQ.0) GOTO 111
      NVP=NVP+1
      A(ML+I)=X(NVP)
 111  CONTINUE
C
      IF (ML.EQ.0) THEN
      FUNK=SQRT(SRMS(N,A,XI,MX,MY,B,WEIGHT,F,D,MODEL,IA,M))  !  NOW SRMS INCLUDES THE NORMALIZATION FACTOR
      RETURN                                                         !    SO I REMOVED THE DIVISION OF SRMS BY N
      ENDIF
C
      NL=0
      DO I=1,ML    !   CALCULATE THE NUMBER NL OF VARIABLE LINEAR PARAMETERS
      IF (IA(I).EQ.1) NL=NL+1
      END DO
C
      ME=MX+MY


      D(1:MY,1:NTOT)=0.d0
      RH(1:ML)=0.d0
      FM(1:ML,1:ML)=0.d0
C 
        DO  16  K = 1, N    !   SUMMATION OVER ALL N DATA POINTS
	   LS = ME * (K - 1)

	   DO  9  L = 1, MX
    9        XI(L) = B(LS+L) !   FILL THE ARRAY XI WITH ARGUMENTS FOR THE CURRENT DATA POINT
         LS = LS + MX        !   PREPARE LS FOR THE NEXT DATA POINT

!         print *, 't', 5, a(NTOT-1)
           CALL MODEL (K,A,XI,F,D,IA,M,ML,NNON,MX,MY,1)
             NM = 0

             II=0

           DO  15  I = 1, ML

c              PRINT *, I,M,ML

             IF (IA(I).EQ.0)  GOTO  15
             NM = NM + 1                 !    NUMBER OF CURRENT ROW IN THE MATRIX TO BE INVERTED
             II = II + 1

	     LPI = (I - 1) * MY
c
c  Calculation of right-hand-side terms in the linear system of fitting eqs:
c
	     DO  10  L = 1, MY
   10          RH(NM) = RH(NM) + B(LS+L) * D(L,I)
c

             JJ=0

             DO  14  J = 1, ML
	       LPJ = (J - 1) * MY

               IF (IA(J).EQ.1)  GOTO  12
c  in the case of a fixed parameter, move the term to the right-hand side:
               DO  11  L = 1, MY
                 RH(NM) = RH(NM) - D(L,I) * D(L,J) * A(J)
   11            CONTINUE

           GOTO  14

   12          CONTINUE

            JJ=JJ+1

c  Calculation of left-hand-side terms in the linear system of fitting eqs:

             DO  13  L = 1, MY
   13            FM(II,JJ) = FM(II,JJ) + D(L,I) * D(L,J)

   14        CONTINUE

   15      CONTINUE
   16    CONTINUE
C
C      Solve the normal system of simultaneous linear
C  equations with symmetric coefficient matrix:
C
        DO 25 I=1,ML
  25    BB(I)=RH(I)

      print *, 'test5'
         CALL SINGVALDE(FM,RH,XX,U,V,W,BB,NL)
C
         MI = 0
         DO  18  I = 1, ML
           IF (IA(I).EQ.0)  GOTO  18
           MI = MI + 1
           A(I) = XX (MI)
   18    CONTINUE
C
c      SSRMS=SRMS(N,A,XI,MX,MY,B,WEIGHT,F,D,MODEL,IA,M)

         FUNK=SQRT(SRMS(N,A,XI,MX,MY,B,WEIGHT,F,D,MODEL,IA,M))    !  NOW SRMS INCLUDES THE NORMALIZATION FACTOR
C                                                                         !    SO I REMOVED THE DIVISION OF SRMS BY N

         RETURN
         END
C
coooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

      SUBROUTINE SINGVALDE(FM,RH,X,U,V,W,BB,NL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FM(150,150),RH(150),W(NL),U(NL,NL),V(NL,NL),X(150),
     * BB(150)
      COMMON /TOLERANCE/ TOL

      DO I=1,NL              !   NUMERICAL RECIPES, P.57
        DO J=1,NL
          U(I,J)=FM(I,J)
        ENDDO
      ENDDO

      CALL SVDCMP(U,NL,NL,NL,NL,W,V)

      WMAX=0.D0

      DO J=1,NL
        IF (W(J).GT.WMAX) WMAX=W(J)
      ENDDO

      WMIN=WMAX*TOL  !    EXPERIMENT WITH THE TOLERANCE CONSTANT

      DO J=1,NL
        IF (W(J).LT.WMIN) W(J)=0.D0
      ENDDO

      CALL SVBKSB (U,W,V,NL,NL,NL,NL,BB,X)

      RETURN
      END

coooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
C
C =========================================================
      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension a(mp,np),v(np,np),w(np)
      PARAMETER (NMAX=500)
CU    USES pythag
      dimension rv1(NMAX)

      g=0.0d0
      scale=0.0d0
      anorm=0.0d0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if(i.le.m)then
          do 11 k=i,m
            scale=scale+dabs(a(k,i))
11        continue
          if(scale.ne.0.0d0)then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale
              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-dsign(dsqrt(s),f)
            h=f*g-s
            a(i,i)=f-g
            do 15 j=l,n
              s=0.0d0
              do 13 k=i,m
                s=s+a(k,i)*a(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
14            continue
15          continue
            do 16 k=i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n
            scale=scale+dabs(a(i,k))
17        continue
          if(scale.ne.0.0d0)then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-dsign(dsqrt(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0d0
              do 21 k=l,n
                s=s+a(j,k)*a(i,k)
21            continue
              do 22 k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n
              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=max(anorm,(dabs(w(i))+dabs(rv1(i))))
25    continue

      do 32 i=n,1,-1
        if(i.lt.n)then
          if (g.ne.0.0d0) then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0d0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue
29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0d0
            v(j,i)=0.0d0
31        continue
        endif
        v(i,i)=1.0d0
        g=rv1(i)
        l=i
32    continue

      do 39 i=min(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          a(i,j)=0.0d0
33      continue
        if(g.ne.0.0d0)then
          g=1.0d0/g
          do 36 j=l,n
            s=0.0d0
            do 34 k=l,m
              s=s+a(k,i)*a(k,j)
34          continue
            f=(s/a(i,i))*g
            do 35 k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
35          continue
36        continue
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0d0
38        continue
        endif
        a(i,i)=a(i,i)+1.0d0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((dabs(rv1(l))+anorm).eq.anorm)  goto 2
            if((dabs(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0d0
          s=1.0d0
          do 43 i=l,k
            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((dabs(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0d0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w(k)
          if (l.eq.k) then
            if(z.lt.0.0d0)then
              w(k)=-z
              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            goto 3
          endif
c          if(its.eq.30) pause 'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
          g=pythag(f,1.0d0)
          f=((x-z)*(x+z)+h*((y/(f+dsign(g,f)))-h))/x
          c=1.0d0
          s=1.0d0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)
            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
45          continue
            z=pythag(f,h)
            w(j)=z
            if(z.ne.0.0d0)then
              z=1.0d0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)
            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0d0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .)+1YX39'=K+1.

C ============================================================


      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension b(mp),u(mp,np),v(np,np),w(np),x(np)
      PARAMETER (NMAX=500)
      dimension tmp(NMAX)
      do 12 j=1,n
        s=0.d0
        if(w(j).ne.0.d0)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.d0
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software .)+1YX39'=K+1.
C ============================================================
      DOUBLE PRECISION FUNCTION pythag(a,b)
      IMPLICIT REAL*8 (A-H,O-Z)

      absa=dabs(a)
      absb=dabs(b)
      if(absa.gt.absb) then
        pythag=absa*dsqrt(1.d0+(absb/absa)**2)
      else
        if(absb.eq.0.d0)then
          pythag=0.d0
        else
          pythag=absb*dsqrt(1.d0+(absa/absb)**2)
        endif
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ?421.1.9.
C ===========================================================
C
      SUBROUTINE SIMPLEX(N,M,ML,A,IA,MX,MY,B,WEIGHT,NI,MODEL,P,MP,NP)
C
C     SEE COMMENTS FOR THE FUNCTION FUNK ABOVE FOR A DETAILED EXPLANATION
C       OF THE PARAMETERS FROM N TO MODEL.
C
C  THE INTEGER PARAMETER NP IS THE NUMBER OF VARIABLE NONLINEAR PARAMETERS
C     (NP<=M-ML) AND MP=NP+1
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
      EXTERNAL MODEL,FUNK
      CHARACTER*80 NAMEPARAM, NAMESIMPLEX
C
c      PARAMETER (NPNT_UPP=1000000)  ! UPPER LIMIT ON THE NUMBER OF DATA POINTS
c      PARAMETER (N_ELEM_B=8000000) ! =NUMVTOT*NPNT_UPP  -

c      DIMENSION B(N_ELEM_B),WEIGHT(NPNT_UPP)
      DIMENSION B(*),WEIGHT(*)
      DIMENSION A(M),IA(M),X(NP),Y(MP),ANEL(NP),P(MP,NP)
      COMMON /NAMES/ NAMEPARAM,NAMESIMPLEX
C
      cur_Q=.0d0
      L = 0

      PRINT *,
     *  ' PARAMETERS WILL BE CALCULATED BY SIMPLEX METHOD;'
C
C  INITIALIZE THE SIMPLEX BY SPECIFYING NVP-COMPONENT VECTOR OF VARIABLE
C     NONLINEAR PARAMETERS IN (NVP+1) POINTS OF NVP-DIMENSIONAL SPACE
C     HERE NVP IS LESS THAN TOTAL NUMBER OF NONLINEAR PARAMETERS,
C     SINCE SOME OF THEM ARE FIXED.
C  FIRST OF ALL, FIND THE NUMBER OF VARIABLE PARAMETERS, N, AND PUT THEIR
C     INITIAL VALUES IN THE ARRAY ANEL:
C
      NNON=M-ML
 333  prev_Q=cur_Q
      NVP=0

      do I=1,NNON
         if (IA(ML+I).NE.0) then ! non linear parameters
            NVP=NVP+1
            ANEL(NVP)=A(ML+I)
         end if
      end do

      NVP1=NVP+1
      
      PRINT *, NVP,'  NONLINEAR PARAMETERS WILL BE VARIED'

C  NOW, INITIALIZE THE VECTORS X AND Y FOR THE AMOEBA:
C
c***** !!! SIMPLEX OPTION IS NOW FROZEN !!! **************************
c      PRINT *, '  INITIALIZE THE SIMPLEX FROM SCRATCH (1) '
c      PRINT *, '       OR READ THE Y-VECTOR FROM DISK (2)  ?'
c      READ *, IREAD
c*********************************************************************
      IREAD=1

      if (IREAD.EQ.1) THEN
         do I=1,NVP1
            P(I,1:NVP)=ANEL(1:NVP)
            IF (ANEL(I).NE.0.D0) THEN
               P(I,I)=P(I,I)+0.15*ANEL(I)
            ELSE
               P(I,I)=P(I,I)+0.15
            ENDIF
            X(1:NVP)=P(I,1:NVP)
            Y(I)=FUNK(X,N,M,ML,A,IA,MX,MY,B,WEIGHT,MODEL)
            cur_Q=Y(I)
            PRINT *, '  INITIALIZING SIMPLEX: I=',I,' OF',NVP1,
     _           '  Q=',Y(I)
         end do
      else
         OPEN (UNIT=1,FILE=NAMESIMPLEX)
         READ (1,4444) ((P(I,K),K=1,NVP),Y(I),I=1,NVP1)
 4444    FORMAT(G17.8)
         CLOSE (1)
      end if

      FTOL=1.d-10
C
      CALL AMOEBA(P,Y,NVP+1,NVP,NVP,FTOL,FUNK,N,M,ML,A,IA,MX,MY,B,
     *     WEIGHT,NI,MODEL)

C****** !!! OPTION OF ANOTHER S-RUN IS NOW FROZEN !!! *********************
c       PRINT *, '  ANOTHER SIMPLEX RUN  ? (1-YES, 0-NO)'
c       READ *, IANOT
c**************************************************************************

       err_Q=abs(prev_Q-cur_Q)/abs(cur_Q)
       L = L +1
       print *, 'Error Q=', err_Q
!       if ((err_Q.lt.1d-7).or.(L.gt.130)) then
       if ((L.gt.60)) then
          print *, 'Iterations= ',L
          IANOT=0
       else
          IANOT=1
       endif

       IF (IANOT.EQ.1) GOTO 333
C
      RETURN
      END
C#####################################################################
c
c
      SUBROUTINE amoeba(p,y,mp,np,ndim,ftol,funk,
     *  NN,MM,ML,A,IA,MX,MY,B,WEIGHT,NI,MODEL)

      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*80 NAMEPARAM, NAMESIMPLEX

c      PARAMETER (NPNT_UPP=1000000)  ! UPPER LIMIT ON THE NUMBER OF DATA POINTS
c      PARAMETER (N_ELEM_B=8000000) ! =NUMVTOT*NPNT_UPP  -

c      DIMENSION B(N_ELEM_B),WEIGHT(NPNT_UPP)
      DIMENSION B(*),WEIGHT(*)

      dimension p(mp,np),y(mp), A(MM),IA(MM)
      PARAMETER (NMAX=20)
      dimension psum(nmax)
      common /q/ Q,Q_old,B_rms

      COMMON /NAMES/ NAMEPARAM,NAMESIMPLEX

      EXTERNAL funk, MODEL
c
C            USES: amotry,funk, MODEL
c
      iter=0
1     do 12 n=1,ndim
        sum=0.
        do 11 m=1,ndim+1
          sum=sum+p(m,n)
11      continue
        psum(n)=sum
12    continue
2     ilo=1
      if (y(1).gt.y(2)) then
        ihi=1
        inhi=2
      else
        ihi=2
        inhi=1
      endif
      do 13 i=1,ndim+1
        if(y(i).le.y(ilo)) ilo=i
        if(y(i).gt.y(ihi)) then
          inhi=ihi
          ihi=i
        else if(y(i).gt.y(inhi)) then
          if(i.ne.ihi) inhi=i
        endif
13    continue
      rtol=2.d0*dabs(y(ihi)-y(ilo))/(dabs(y(ihi))+dabs(y(ilo)))
      if (rtol.lt.ftol) then
        swap=y(1)
        y(1)=y(ilo)
        y(ilo)=swap
        do 14 n=1,ndim
          swap=p(1,n)
          p(1,n)=p(ilo,n)
          p(ilo,n)=swap
14      continue
        return
      endif
c
         PRINT *,
     *                                 '            IT=',iter
C
      if (iter.ge.NI) THEN
            PRINT *, ' ITER>ITMAX'
            RETURN
      endif
c
      iter=iter+2
      ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,-1.0d0,
     *  NN,MM,ML,A,IA,MX,MY,B,WEIGHT,MODEL)
      if (ytry.le.y(ilo)) then
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,2.0d0,
     *  NN,MM,ML,A,IA,MX,MY,B,WEIGHT,MODEL)
      else if (ytry.ge.y(inhi)) then
        ysave=y(ihi)
        ytry=amotry(p,y,psum,mp,np,ndim,funk,ihi,0.5d0,
     *  NN,MM,ML,A,IA,MX,MY,B,WEIGHT,MODEL)
        if (ytry.ge.ysave) then
          do 16 i=1,ndim+1
            if(i.ne.ilo)then
              do 15 j=1,ndim
                psum(j)=0.5d0*(p(i,j)+p(ilo,j))
                p(i,j)=psum(j)
15            continue
              y(i)=funk(psum, NN,MM,ML,A,IA,MX,MY,B,WEIGHT,MODEL)
C
       Q=Y(I)
       if (Q.lt.Q_old) then
            PRINT *, '   Q=',Y(I)  ! THIS IS A TEST PRINT, TO KEEP TRACK OF Q
            Q_old=Q
            open (unit=1,file=NAMEPARAM)
            write (1,111) A
            write (1,112) Q
            write (1,113) B_rms
111         format(G15.6)
112         format('     Q=',G15.6)
113         format(' B_rms=',G15.6)
            close(1)

          OPEN (UNIT=1,FILE=NAMESIMPLEX)
          WRITE (1,4444) ((P(II,KK),KK=1,NP),Y(II),II=1,MP),Q
4444      FORMAT(G17.8)
          CLOSE (1)

       endif
C
            endif
16        continue
          iter=iter+ndim
          goto 1
        endif
      else
        iter=iter-1
      endif
      goto 2
      END
C  (C) Copr. 1986-92 Numerical Recipes Software v1]r1wj,W.
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      FUNCTION amotry(p,y,psum,mp,np,ndim,funk,ihi,fac,
     * N,M,ML,A,IA,MX,MY,B,WEIGHT,MODEL)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*80 NAMEPARAM,NAMESIMPLEX
      COMMON /NAMES/ NAMEPARAM,NAMESIMPLEX

c      PARAMETER (NPNT_UPP=1000000)  ! UPPER LIMIT ON THE NUMBER OF DATA POINTS
c      PARAMETER (N_ELEM_B=8000000) ! =NUMVTOT*NPNT_UPP  -

c      DIMENSION B(N_ELEM_B),WEIGHT(NPNT_UPP)
      DIMENSION B(*),WEIGHT(*)
      dimension p(mp,np),psum(np),y(mp), A(M),IA(M)
      PARAMETER (NMAX=25)
      dimension ptry(nmax)

      common /q/ Q,Q_old,B_rms

      EXTERNAL funk, MODEL
c
C            USES: funk , MODEL
c
      fac1=(1.d0-fac)/ndim
      fac2=fac1-fac
      do 11 j=1,ndim
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
11    continue
      ytry=funk(ptry,N,M,ML,A,IA,MX,MY,B,WEIGHT,MODEL)
C

       Q=ytry
       if (Q.lt.Q_old) then
            PRINT *, ' Amotry:  Q=',YTRY  ! THIS IS A TEST PRINT, TO KEEP TRACK OF Q
!         print *, 't', 3
            Q_old=Q
            open (unit=1,file=NAMEPARAM)
            write (1,111) A
            write (1,112) Q
            write (1,113) B_rms
111         format(G15.6)
112         format('     Q=',G15.6)
113         format(' B_rms=',G15.6)
            close(1)

         OPEN (UNIT=1,FILE=NAMESIMPLEX)
         WRITE (1,4444) ((P(II,KK),KK=1,NP),Y(II),II=1,MP),Q
4444     FORMAT(G17.8)
         CLOSE (1)

       endif
C
      if (ytry.lt.y(ihi)) then
        y(ihi)=ytry
        do 12 j=1,ndim
          psum(j)=psum(j)-p(ihi,j)+ptry(j)
          p(ihi,j)=ptry(j)
12      continue
      endif
      amotry=ytry
      return
      END


