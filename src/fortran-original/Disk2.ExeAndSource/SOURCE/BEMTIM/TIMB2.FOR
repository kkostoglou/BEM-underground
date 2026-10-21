C------------------------------------------------
        SUBROUTINE KERNEL(NCOUNT,Y1,Y2,Y3,X1,X2,X3,RN,ZU,ZT,N9)
C
C       SUPPLIES VALUES OF SINGULAR KERNELS-STOKE'S SOLUTION FROM
C       ERINGEN-LINEAR TIME VARIATION
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RN(3),ZU(2,2),ZT(2,2),SU(3,3),ST(3,3)
        DIMENSION H(4),HB(4)
        COMMON/CONTR/NNODE,NELEM,NP
        COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM
C
        NTINT=N9
        CC1=C1**2
        CC2=C2**2
C
        CALL CLRARR(ZU,4)
        CALL CLRARR(ZT,4)
        DT=TIM/DFLOAT(NTINT)
        TT=DT*DFLOAT(NCOUNT)
        R=DSQRT((X1-Y1)**2+(X2-Y2)**2+(X3-Y3)**2)
        RR=R*R
        IF (TT .LT. (R/C1)) RETURN
C
C       ONE REGULAR TIME INTERVAL PLUS ONE BACKSTEP
C
        DO 7 L=1,2
        TT=DT*DFLOAT(NCOUNT+1-L)
        TTB=DT*DFLOAT(NCOUNT-L)
        CALL CLRARR(H,4)
        CALL CLRARR(HB,4)
        IF (TT .LE. 1.D-15) GO TO 14
        IF (TT .LT. (R/C1)) GO TO 14
        H(1)=1.0D0
        IF (TT .GE. (R/C2)) H(2)=1.0D0
        H(3)=1.0D0/(2.0D0*CC1)
        H(4)=1.0D0/(2.0D0*CC2)
        IF (TT .LT. (R/C2)) H(4)=0.5D0*(TT/R)**2
        C12=C2
        IF (TT .LT. (R/C2)) C12=R/TT
        RRCAP1=(TT*C1)**2
        ZCAP1=DSQRT(RRCAP1-RR)
        RRCAP2=(TT*C2)**2
        ZCAP2=0.0D0
        IF (RRCAP2 .GT. RR) ZCAP2=DSQRT(RRCAP2-RR)
        RRCAP9=(TT*C12)**2
        ZCAP12=0.0D0
        IF (RRCAP9 .GT. RR) ZCAP12=DSQRT(RRCAP9-RR)
        ABC=0.0D0
        CALL KERNA(Y1,Y2,ABC,X1,X2,RN,ZCAP1,ZCAP2,ZCAP12,SU,ST,H,L,N9)
        DO 8 I=1,2
        DO 8 J=1,2
          ZU(I,J)=ZU(I,J)+SU(I,J)
          ZT(I,J)=ZT(I,J)+ST(I,J)
  8     CONTINUE
  14    CONTINUE
        IF (TTB .LE. 1.D-15) GO TO 15
        IF (TTB .LT. (R/C1)) GO TO 15
        HB(1)=1.0D0
        IF (TTB .GE. (R/C2)) HB(2)=1.0D0
        HB(3)=1.0D0/(2.0D0*CC1)
        HB(4)=1.0D0/(2.0D0*CC2)
        IF (TTB .LT. (R/C2)) HB(4)=0.5D0*(TTB/R)**2
        C12=C2
        IF (TTB .LT. (R/C2)) C12=R/TTB
        RRCAP1=(TTB*C1)**2
        ZCAP1=DSQRT(RRCAP1-RR)
        RRCAP2=(TTB*C2)**2
        ZCAP2=0.0D0
        IF (RRCAP2 .GT. RR) ZCAP2=DSQRT(RRCAP2-RR)
        RRCAP9=(TTB*C12)**2
        ZCAP12=0.0D0
        IF (RRCAP9 .GT. RR) ZCAP12=DSQRT(RRCAP9-RR)
        DO 6 I=1,4
        HB(I)=-HB(I)
  6     CONTINUE
        DEF=0.0D0
        CALL KERNA(Y1,Y2,DEF,X1,X2,RN,ZCAP1,ZCAP2,ZCAP12,SU,ST,HB,L,N9)
        DO 9 I=1,2
        DO 9 J=1,2
          ZU(I,J)=ZU(I,J)+SU(I,J)
          ZT(I,J)=ZT(I,J)+ST(I,J)
  9     CONTINUE
  15    CONTINUE
  7     CONTINUE
        RETURN
        END
C------------------------------------------------
        SUBROUTINE KERNA(Y1,Y2,Y3,X1,X2,RN,ZCAP1,ZCAP2,ZCAP12,
     &                   SU,ST,H,L,N9)
C
C       INTEGRATION OF 3-D KERNELS ALONG Z-DIRECTION
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION SU(3,3),ST(3,3),SS(3,3,3)
        DIMENSION AA(3,3),AB(3,3),AC(3,3)
        DIMENSION BA(3,3,3),BB(3,3,3),BC(3,3,3),BD(3,3,3),BE(3,3,3)
        DIMENSION RN(3),H(4),D(3,3),RC(3),S(8)
        COMMON/CONTR/NNODE,NELEM,NP
        COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM
        COMMON/GAUSS/GPN(10,10),WTN(10,10),GPS(10,10),WTS(10,10)
        DATA D/1.0D0,0.0D0,0.0D0,0.0D0,1.0D0,0.0D0,0.0D0,0.0D0,1.0D0/
C
        NTINT=N9
        PI=3.14159265359D0
        DT=TIM/DFLOAT(NTINT)
        CC1=C1**2
        CC2=C2**2
        CCC1=C1**3
        CCC2=C2**3
        CON=1.0D0/(4.0D0*PI*RHO)
        CAN=1.0D0/(4.0D0*PI)
C
        SHP=1.0D0/DT
        CAI=H(4)-H(3)
        CB=H(1)
        CC=H(2)
        DO 10 I=1,6
        S(I)=1.0D+00
  10    CONTINUE
        S(7)=SHP
        S(8)=SHP
        IF (L .EQ. 1) THEN
          S(7)=-SHP
          S(8)=-SHP
          GO TO 11
        ENDIF
        DO 12 I=1,6
        S(I)=0.0D+00
  12    CONTINUE
  11    CONTINUE
C
        CALL CLRARR(AA,9)
        CALL CLRARR(AB,9)
        CALL CLRARR(AC,9)
        CALL CLRARR(BA,27)
        CALL CLRARR(BB,27)
        CALL CLRARR(BC,27)
        CALL CLRARR(BD,27)
        CALL CLRARR(BE,27)
C
        NGP=8
        DO 1 IX=1,NGP
        EK=GPN(NGP,IX)
        WT=WTN(NGP,IX)
C
C       THE C1 WAVE CONTRIBUTIONS
C
        X3=0.5D0*ZCAP1*(1.0D0+EK)
        R=DSQRT((X1-Y1)**2+(X2-Y2)**2+(X3-Y3)**2)
        RR=R*R
        RC(1)=(X1-Y1)/R
        RC(2)=(X2-Y2)/R
        RC(3)=(X3-Y3)/R
        AMUL=WT*ZCAP1
        DO 25 I=1,3
        DO 25 J=1,3
        A2=RC(I)*RC(J)/R
        AB(I,J)=AB(I,J)+(A2/CC1)*AMUL
        DO 25 K=1,3
        B1=RC(I)*RC(J)*RC(K)/RR
        B2=(D(I,J)*RC(K)+D(I,K)*RC(J)+D(J,K)*RC(I))/RR
        B3=RC(K)*D(I,J)*(1.0D0-2.0D0*CC2/CC1)/RR
        B4=2.0D0*B1*R/C2
        BB(I,J,K)=BB(I,J,K)+(-2.0D0*(6.0D0*B1-B2)*(CC2/CC1)-B3)*AMUL
        BD(I,J,K)=BD(I,J,K)+(-B4*(CCC2/CCC1)-B3*(R/C1))*AMUL
  25    CONTINUE
C
C       THE C2 WAVE CONTRIBUTIONS
C
        X3=0.5D0*ZCAP2*(1.0D0+EK)
        R=DSQRT((X1-Y1)**2+(X2-Y2)**2+(X3-Y3)**2)
        RR=R*R
        RC(1)=(X1-Y1)/R
        RC(2)=(X2-Y2)/R
        RC(3)=(X3-Y3)/R
        AMUL=WT*ZCAP2
        DO 26 I=1,3
        DO 26 J=1,3
        A2=RC(I)*RC(J)/R
        A3=D(I,J)/R
        AC(I,J)=AC(I,J)+((-A2+A3)/CC2)*AMUL
        DO 26 K=1,3
        B1=RC(I)*RC(J)*RC(K)/RR
        B2=(D(I,J)*RC(K)+D(I,K)*RC(J)+D(J,K)*RC(I))/RR
        B4=2.0D0*B1*R/C2
        B5=(D(I,K)*RC(J)+D(J,K)*RC(I))/RR
        BC(I,J,K)=BC(I,J,K)+(2.0D0*(6.0D0*B1-B2)-B5)*AMUL
        BE(I,J,K)=BE(I,J,K)+(B4-B5*(R/C2))*AMUL
  26    CONTINUE
C
C       THE C1,C2 WAVE CONTRIBUTIONS
C
        X3=0.5D0*ZCAP12*(1.0D0+EK)
        R=DSQRT((X1-Y1)**2+(X2-Y2)**2+(X3-Y3)**2)
        RR=R*R
        RC(1)=(X1-Y1)/R
        RC(2)=(X2-Y2)/R
        RC(3)=(X3-Y3)/R
        AMUL=WT*ZCAP12
        DO 27 I=1,3
        DO 27 J=1,3
        A3=D(I,J)/R
        A2=RC(I)*RC(J)/R
        A1=3.0D0*A2-A3
        AA(I,J)=AA(I,J)+A1*AMUL
        DO 27 K=1,3
        B1=RC(I)*RC(J)*RC(K)/RR
        B2=(D(I,J)*RC(K)+D(I,K)*RC(J)+D(J,K)*RC(I))/RR
        BA(I,J,K)=BA(I,J,K)-6.0D0*CC2*(5.0D0*B1-B2)*AMUL
  27    CONTINUE
  1     CONTINUE
C
        DO 28 I=1,3
        DO 28 J=1,3
        SU(I,J)=CON*(AA(I,J)*S(1)*CAI+
     &               AB(I,J)*S(2)*CB+
     &               AC(I,J)*S(3)*CC)
        DO 28 K=1,3
        SS(I,J,K)=CAN*(BA(I,J,K)*S(4)*CAI+
     &                 BB(I,J,K)*S(5)*CB+
     &                 BC(I,J,K)*S(6)*CC+
     &                 BD(I,J,K)*S(7)*CB+
     &                 BE(I,J,K)*S(8)*CC)
  28    CONTINUE
C
        DO 29 I=1,3
        DO 29 K=1,3
        ST(I,K)=SS(I,1,K)*RN(1)+SS(I,2,K)*RN(2)+SS(I,3,K)*RN(3)
  29    CONTINUE
        RETURN
        END
C------------------------------------------------
        SUBROUTINE ZERNEL(Y1,Y2,X1,X2,RN,ZU,ZT)
C
C       2-D ELASTOSTATIC KERNELS
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION RN(2),ZU(2,2),ZT(2,2),D(2,2),RC(2)
        COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM
        DATA D/1.0D0,0.0D0,0.0D0,1.0D0/
C
        CALL CLRARR(ZU,4)
        CALL CLRARR(ZT,4)
C
        PI=3.14159265359D0
        R=DSQRT((X1-Y1)**2+(X2-Y2)**2)
        RC(1)=(X1-Y1)/R
        RC(2)=(X2-Y2)/R
        DRN=RC(1)*RN(1)+RC(2)*RN(2)
C
        DO 1 I=1,2
        DO 1 J=1,2
        ZU(I,J)=(-(3.0D0-4.0D0*RNU)*DLOG(R)*D(I,J)+RC(I)*RC(J))/
     &              (8.0D0*PI*RMU*(1.0D0-RNU))
        ZT(I,J)=-(DRN*((1.0D0-2.0D0*RNU)*D(I,J)+2.0D0*RC(I)*RC(J))-
     &    (1.0D0-2.0D0*RNU)*(RC(I)*RN(J)-RC(J)*RN(I)))/
     &    (4.0D0*PI*(1.0D0-RNU)*R)
  1     CONTINUE
        RETURN
        END
C---------------------------------------------------------------------
C
        SUBROUTINE GENR(NCOUNT,KONN,MBCU,ZGM,ZFM,ZLOAD,ZTEMP,ZCOL,
     &         ZB,ZPQ,N1,N2,N3,N4,SFDIAG,LBCU,NDTYPE,LINV,MINV)
C
C       READING OF LOADING AND DEFORMATION OF SYSTEM MATRICES
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION KONN(N3,N4),MBCU(N2)
        DIMENSION ZGM(N2,N2),ZFM(N2,N2),ZLOAD(N2),ZTEMP(N2)
        DIMENSION ZCOL(N2),ZB(N2),ZPQ(N2)
        DIMENSION ZKU(3,2,2),ZKT(3,2,2),DJUMP(2,2)
        DIMENSION SFDIAG(N1,2,2)
        DIMENSION LBCU(N3),NDTYPE(N1)
        DIMENSION LINV(N2),MINV(N2)
        COMMON/CONTR/NNODE,NELEM,NP
        COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM
        COMMON/FREE/NFFS,NFRE
C
        N5=N2*N2
        CALL CLRARR(ZGM,N5)
        CALL CLRARR(ZFM,N5)
        DATA DJUMP/1.0D0,0.0D0,0.0D0,1.0D0/
C
C       READ AND WRITE LOAD IN TIME DOMAIN DOMAIN
C
        READ(3,*) ZLOAD
C       TAKE THE LOADS ON THE ACTUAL FREE SURFACE NODES ZERO
        DO 68 II=1,N1
          JREF=2*(II-1)
          IF(NDTYPE(II).EQ.2) THEN
            ZLOAD(JREF+1)=DCMPLX(0.0D0,0.0D0)
            ZLOAD(JREF+2)=DCMPLX(0.0D0,0.0D0)
          ENDIF
  68    CONTINUE
        WRITE(6,100) NCOUNT
        WRITE(6,103) (ZLOAD(I),I=1,N2)
C
        NOD=3
        GSCALE=RMU
        REWIND 7
        REWIND 8
        DO 2 ICL=1,NNODE
        IROW=2*(ICL-1)
        DO 3 JCL=1,NELEM
        DO 4 K1=1,NOD
        LN=KONN(JCL,K1)
        JCOL=2*(LN-1)
        READ(7,*) ((ZKU(K1,K2,K3),K2=1,2),K3=1,2)
        READ(8,*) ((ZKT(K1,K2,K3),K2=1,2),K3=1,2)
        DO 6 I=1,2
        DO 6 J=1,2
        ZGM(IROW+I,JCOL+J)=ZGM(IROW+I,JCOL+J)+
     &         ZKU(K1,I,J)*GSCALE
        ZFM(IROW+I,JCOL+J)=ZFM(IROW+I,JCOL+J)+
     &         ZKT(K1,I,J)
  6     CONTINUE
  4     CONTINUE
  3     CONTINUE
C
C       ADJUST DIAGONAL TERMS USING CONCEPT OF SUMMAT. OF OFF-DIAGONALS
C       DJUMP ARE JUMP TERMS FOR EXTERIOR TYPE PROBLEM
C
        IF (NCOUNT .GT. 1) GO TO 2
        DO 9 I=1,2
        DO 9 J=1,2
C       NODE ON THE FREE SURFACE (SEMINFINITY)
        IF(NDTYPE(ICL).EQ.2) THEN
          ZFM(IROW+I,IROW+J)=ZFM(IROW+I,IROW+J)+
     &       SFDIAG(ICL,I,J)+0.5D0*DJUMP(I,J)
          GO TO 9
        ENDIF
        ZFM(IROW+I,IROW+J)=ZFM(IROW+I,IROW+J)+
     &         SFDIAG(ICL,I,J)+DJUMP(I,J)
  9     CONTINUE
  2     CONTINUE
        REWIND 7
        REWIND 8
C
C       SOLVE FOR UNKNOWN DISPLACEM. AND/OR LOADS-WATCH SCALING FACTOR
C
        DO 10 I=1,N2
        LC=MBCU(I)
        IF (LC .EQ. 1) GO TO 10
        ZLOAD(I)=ZLOAD(I)/GSCALE
  10    CONTINUE
C
        WRITE(9) ZLOAD
C
        CALL BOUND(NCOUNT,ZGM,ZFM,ZLOAD,ZTEMP,ZCOL,ZB,ZPQ,MBCU,
     &       LINV,MINV,N2)
C
C       SORT OUT DISPLACEMENTS AND TRACTIONS-STORE IN ZTEMP,ZLOAD RESP.
C
        BACKSPACE 4
        BACKSPACE 9
        READ(4) ZTEMP
        READ(9) ZLOAD
        DO 12 I=1,N2
        LC=MBCU(I)
        IF (LC .EQ. 0) GO TO 11
        ZSTUMP=ZTEMP(I)
        ZTEMP(I)=ZLOAD(I)
        ZLOAD(I)=ZSTUMP*GSCALE
        GO TO 12
  11    ZLOAD(I)=ZLOAD(I)*GSCALE
  12    CONTINUE
C
C       PRINTING OF TRACTIONS AND DISPLACEM. IN TIME DOMAIN DOMAIN
C
        WRITE(6,101) NCOUNT
        WRITE(6,103) (ZLOAD(I),I=1,N2)
        WRITE(6,102) NCOUNT
        WRITE(6,103) (ZTEMP(I),I=1,N2)
  100   FORMAT(//2X,' LOADING AT STEP=',I3)
  101   FORMAT(//2X,' TRACTIONS AT STEP=',I3)
  102   FORMAT(//2X,' DISPLACEMENTS AT STEP=',I3)
  103   FORMAT(8(2X,E12.6))
        RETURN
        END
C------------------------------------------------
        SUBROUTINE BOUND(NCOUNT,ZGM,ZFM,ZLOAD,ZTEMP,
     &               ZCOL,ZB,ZPQ,MBCU,LINV,MINV,N2)
C
C       MATRIX RE-ARRANGEM. ACCORDING TO ROUNDARY CONDIT. AND SOLUTION
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION ZGM(N2,N2),ZFM(N2,N2),ZLOAD(N2),ZTEMP(N2)
        DIMENSION ZCOL(N2),ZB(N2),ZPQ(N2),MBCU(N2)
        DIMENSION LINV(N2),MINV(N2)
C
C       INTERCHANGE COLUMNS OF ZGM, ZFM MATRICES
C
        NKA=N2
        DO 10 J=1,NKA
        IF (MBCU(J) .EQ. 0) GO TO 10
        DO 11 I=1,NKA
        ZTEMP(I)=ZFM(I,J)
        ZFM(I,J)=-ZGM(I,J)
        ZGM(I,J)=-ZTEMP(I)
  11    CONTINUE
  10    CONTINUE
C
C       WRITE MATRICES ZGM,ZFM, ON TAPE-PIECEMEAL WRITING
C
        DO 12 I=1,NKA
        DO 13 J=1,NKA
        ZPQ(J)=ZGM(I,J)
  13    CONTINUE
        WRITE(12,*) ZPQ
  12    CONTINUE
        IF (NCOUNT .EQ. 1) GO TO 15
        DO 16 I=1,NKA
        DO 17 J=1,NKA
        ZPQ(J)=ZFM(I,J)
  17    CONTINUE
        WRITE(14,*) ZPQ
  16    CONTINUE
  15    CONTINUE
C
C       SOLUTION FOR RESPONSE
C
        IF (NCOUNT .GT. 1) GO TO 18
        CALL RMINV(ZFM,NKA,LINV,MINV)
        WRITE(11,*) ZFM
  18    CONTINUE
        CALL SOLSCE(NCOUNT,ZFM,ZLOAD,ZTEMP,ZCOL,ZB,ZPQ,N2)
        RETURN
        END
C------------------------------------------------
        SUBROUTINE RMINV(A,N,L,M)
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(N*N),L(N),M(N)
C
        DO 5 I=1,N
          L(N)=0
          M(N)=0
  5     CONTINUE
        NLIN=6
        NK=-N
        DO 180 K=1,N
        NK=NK+N
        L(K)=K
        M(K)=K
        KK=NK+K
        BIGA=A(KK)
        DO 20 J=K,N
        IZ=N*(J-1)
        DO 20 I=K,N
        IJ=IZ+I
        CA=A(IJ)
        RA=DABS(CA)
        RB=DABS(BIGA)
        IF (RB-RA) 10,20,20
  10    BIGA=A(IJ)
        L(K)=I
        M(K)=J
  20    CONTINUE
C
C       INTERCHANGE ROWS
C
        J=L(K)
        IF (J-K) 50,50,30
  30    KI=K-N
        DO 40 I=1,N
        KI=KI+N
        HOLD=-A(KI)
        JI=KI-K+J
        A(KI)=A(JI)
  40    A(JI)=HOLD
C
C       INTERCHANGE COLUMNS
C
  50    I=M(K)
        IF (I-K) 80,80,60
  60    JP=N*(I-1)
        DO 70 J=1,N
        JK=NK+J
        JI=JP+J
        HOLD=-A(JK)
        A(JK)=A(JI)
  70    A(JI)=HOLD
  80    IF (DABS(BIGA)-1.0D-15) 90,100,100
  90    CONTINUE
        WRITE (NLIN,95) K
  95    FORMAT(20X,'DURING A MATRIX INVERSION OPERATION A ZERO',/
     &       20X,' DETERMINANT HAS BEEN DETECTED. THE RUN HAS
     &       BEEN',/' TERMINATED. ROW OR COLUMN NO. WAS ',I4//)
        STOP
  100   DO 120 I=1,N
        IF (I-K) 110,120,110
  110   IK=NK+I
        A(IK)=A(IK)/(-BIGA)
  120   CONTINUE
C
C       REDUCE MATRIX
C
        DO 150 I=1,N
        IK=NK+I
        HOLD=A(IK)
        IJ=I-N
        DO 150 J=1,N
        IJ=IJ+N
        IF (I-K) 130,150,130
  130   IF (J-K) 140,150,140
  140   KJ=IJ-I+K
        A(IJ)=HOLD*A(KJ)+A(IJ)
  150   CONTINUE
C
C       DEVIDE BY PIVOT
C
        KJ=K-N
        DO 170 J=1,N
        KJ=KJ+N
        IF (J-K) 160,170,160
  160   A(KJ)=A(KJ)/BIGA
  170   CONTINUE
C
C       PRODUCT OF PIVOTS
C       REPLACE PIVOT BY RECIPROCAL
C
        A(KK)=1.0D0/BIGA
  180   CONTINUE
C
C       FINAL ROW AND COLUMN INTERCHANGE
C
        K=N
  190   K=(K-1)
        IF (K) 260,260,200
  200   I=L(K)
        IF (I-K) 230,230,210
  210   JQ=N*(K-1)
        JR=N*(I-1)
        DO 220 J=1,N
        JK=JQ+J
        HOLD=A(JK)
        JI=JR+J
        A(JK)=-A(JI)
  220   A(JI)=HOLD
  230   J=M(K)
        IF (J-K) 190,190,240
  240   KI=K-N
        DO 250 I=1,N
        KI=KI+N
        HOLD=A(KI)
        JI=KI-K+J
        A(KI)=-A(JI)
  250   A(JI)=HOLD
        GO TO 190
  260   RETURN
        END
C------------------------------------------------
        SUBROUTINE VPROD(ZA,ZB,N,ZC)
C
C       PRODUCT OF REAL MATRIX TIMES REAL VECTOR
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION ZA(N,N),ZB(N),ZC(N)
C
        DO 20 I=1,N
        ZC(I)=0.0D0
        DO 10 J=1,N
        ZC(I)=ZC(I)+ZA(I,J)*ZB(J)
  10    CONTINUE
  20    CONTINUE
        RETURN
        END
C------------------------------------------------
        SUBROUTINE HOOP(KONN,NONN,ZLOAD,ZTEMP,ZLOAD1,
     &    XCR,YCR,N1,N2,N3,N4,N5,ZZ1,ZZ2,ZZ3)
C
C       HOOP STRESSES IN THE LAPLACE DOMAIN
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION KONN(N3,N4),ZLOAD(N2),ZTEMP(N2),XCR(N1),YCR(N1)
        DIMENSION ZMA(7,7),ZV(7),ZW(7)
        DIMENSION ZTR(2),ZDI(2)
        DIMENSION LDIR(3),XY(3,2),SPOINT(3),SHP(3),DSHP(3),RN(2),RQ(2)
        DIMENSION ZZ1(N4,N3),ZZ2(N4,N3),ZZ3(N4,N3)
        DIMENSION ZLOAD1(N5)
        DIMENSION NONN(N1,N4)
        DIMENSION LINV(7),MINV(7)
        COMMON/CONTR/NNODE,NELEM,NP
        COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM
        COMMON/GAUSS/GPN(10,10),WTN(10,10),GPS(10,10),WTS(10,10)
        COMMON/FREE/NFFS,NFRE
C
        NOD=3
        DATA SPOINT/-1.0D0,0.0D0,1.0D0/
C
C       LOOP OVER ELEMENTS
        READ(2,*) ZLOAD1
        DO 3 JCL=1,NELEM
        DO 4 I=1,NOD
        LN=KONN(JCL,I)
        LDIR(I)=LN
        XY(I,1)=XCR(LN)
        XY(I,2)=YCR(LN)
  4     CONTINUE
C       FIND THE LENGTH OF THE ELEMENT WITH NODES PREVIUS,CURRENT AND NEXT
        NGP=6
        CALL ELEMLEN(XY,NGP,NOD,ELEG)
        DO 6 I=1,NOD
        CALL CLRARR(ZMA,49)
        CALL CLRARR(ZV,7)
        CALL CLRARR(ZW,7)
        EKS=SPOINT(I)

        CALL CSHAPE(EKS,NOD,XY,SHP,DSHP,RN,DETJ)
        RQ(1)=-RN(2)
        RQ(2)=RN(1)
        DO 7 J=1,NOD
        JREF=2*(LDIR(J)-1)
        ZTR(1)=ZLOAD(JREF+1)
        ZTR(2)=ZLOAD(JREF+2)
        ZDI(1)=ZTEMP(JREF+1)
        ZDI(2)=ZTEMP(JREF+2)
        ZV(2)=ZV(2)+SHP(J)*ZTR(1)
        ZV(3)=ZV(3)+SHP(J)*ZTR(2)
        ZV(5)=ZV(5)+DSHP(J)*ZDI(1)
        ZV(7)=ZV(7)+DSHP(J)*ZDI(2)
  7     CONTINUE
        ZV(5)=ZV(5)/(0.5D0*ELEG)
        ZV(7)=ZV(7)/(0.5D0*ELEG)
        ZMA(1,1)=-1.0D0
        ZMA(1,4)=RLA+2.0D0*RMU
        ZMA(1,7)=RLA
        ZMA(2,1)=RN(1)
        ZMA(2,3)=RN(2)
        ZMA(3,2)=RN(2)
        ZMA(3,3)=RN(1)
        ZMA(4,3)=-1.0D0
        ZMA(4,5)=RMU
        ZMA(4,6)=RMU
        ZMA(5,4)=RQ(1)
        ZMA(5,6)=RQ(2)
        ZMA(6,2)=-1.0D0
        ZMA(6,4)=RLA
        ZMA(6,7)=RLA+2.0D0*RMU
        ZMA(7,5)=RQ(1)
        ZMA(7,7)=RQ(2)
        CALL RMINV(ZMA,7,LINV,MINV)
        CALL VPROD(ZMA,ZV,7,ZW)
        Z11=ZW(1)
        Z22=ZW(2)
        Z12=ZW(3)
        ZSRR=Z11*RN(1)**2+Z22*RN(2)**2+
     &             2*Z12*RN(1)*RN(2)
        ZSTT=Z11*RN(2)**2+Z22*RN(1)**2-
     &             2*Z12*RN(1)*RN(2)
        ZSRT=(Z22-Z11)*RN(1)*RN(2)+
     &             Z12*(RN(1)**2-RN(2)**2)
        JREF=3*(LDIR(I)-1)
        ZZ1(I,JCL)=ZSRR+ZLOAD1(JREF+1)
        ZZ2(I,JCL)=ZSTT+ZLOAD1(JREF+2)
        ZZ3(I,JCL)=ZSRT+ZLOAD1(JREF+3)
  6     CONTINUE
  3     CONTINUE
        RETURN
        END
C -------------------------------------------------------------------------
      SUBROUTINE ELEMLEN(XY,NGP,NOD,ELEG)
C      
C     CALCULATES THE LENGTH OF THE ELEMENT, BY GAUSSIAN INTERPOLATION
C     NOD: 3 FOR QUADRATIC INTERPOLATION
C          2 FOR LINEAR INTERPOLATION
C     XY : COORDINATES OF THE NODES OF THE ELEMENT
C     NGP: NUMBER OF POINTS FOR GAUSSIAN INTERPOLATION
C     ELEG: RETURN VALUE OF THE LENGTH OF THE ELEMENT
C
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT REAL*8 (A-H,O-Y)
      DIMENSION XY(3,2),SHP(3),DSHP(3),RN(2)
      COMMON/GAUSS/GPN(10,10),WTN(10,10),GPS(10,10),WTS(10,10)

C
      ELEG=0.0D0
      DO 8 IX=1,NGP
        EK=GPN(NGP,IX)
        WT=WTN(NGP,IX)
        CALL CSHAPE(EK,NOD,XY,SHP,DSHP,RN,DETJ)
        ELEG=ELEG+DETJ*WT
 8    CONTINUE
      RETURN
      END
C------------------------------------------------
C
        SUBROUTINE SOLSCE(NCOUNT,ZFM,ZLOAD,ZTEMP,ZCOL,ZB,ZPQ,N2)
C
C       TIME CONVOLUTIONS
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION ZFM(N2,N2),ZLOAD(N2),ZTEMP(N2)
        DIMENSION ZCOL(N2),ZB(N2),ZPQ(N2)
C
        NKA=N2
C
        REWIND 12
        IF (NCOUNT .GT. 1) REWIND 14
        IF (NCOUNT .GT. 1) BACKSPACE 4
        BACKSPACE 9
C
C       FORM PRODUCT OF ZGM*ZLOAD
C
        CALL CLRARR(ZB,NKA)
        DO 10 L=1,NCOUNT
C       ^^^^^^^^^^^^^^^^^^^^^^
        CALL CLRARR(ZTEMP,N2)
        IF(L.EQ.NCOUNT) GO TO 25
        BACKSPACE 9
        READ(9) ZTEMP
  25    CONTINUE
        READ(9) ZLOAD
        DO 26 II=1,N2
          ZLOAD(II)=(ZLOAD(II)+ZTEMP(II))/2.0D0
  26    CONTINUE
C       ^^^^^^^^^^^^^^^^^^^^^^.
        BACKSPACE 9
        CALL GTPROD(ZPQ,ZLOAD,ZCOL,NKA)
        CALL ADD(ZCOL,ZB,NKA)
        IF (L .NE. NCOUNT) BACKSPACE 9
  10    CONTINUE
C
C       RE-SET UNIT 9
C
        DO 11 L=1,NCOUNT
        READ(9) ZLOAD
  11    CONTINUE
C
C       FORM PRODUCT OF ZFM*ZTEMP
C
        MCOUNT=NCOUNT-1
        IF (MCOUNT .EQ. 0) GO TO 30
        DO 20 L=1,MCOUNT
C       ^^^^^^^^^^^^^^^^^^^^^^
        CALL CLRARR(ZLOAD,N2)
        IF(L.EQ.MCOUNT) GO TO 29
        BACKSPACE 4
        READ(4) ZLOAD
  29    CONTINUE
        READ(4) ZTEMP
C        READ(4,*) ZTEMP
        DO 36 II=1,N2
          ZTEMP(II)=(ZTEMP(II)+ZLOAD(II))/2.0D0
  36    CONTINUE
C       ^^^^^^^^^^^^^^^^^^^^^^.
        BACKSPACE 4
        CALL FUPPRO(ZPQ,ZTEMP,ZCOL,NKA)
        CALL DIFF(ZCOL,ZB,NKA)
        IF (L .NE. MCOUNT) BACKSPACE 4
  20    CONTINUE
C
C       RE-SET UNIT 4
C
        DO 21 L=1,MCOUNT
        READ(4) ZTEMP
  21    CONTINUE
  30    CONTINUE
C
C       POST-MULTIPLY RHS BY ZFM INVERSE
C
        REWIND 11
        READ(11,*) ZFM
        CALL VPROD(ZFM,ZB,NKA,ZTEMP)
C
        WRITE(4) ZTEMP
C
        RETURN
        END
C------------------------------------------------
        SUBROUTINE GTPROD(A,B,C,N)
C
C       PRODUCT OF A*B RETURNED IN C-PIECEMEAL READING OF A
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(N),B(N),C(N)
C
        DO 10 I=1,N
        C(I)=0.0D0
        READ(12,*) A
  9     DO 10 J=1,N
        C(I)=C(I)+A(J)*B(J)
  10    CONTINUE
        RETURN
        END
C------------------------------------------------
        SUBROUTINE FUPPRO(A,B,C,N)
C
C       PRODUCT OF A*B RETURNED IN C - PIECEMEAL READING OF A
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(N),B(N),C(N)
C
        DO 10 I=1,N
        C(I)=0.0D0
        READ(14,*) A
  9     DO 10 J=1,N
        C(I)=C(I)+A(J)*B(J)
  10    CONTINUE
        RETURN
        END
C------------------------------------------------
        SUBROUTINE ADD(A,B,N)
C
C       SUM OF VECTORS B AND A RETURNED IN B
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(N),B(N)
        DO 10 I=1,N
        B(I)=B(I)+A(I)
  10    CONTINUE
        RETURN
        END
C------------------------------------------------
        SUBROUTINE DIFF(A,B,N)
C
C       DIFFERENCE OF VECTORS A AND B RETURNED IN B
C
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(N),B(N)
        DO 10 I=1,N
        B(I)=B(I)-A(I)
  10    CONTINUE
        RETURN
        END
