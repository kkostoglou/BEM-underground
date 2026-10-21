C
$DEBUG
C -----------------------------------------------------------------------
      SUBROUTINE REFLCOND(FAI1,THN1,KREFL,KANGL)
C
C     EXAMINES THE REFLECTION CONDITIONS OF P OR SV PLANE WAVE
C
C     FAI1: ANGLE OF AXIS OF WAVEFRONT WITH X-AXIS (0<=FAI<=PI)
C     THN1: ANGLE OF THE OUTWARD TO THE FREE SURFACE UNIT NORMAL 
C     WITH THE X AXIS
C     KREFL: RETURN VALUE
C            0: NO REFLECTION
C            1: REFLECTION
C            2: CASE THAT THE INCIDENT WAVE IS PARALLEL TO THE FREE SURFACE,
C               IT CANOT BE REFLECTED.
C     KANGL : INDEX FOR CALCULATING THE ANGLES OF INCIDENT AND REFLECTED,
C             RETURN VALUE
C            1: TH0=TH0,     TH2=TH2
C            2: TH0=2PI-TH0, TH2=2PI-TH2
C
C      COMMON/GEOM3/FAI,PSAI
C
C      DIMENSION PSAI(5)
C
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT REAL*8 (A-H,O-Y)
C
      PI=3.14159265359D0
C     TOLERANCE 
      TOL=1.0E-11

C
C             ATTENTIONNNNNNNNNNNNNNNNNNNNN
C                   NO REFLECTIONS FOR WAVE DIRECTION PARALLEL TO FREE SURFACE
C             CHECK IF THE REFLECTION CONDITION IS FULLFILED
C             FOR 0<=THN1<PI/2
C               WE HAVE REFLECTION FOR PI/2<FAI-THN1<3*PI/2          [KREFL=1]
C                 ANGLE TH0 AND TH2 IN LOCAL COORDINATE SYSTEM WITH ORIGIN
C                 AT A NODE AND AXES PARALLEL TO N,Q (N:NORMAL UNIT VECTOR
C                                                     Q:TANGENT UNIT VECTOR)
C                 FOR PI/2<FAI-THN1<PI TH0=TH0, TH2=TH2              [KANGL=1]
C                 FOR PI<FAI-THN1<3PI/2 TH0=2PI-TH0, TH2=2PI-TH2     [KANGL=1]
C             FOR PI/2 <= THN1 <= PI 
C               WE HAVE REFLECTION FOR PI/2 < FAI-THN1 < 2PI-THN1  
C               OR WHEN -PI/2<FAI-THN1<0                             [KREFL=1]
C                 ANGLE TH0 AND TH2 IN LOCAL COORDINATE SYSTEM WITH ORIGIN
C                 AT A NODE AND AXES PARALLEL TO N,Q (N:NORMAL UNIT VECTOR
C                                                     Q:TANGENT UNIT VECTOR)
C                 FOR PI/2<FAI-THN1<PI TH0=TH0, TH2=TH2              [KANGL=1]
C                 FOR PI<FAI-THN1<2PI-THN1, 
C                 OR -PI/2<FAI-THN1<0    TH0=2PI-TH0, TH2=2PI-TH2    [KANGL=2]
C             FOR PI<= THN1 < 3*PI/2
C               WE  HAVE REFLECTION WHEN  PI/2<FAI-THN1<2PI-THN1
C               OR WHEN -PI/2<FAI-THN1<0                             [KREFL=1]
C                 ANGLE TH0 AND TH2 IN LOCAL COORDINATE SYSTEM WITH ORIGIN
C                 AT A NODE AND AXES PARALLEL TO N,Q (N:NORMAL UNIT VECTOR
C                                                     Q:TANGENT UNIT VECTOR)
C                 FOR PI/2<FAI-THN1<2PI-THN1, 
C                 OR -THN1<FAI-THN1<-PI  TH0=TH0, TH2=TH2            [KANGL=1]
C                 FOR -PI<FAI-THN1<0     TH0=2PI-TH0, TH2=2PI-TH2    [KANGL=2]
C             FOR 3PI/2<= THN1 < 2PI
C               WE  HAVE REFLECTION WHEN  -3PI/2<FAI-THN1<-PI/2      [KREFL=1]
C                 ANGLE TH0 AND TH2 IN LOCAL COORDINATE SYSTEM WITH ORIGIN
C                 AT A NODE AND AXES PARALLEL TO N,Q (N:NORMAL UNIT VECTOR
C                                                     Q:TANGENT UNIT VECTOR)
C                 FOR -3PI/2<FAI-THN1<-PI   TH0=TH0, TH2=TH2         [KANGL=1]
C                 FOR -PI<FAI-THN1<-PI/2    TH0=2PI-TH0, TH2=2PI-TH2 [KANGL=2]
C
              DIF1=FAI1-THN1
C
C             CASE THAT THE INCIDENT WAVE IS PARALLEL TO THE FREE SURFACE,
C             IT CANOT BE REFLECTED.
              IF((DABS(DABS(DIF1)-PI/2.0D0).LT.1.0D-6).OR.
     &           (DABS(DABS(DIF1)-3.0D0*PI/2.0D0).LT.1.0D-6))KREFL=2
CCCCCCCCCCCCCCCCCCCCCCC
              IF(KREFL.EQ.2) GO TO 49

C
C
C             0<=THN1<PI/2
              IF((THN1.GE.0.0D0).AND.(THN1.LT.PI/2.0D0)) THEN
                IF((DIF1.GT.PI/2.0D0).AND.(DIF1.LT.3.0D0*PI/2.0D0))
     &                  KREFL=1
                IF((DIF1.GT.PI/2.0D0).AND.(DIF1.LT.PI))
     &                  KANGL=1
                IF((DIF1.GT.PI).AND.(DIF1.LT.3*PI/2.0D0))
     &                  KANGL=2
              ENDIF
C             PI/2<=THN1<PI
              IF((THN1.GE.PI/2.0D0).AND.(THN1.LT.PI)) THEN
                IF((DIF1.GT.PI/2.0D0).AND.(DIF1.LT.(2.0D0*PI-THN1)))
     &               KREFL=1
                IF((DIF1.GT.(-PI/2.0D0)).AND.(DIF1.LT.0.0D0))
     &               KREFL=1
                IF((DIF1.GT.PI/2.0D0).AND.(DIF1.LT.PI))
     &                  KANGL=1
                IF((DIF1.GT.PI).AND.(DIF1.LT.(2.0D0*PI-THN1)))
     &                  KANGL=2
                IF((DIF1.GT.-PI/2.0D0).AND.(DIF1.LT.0.0D0))
     &                  KANGL=2
              ENDIF
C             PI<=THN1<3PI/2
              IF((THN1.GE.PI).AND.(THN1.LT.3.0D0*PI/2.0D0)) THEN
                IF((DIF1.GT.PI/2.0D0).AND.(DIF1.LT.(2.0D0*PI-THN1)))
     &                  KREFL=1
                IF((DIF1.GT.(-PI/2.0D0)).AND.(DIF1.LT.0.0D0))
     &               KREFL=1
                IF((DIF1.GT.PI/2.0D0).AND.(DIF1.LT.(2*PI-THN1)))
     &                  KANGL=1
                IF((DIF1.GT.-THN1).AND.(DIF1.LT.-PI))
     &                  KANGL=1
                IF((DIF1.GT.-PI).AND.(DIF1.LT.0.0D0))
     &                  KANGL=2
              ENDIF
C             3PI/2<=THN1<2PI
              IF((THN1.GE.3.0D0*PI/2.0D0).AND.(THN1.LT.2.0D0*PI))
     &        THEN
                IF((DIF1.GT.(-3.0*PI/2.0D0)).AND.(DIF1.LT.
     &              (-PI/2.0D0))) KREFL=1
                IF((DIF1.GT.-3.0D0*PI/2.0D0).AND.(DIF1.LT.-PI))
     &                  KANGL=1
                IF((DIF1.GT.-PI).AND.(DIF1.LT.-PI/2.0D0))
     &                  KANGL=2
              ENDIF
  49  CONTINUE
      RETURN
      END
