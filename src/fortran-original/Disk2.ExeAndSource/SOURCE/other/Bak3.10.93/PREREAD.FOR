C$DEBUG
C   --------------------------------------------------------------------
      SUBROUTINE TITLES
C
C     READS AND PRINTS TITLE AND DISCRIPTION OF THE PROBLEM
C
      COMMON/BLK20/DATE,TIME,SNAME,FNAMEIN,FNAMEPRE,
     &             FNAMEPRO,FNAMEPOS
      CHARACTER DATE*35,TIME*15,SNAME*35,FNAMEIN*14,FNAMEPRE*14
      CHARACTER FNAMEPRO*14,FNAMEPOS*14
      CHARACTER TITLE*80
C
C     PRINT TITLE, DATE, TIME, AND PERSON RUNNING PROGRAM
   28 WRITE(6,110)
  110 FORMAT(30X,'** INDATA **',//
     &	     36X,' BY',/
     &	     29X,'KOSTOGLOU KONSTANTINOS')
C *** WRITE(*,111)
      WRITE(6,111)
C  111 FORMAT(1H1,//,10X,'   --INPUT DATA FOR BOUNDARY ',
C     &                  'ELEMENT PROGRAM --',//)
  111 FORMAT(//,10X,'   --INPUT DATA FOR BOUNDARY ',
     &                  'ELEMENT PROGRAM --',//)
       WRITE(6,200)DATE,TIME,SNAME,FNAMEIN,FNAMEPRE,
     &             FNAMEPRO,FNAMEPOS
  200 FORMAT(10X,'RUN DATE:                       ',A,/
     &	     10X,'TIME OF RUN:                    ',A,/
     &	     10X,'RUN BY:                         ',A,/
     &	     10X,'INPUT DATA FILENAME:            ',A,/
     &	     10X,'OUTPUT PREPROCESSING FILENAME:  ',A,/
     &	     10X,'OUTPUT PROCESSING FILENAME:     ',A,/
     &       10X,'OUTPUT POSTPROCESSING FILENAME: ',A)
       WRITE(6,300)
  300 FORMAT(//,10X,'DATA IN INTERNATIONAL UNIT SYSTEM')

C      IF(PLOTFILE.NE.'NONE') WRITE(6,201)PLOTFILE
C 201  FORMAT(10X,'PLOTTED OUTPUT FILENAME:  ',A)
C ------------------------------------
C   READ AND PRINT PROBLEM DESCRIPTION
C ------------------------------------
      READ(7,113) TITLE
C 113  FORMAT(20A4)
 113  FORMAT(A)
C *** WRITE(*,112) TITLE
      WRITE(6,112) TITLE
C  112 FORMAT(////,2X,'PROBLEM DISCRIPTION ',10A4,/32X,10A4,//)
C  112 FORMAT(////,'TITLE',/20A4)
 112  FORMAT(////,'TITLE',/A)


      RETURN
      END
C   --------------------------------------------------------------------
      SUBROUTINE READER(SUM,ISUM,IREAD,IR,I)
C
C   FUNCTIONS -
C
C        ENABLES PROGRAM TO READ NUMERICAL DATA, CONSISTING OF POSITIVE
C        INTEGERS AND BOTH NEGATIVE AND POSITIVE REAL NUMBERS
C        FREE-FORM.
C          IR=0 : READ ONLY INTEGER POSITIVE NUMBERS
C                 RETURN VALUE (ISUM) IS INTEGER
C          IR=1 : READ POSSITIVE OR NEGATIVE INTEGERS OR REALS,
C                 RETURN VALUE (SUM) IS REAL
C          IR=2 : READ ONLY REALS,
C                 RETURN VALUE (SUM) IS REAL
C          IR=3 : READ ONLY INTEGER POSITIVE AND NEGATIVE NUMBERS,
C                 RETURN VALUE (ISUM) IS INTEGER
C
C          IREAD=0 : READS A NEW  WHOLE LINE FROM THE DATA FILE
C          IREAD=1 : READS FROM THE SAME LINE
C                   WHEN THE PROGRAM IS CALLING WITH THE VALUE 
C                   OF IREAD 1 OR 0, IT RETURNS BACK THE VALUE 1,
C                   SO IN THE NEXT CALLING OF THE SUBROUTINE  IT WILL 
C                   READ AGAIN FROM THE SAME LINE, UNLESS WE GIVE 
C                   IREAD=0 BEFORE CALLING THE SUBROUTINE.
C          IEXIT   : IF AN ERROR OCCURS THEN IEXIT IS TAKEN THE VALUE
C                  1: THE PROGRAM CONTINUES EXECUTING
C                  2: THE PROGRAM TERMINATES WITH AN ERROR MESSAGE
C          I:  A RETURNED INTEGER THAT SHOWS IN WHICH COLUMN OF THE ROW
C              M(I) THE PROCESSING STOPS. THIS IS USED IN THE CASE THAT
C              WE HAVE IREAD=1, READING FROM THE SAME LINE.
C
C
      COMMON/BLKEX/IEXIT
C      DIMENSION X(200)
C      DIMENSION X(80)



C      DIMENSION M(80),N(14),A(10),IA(10),ERROR(7)
C     ATTENTION    !!!!!!!!!!!!!! LETTER 'E' IS INTRODUCED
      DIMENSION M(80),N(16),ERROR(7)
      DIMENSION AA(10), IA(10)
      CHARACTER N*1,ERROR*50
      CHARACTER M*1
      REAL*4 SUM,SUM1,AA
C      EQUIVALENCE (M(1),X(1))



C          ATTENTION   !!!!!!!!! LETTER 'E' IS INTRODUCED
      DATA N/'1','2','3','4','5','6','7','8','9','0',' ','.','O',
     &       '-','E','*'/
      DATA AA/1.,2.,3.,4.,5.,6.,7.,8.,9.,0./
      DATA IA/1,2,3,4,5,6,7,8,9,0/
C      DATA ERROR/'FR01','FR02','FR03','FR04','FR05','FR06','FR07'/
C------------------------------------------------------------
C     NOT VALID ERRORS
C------------------------------------------------------------
      ERROR(1)='  NOT VALID - LINE OF DATA IN FILE CORRUPTED '
      ERROR(2)='  NOT VALID - FIRST COLUMN EMPTY '
      ERROR(3)='  NOT VALID - DECIMAL POINT IN AN INTEGER '
      ERROR(4)='  NOT VALID - NEGATIVE INTEGERS NOT PERMITED '
      ERROR(5)='  NOT VALID - ILLEGAL CHARACTER  IN AN INTEGER '
      ERROR(6)='  NOT VALID - DECIMAL POINT MUST BE IN A REAL '
      ERROR(7)='  NOT VALID - ILLEGAL CHARACTER IN A REAL '
      SSIGN=1.
      ISIGN=1
      ID=0
      ISUM=0
      SUM=0.0
      IF(I.GE.80)GO TO 66
      IF((IREAD.EQ.1).AND.(IR.EQ.1))GO TO 3
      IF((IREAD.EQ.1).AND.(IR.EQ.2))GO TO 3
      IF(IREAD.EQ.1)GO TO 11
C -------------------------------------------------
C   READ ENTIRE LINE OF DATA - ONE COLUMN AT A TIME
C -------------------------------------------------
C     CLEAR ARRAY M
C      DO II=1,80
C	M(II)=' '
C      END DO

    5 READ(7,101,END=23)(M(I),I=1,80)
C  101 FORMAT(80A1)
  101 FORMAT(80A)


C     IF THE CHARACTER '*' IS FOUND IN THE BEGINING OF THE LINE,
C     ACCEPT THIS LINE AS COMMENT, RETURN
      IF(M(1).EQ.N(16)) GO TO 66
      GO TO 24
   23 WRITE(6,105)ERROR(1)
C *** WRITE(*,105)ERROR(1)
  105 FORMAT(/,10X,'****  ERROR - ',A,' ****',/)
      CALL QUIT
C   24 IF(M(1).NE.N(11))GO TO 22
CC *** WRITE(*,105)ERROR(2)
C      WRITE(6,105)ERROR(2)
CC *** WRITE(*,107)(M(II),II=1,80)
C      WRITE(6,107)(M(II),II=1,80)


  107 FORMAT(1X,'****  ERROR IN FOLLOWING LINE OF DATA  ****',//,
     &	 1X,1HC,/,1X,80A1,//)

  24   GO TO 22
CC ------------------------------------
CC   TERMINATE EXECUTION IF DATA FAULTY
CC ------------------------------------
C      CALL QUIT
   22 IREAD=1
      K=1
      I=0
C      IF(IR.EQ.1)GO TO 4
      IF(IR.EQ.1)GO TO 64
C      GO TO 12
C   INTEGER NUMBERS
C   =========================
C -------------------------
C   ASSEMBLE INTEGER NUMBER
C ------------------------- 
C   11 IF(M(I).EQ.N(11).AND.M(I+1).EQ.N(11)) THEN
C        I=I+2
C        IF(I.GE.80)GO TO 13
C        GO TO 11
C       ENDIF
      I=I+1
C     IGNORE THE EMPTY CHARACTERS IN THE BEGINING OF THE LINE
   11 IF(M(I).EQ.N(11)) THEN
        I=I+1

C^^^^^^^^^^^^^
C        IF(I.GE.80)GO TO 13
        IF(I.GE.80)GO TO 33


        GO TO 11
       ENDIF
C      K=I+1
      K=I
      GO TO 71
   12 I=I+1
   71 CONTINUE
C -------------------------------------------
C   IF LETTER (O) READ, SUBSTITUTE NUMBER (0)
C -------------------------------------------
      IF(M(I).EQ.N(13))M(I)=N(10)
C -------------------------
C   CHECK FOR DECIMAL POINT
C -------------------------
      IF(M(I).NE.N(12))GO TO 16
C *** WRITE(*,105)ERROR(3)
      WRITE(6,105)ERROR(3)
      IEXIT=1
C -----------------------
C   CHECK SIGN OF INTEGER
C -----------------------
   16 IF(M(I).NE.N(14))GO TO 20
C     ATTENTION      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     CASE POSSITIVE AND NEGATIVE INTEGERS
      IF ((M(I).EQ.N(14)).AND.(IR.EQ.3)) THEN
C        I=I+1
        ISIGN=-1
        GO TO 12
       ENDIF
C *** WRITE(*,105)ERROR(4)
      WRITE(6,105)ERROR(4)
      IEXIT=1
C ------------------------------------
C   TERMINATE EXECUTION IF DATA FAULTY
C ------------------------------------
   20 IF(M(I).NE.N(11))GO TO 12
      IF(IEXIT.LE.1) GOTO 26
C *** WRITE(*,107)(M(II),II=1,80)
      WRITE(6,107)(M(II),II=1,80)
CCC     ENDIF
      CALL QUIT
   26 KI=K
      K=I+1
      KF=I-1
      ISUM=0
C ------------------------------------
C   DETERMINE DIGITS OF INTEGER NUMBER
C ------------------------------------
      DO 13 II1=KI,KF
C        ATTENTION    !!!!!!!!!!!!!!!!
C     CASE POSSITIVE AND NEGATIVE INTEGERS
      IF((IR.EQ.3).AND.(M(II1).EQ.N(14)))GO TO 13
      DO 14 J=1,10
      L=J
      IF(M(II1).EQ.N(J))GO TO 15
   14 CONTINUE
C ------------------------------------
C   TERMINATE EXECUTION IF DATA FAULTY
C ------------------------------------
C *** WRITE(*,105)ERROR(5)
      WRITE(6,105)ERROR(5)
C *** WRITE(*,107)(M(IJ),IJ=1,80)
      WRITE(6,107)(M(IJ),IJ=1,80)
      CALL QUIT
   15 ISUM=ISUM+IA(L)*10**(KF-II1)*ISIGN
   13 CONTINUE

C^^^^^^^^^^^^^
   33 CONTINUE



      RETURN
C   REAL OR INTEGER (POSITIVE OR NEGATIVE) NUMBERS
C=================================================
C ----------------------
C   ASSEMBLE REAL OR INTEGER NUMBER
C ----------------------
C    3 IF(M(I).EQ.N(11).AND.M(I+1).EQ.N(11)) THEN
C        I=I+2
C       IF(I.GE.80)GO TO 66
C        GO TO 3
C       ENDIF


    3 I=I+1


C      IGNORE THE EMPTY CHARACTERS IN THE BEGINING OF THE LINE
   55  IF(M(I).EQ.N(11)) THEN
        I=I+1
       IF(I.GE.80)GO TO 66
        GO TO 55
       ENDIF

       K=I
C      INITIALIZE COUNTER FOR LETTER E 
       KES=0
       KEF=0


      GO TO 73
   64 I=I+1 

C     IGNORE THE EMPTY CHARACTERS IN THE BEGINING OF THE LINE
   81 IF(M(I).EQ.N(11)) THEN
        I=I+1
        IF(I.GE.80)GO TO 66
        GO TO 81
      ENDIF

C      K=I+1
      K=I
      GO TO 73
    4 I=I+1
   73 CONTINUE
C -------------------------------------------
C   IF LETTER (O) READ, SUBSTITUTE NUMBER (0)
C -------------------------------------------
      IF(M(I).EQ.N(13))M(I)=N(10)
C     IF DECIMAL POINT FOUND
      IF(M(I).EQ.N(12))ID=I
C     LETTER E IS ENCOUNTERED
C     STARTING POINTER OF LETTER E
      IF (M(I).EQ.N(15)) KES=I
C     FINISHING POINTER OF EXPONENTIAL NUMBER
      IF (KES.GT.0) KEF=I
C
      IF(M(I).NE.N(11))GO TO 4
      KI=K
      K=I+1
C     FINISHING POINTER OF EXPONENTIAL NUMBER
      IF (KEF.GT.0) KEF=KEF-1
C ------------
C   CHECK SIGN
C ------------
      IF(M(KI).NE.N(14))GO TO 19
      KI=KI+1
      SSIGN=-1.
   19 KF=I
      IF (IR.EQ.1) THEN
C       DECIMAL POINT NOT FOUND
        IF (ID.EQ.0)  THEN
          ID=I
C         LETTER E HAD ENCOUNTERED
          IF (KES.GT.0) ID=KES-1
        ENDIF
      ENDIF
      ID1=ID-1
      SUM=0.
      SUM1=1.0
C -------------------------
C   CHECK FOR DECIMAL POINT
C -------------------------
       DO 17 II1=KI,KF 
C        REAL OR INTEGER NUMBER
         IF (IR.EQ.1) GO TO 18
C         REAL NUMBER ONLY
         IF (IR.EQ.2) THEN
C          SEEK UNTIL FIND DECIMAL POINT
           IF(M(II1).EQ.N(12))GO TO 18
         ENDIF
   17  CONTINUE
C  *** WRITE(*,105)ERROR(6)
       WRITE(6,105)ERROR(6)
C ------------------------------ 
C   CHECK FOR ILLEGAL CHARACTERS
C ------------------------------
      DO 25 II1=KI,KF
C      DO 21 J=1,14 
      DO 21 J=1,15
        IF(M(II1).EQ.N(J))GO TO 25
   21 CONTINUE
C *** WRITE(*,105)ERROR(7)
      WRITE(6,105)ERROR(7)
   25 CONTINUE
C *** WRITE(*,107)(M(II),II=1,80)
      WRITE(6,107)(M(II),II=1,80)
C ------------------------------------
C   TERMINATE EXECUTION IF DATA FAULTY
C ------------------------------------
      CALL QUIT
C ---------------------------------
C   DETERMINE DIGITS OF REAL NUMBER
C ---------------------------------
   18 DO 2 II1=KI,KF-1
C       COUNTER EQUAL DECIMAL POINT COUNTER
C       AND NOT EQUAL EXPONENTIAL LETTER E
        IF((II1.GE.ID).AND.(KES.EQ.0))GO TO 6
        IF(KES.GT.0) THEN
C         COUNTER EQUALS EXPONENTIAL POINT E COUNTER
          IF((II1.GE.ID).AND.(II1.LT.KES))GO TO 6 
          IF(II1.EQ.KES) GO TO 2
          IF(II1.GT.KES)THEN
            KESF=KEF-II1 
            DO 62 KK=1,10
              IF(M(II1).EQ.N(KK))THEN
                POW=AA(KK)*10.**KESF
                SUM1=SUM1*10.**POW
              ENDIF
   62       CONTINUE
            IF(II1.EQ.KEF) SUM=SUM*SUM1
            GO TO 2
          ENDIF
        ENDIF
        DO 7 J=1,10
          L=J
          IF(M(II1).EQ.N(J))GO TO 8
    7   CONTINUE
C       ------------------------------------
C       TERMINATE EXECUTION IF DATA FAULTY
C       ------------------------------------
C ***   WRITE(*,105)ERROR(7)
        WRITE(6,105)ERROR(7)
C ***   WRITE(*,107)(M(IJ),IJ=1,80)
        WRITE(6,107)(M(IJ),IJ=1,80)
        CALL QUIT
C    8   SUM=SUM+AA(L)*10.**(ID1-II1)
    8   CONTINUE
        IK1=ID1-II1
        RIK1=FLOAT(IK1)
        SUM=SUM+AA(L)*10.0**RIK1



        GO TO 2
    6   IF(II1.EQ.ID)GO TO 2
        DO 9 J=1,10
          L=J
          IF(M(II1).EQ.N(J))GO TO 10
    9   CONTINUE
C   10   SUM=SUM+AA(L)/10.**(II1-ID)
   10   CONTINUE
        IK1=II1-ID
        RIK1=FLOAT(IK1)
        SUM=SUM+AA(L)*10.0**(-RIK1)
    2 CONTINUE
   52 CONTINUE
C     ------------------
C     ACCOUNT FOR SIGN
C     ------------------
      SUM=SUM*SSIGN
  66  CONTINUE
      RETURN
      END
C   -------------------------------------------------------------------
      SUBROUTINE READTITL(SUBTITL,NT)
C
C   FUNCTIONS -
C
C        ENABLES PROGRAM TO READ A LINE CONSISTING OF MANY WORDS.
C        ONLY THE FIRST 4 CHARACTERS OF THE WORDS ARE READING.
C        RETURNS NT SUBTITL OF 4 CHARACTER LENTH EACH
C        E.G. IF THE LINE IS:
C        AXIS.NUMBER STRESS.AXIAL SHEAR.STRESS
C        THE SUBROUTINE WILL RETURN 4 SUBTITL: 'AXIS','STRE','SHEA'
C
      DIMENSION M(80),ERRORS(7)
      DIMENSION M1(4),NC(2)
      DIMENSION SUBTITL(10)
      CHARACTER NC*1,ERRORS*4
      CHARACTER M*1
      CHARACTER M1*1
      CHARACTER SUBTITL*4
      DATA NC/' ','*'/
      DATA ERRORS/'FR01','FR02','FR03','FR04','FR05','FR06','FR07'/
      KDIAK=0
C -------------------------------------------------
C   READ ENTIRE LINE OF DATA
C -------------------------------------------------
    5 READ(7,101,END=23)(M(I),I=1,80)  
  101 FORMAT(80A1)

C     IF THE CHARACTER '*' IS FOUND IN THE BEGINING OF THE LINE,
C     ACCEPT THIS LINE AS COMMENT, RETURN
      IF(M(1).EQ.NC(2)) GO TO 764





C      IF(M(1).NE.NC(1)) THEN




C -------------------------------------------------
C   PROCESS THIS LINE AND EXTRACT THE SUBTITLES, WHITCH ARE DIVIDED BY A SPACE
C -------------------------------------------------
C        IREAD=1
        DO 759 II=1,10
        DO 759 JJ=1,4
            SUBTITL(II)(JJ:JJ)=' '
 759    CONTINUE

        NT=0
C       K: COUNTER OF THE COLUMN OF THE LINE (FROM 1 TO 80)
C       I: COUNTER IN THE SUBELEMENT (FROM 1 TO 4)
        K=0
        I=0
        K=K+1
 742    I=0
        KDIAK=0
        CALL ICLRIRR(M1,4)
C       IGNORE THE EMPTY CHARACTERS IN THE BEGINING OF THE LINE 
C       OR BETWEEN SUBTITLES
   11   IF(M(K).EQ.NC(1)) THEN
          K=K+1
          IF(K.GE.80)GO TO 764
          GO TO 11
        ENDIF
        GO TO 744

 743    K=K+1
 744    I=I+1

        M1(I)=M(K)
C      CASE OF LENGTH OF SUBTITLE LOWER THAN 4
       IF ((M(K).EQ.NC(1)).AND.(I.LT.4)) THEN
         KDIAK=1
         I=I-1
C         GO TO 753
       ENDIF
C
C   41  IF((M(K).EQ.NC(1)).AND.(M(K+1).EQ.NC(1)).AND.(KDIAK.EQ.0)) THEN
C        I=I+2
   41  IF((M(K).EQ.NC(1)).AND.(KDIAK.EQ.0)) THEN
        I=I+1







        IF(I.GE.80)GO TO 762
        GO TO 41
       ENDIF
       IF ((M(K).EQ.NC(1)).AND.(KDIAK.EQ.0)) GO TO 742
  753  IF ((I.GE.4).OR.(KDIAK.EQ.1)) THEN
          NT=NT+1
          DO 758 II=1,I
            SUBTITL(NT)(II:II)=M1(II)
  758     CONTINUE
  763     IF (K.EQ.80) GO TO 762
C          IF ((M(K).EQ.NC(1)).AND.(M(K+1).EQ.NC(1))) GO TO 761
          IF ((M(K).EQ.NC(1)).OR.(KDIAK.EQ.1)) GO TO 742
  761     K=K+1
          GO TO 763
        ENDIF 
  752   CONTINUE
        GO TO 743
C      ENDIF





C      IF (M(1).EQ.NC(1)) THEN
CC   *** WRITE(*,105)ERRORS(2)
C        WRITE(6,105)ERRORS(2)
CC   *** WRITE(*,107)(M(II),II=1,80)
C        WRITE(6,107)(M(II),II=1,80)
C  107   FORMAT(1X,'****  ERROR IN FOLLOWING LINE OF DATA  ****',//,
C     &   1X,1HC,/,1X,80A1,//)
CC ------------------------------------
CC   TERMINATE EXECUTION IF DATA FAULTY
CC ------------------------------------
C        CALL QUIT
C      ENDIF
  762 GO TO 764
   23 WRITE(6,105)ERRORS(1)
C *** WRITE(*,105)ERRORS(1)
  105 FORMAT(/,10X,14H****  ERROR - ,A4,6H  ****,/)
      CALL QUIT 
  764 CONTINUE
      RETURN
      END
C
C   --------------------------------------------------------------------
      SUBROUTINE QUIT
C
C   FUNCTIONS -
C
C        STOPS EXECUTION OF THE PROGRAM
C
C        TERMINATES PLOTTING DEVICE.
C
C *** WRITE(*,102)
      WRITE(6,102)
  102 FORMAT(//,
     &	10X,'******************************************',/,
     &	10X,'****** EXECUTION OF PREBEML ABORTED ******',/,
     &	10X,'******************************************',/)
      STOP
      END
C   -------------------------------------------------------------------
      SUBROUTINE READFIR(TITL,KCOM)
C
C   FUNCTIONS -
C        READS THE COMMAND TITLES
C        ONLY THE FIRST 4 CHARACTERS OF THE TITLE IS READING.
C        E.G. IF THE LINE IS:
C        MATERIAL.PROPERTIES
C        THE SUBROUTINE WILL RETURN THE TITL: 'MATE'
C
      DIMENSION M(80),ERRORS(7)
      DIMENSION N(2)
C        M1(4),
      CHARACTER N*1,ERRORS*4
      CHARACTER M*1
      CHARACTER TITL*4
      DATA N/' ','*'/
      DATA ERRORS/'FR01','FR02','FR03','FR04','FR05','FR06','FR07'/
      KDIAK=0
C -------------------------------------------------
C   READ ENTIRE LINE OF DATA
C -------------------------------------------------
    5 READ(7,101,END=23)(M(I),I=1,80)  
  101 FORMAT(80A1)



C -------------------------------------------------
C   PROCESS THIS LINE AND EXTRACT THE TITLE
C -------------------------------------------------
C        IREAD=1
        DO 759 JJ=1,4
            TITL(JJ:JJ)=' '
 759    CONTINUE
C     IF THE CHARACTER '*' IS FOUND IN THE BEGINING OF THE LINE,
C     ACCEPT THIS LINE AS COMMENT, RETURN
      IF(M(1).EQ.N(2)) THEN
        TITL(1:1)=M(1)
        KCOM=1
        GO TO 764
      ENDIF 


        NT=0
C       K: COUNTER OF THE COLUMN OF THE LINE (FROM 1 TO 80)
CC       I: COUNTER IN THE SUBELEMENT (FROM 1 TO 4)
        K=0
C        I=0
        K=K+1
C       IGNORE THE EMPTY CHARACTERS IN THE BEGINING OF THE LINE
   11   IF(M(K).EQ.N(1)) THEN
C          I=I+1
          K=K+1
C          IF(I.GE.80)GO TO 13
          IF(K.GE.80)GO TO 764
          GO TO 11
        ENDIF
        DO 21 II=1,4
          IF((K+II).GE.80)GO TO 764
          TITL(II:II)=M(K+II-1)
  21   CONTINUE
       GO TO 764
   23 WRITE(6,105)ERRORS(1)
C *** WRITE(*,105)ERRORS(1)
  105 FORMAT(/,10X,14H****  ERROR - ,A4,6H  ****,/)
      CALL QUIT 
  764 CONTINUE
      RETURN
      END
C   -------------------------------------------------------------------
      SUBROUTINE READFILE(TITL,L,KE)
C
C   FUNCTIONS -
C        READS FROM THE SCREEN, THE NAME OF A FILE WITHOUT EXTENSION, 
C        L CHARACTERS LONG
C        TITL: THE NAME OF THE FILE, RETURN VALUE
C        L   : GIVEN LENGTH OF TITL
C        KE  : POSITION OF THE LAST CHARACTER (WITHOUT THE SPACES), IN
C              TITL, RETURN VALUE
C        ONLY THE FIRST L CHARACTERS OF THE NAME OF THE FILE ARE READING.
C        FOR NAMES LESS THAN 8 CHARACTERS, THE SPACES ARE OMITED.
C        ONLY ALPHABETICAL CHARACTERS ARE ACCEPTED.
C        E.G. IF THE LINE IS:
C        MYDATAFILE
C        THE SUBROUTINE WILL RETURN THE TITL: 'MYDATAFI', K=8
C        OR FOR MYFIL, RETURN VALUE 'MYFIL   ', KE=5
C
      DIMENSION M(80)
      DIMENSION N(2)
C        M1(4),
      CHARACTER N*1
      CHARACTER M*1
      CHARACTER TITL*8
      DATA N/' ','.'/
      KDIAK=0
C -------------------------------------------------
C   READ L CHARACTERS FROM THE LINE ON THE SCREEN
C -------------------------------------------------
C    5 READ(*,101,END=23)(M(I),I=1,L)  
    5 READ(*,101,END=23) TITL
  101 FORMAT(A)
C -------------------------------------------------
C   PROCESS THIS LINE AND EXTRACT THE TITLE
C -------------------------------------------------
C        IREAD=1
      DO 6 JJ=1,80
        M(JJ)=' '
 6    CONTINUE
      DO 7 JJ=1,L
        M(JJ)=TITL(JJ:JJ)
        TITL(JJ:JJ)=' '
 7    CONTINUE
      NT=0
C     KE: COUNTER OF THE LAST READING CHARACTER (WITHOUT SPACES)
C     K : COUNTER IN THE REDING ATRING M
      KE=0
      K=0
  11  K=K+1
C     IGNORE THE EMPTY CHARACTERS IN THE LINE, IF 2 EMPTY CHARACTERS
C     ARE ENCOUNTERED, STOP THE READING.
      IF((M(K).EQ.N(1)).AND.(M(K+1).EQ.N(1))) THEN
        K=K-1
        GO TO 34
      ENDIF
      IF(M(K).EQ.N(1)) THEN
C        K=K+1
        IF(K.GT.L)GO TO 34
        GO TO 11
      ENDIF
      KE=KE+1
      TITL(KE:KE)=M(K)
      IF(M(K).EQ.N(2)) GO TO 23
      IF(K.GT.L) GO TO 34
      GO TO 11
  23  WRITE(6,105)
      WRITE(*,105)
      CALL QUIT 
  105 FORMAT(/,10X,'****  ERROR - GIVE THE DATA FILE WITHOUT ',
     &                   'EXTENSION  ****',/)
  34  CONTINUE
      RETURN
      END
