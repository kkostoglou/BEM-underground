C************************************************
C
C
        PROGRAM BEMZONL
C==================================================================
c         50 elements, 100 nodes, 200 time steps
C==================================================================
C       A 2-D BEM PROGRAM IN LAPLACE TRANSFORM DOMAIN
C       WITH INVERSE LAPLACE TRANSFORMATION
c       It calculates the free field stresses on a boundary, by the
c       method of subregions

c       mst : total time steps
c       mnd : total nodes of the mesh
c       mel : total elements of the mesh
c       mbd : total boundaries in the mesh
c       mpl : number of total plot categories
c       mzn : number of total zones of different material
c       mdm : dimensions of the problem (2D here)
c        parameter(mst=200,mnd=100,mel=50,mbd=5,mpl=30,mzn=5,mdm=2)
        parameter(mst=100,mnd=100,mel=50,mbd=5,mpl=30,mzn=5,mdm=2)

C
         IMPLICIT COMPLEX*16 (Z)
         IMPLICIT REAL*8 (A-H,O-Y)
c         DIMENSION IARR(860),RARR(210),ZARR(80500)
c         DIMENSION IARR(1710),RARR(210),ZARR(80500)
         DIMENSION IARR(1+4*mel+6*mnd),RARR(1+2*mnd)
         dimension ZARR(1+8*mnd*mnd+4*mnd)
         COMMON/CONTR/NNODE,NELEM,NP
c                 modification for viscoelesticity
         COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM,fvisc,nmode
         common/free/nffs,nfre
c       needed for fast fourier transformation
c        dimension wr(200),wi(200)
c        dimension jp1(5)
c        dimension dlancz(200)
c        dimension axfft(200),bxfft(200),arfft(200),aifft(200)
c        dimension sreallap(200),timelap(200),cflap(200)
        dimension wr(mst),wi(mst),jp1(5),dlancz(mst)
        dimension axfft(mst),bxfft(mst),arfft(mst),aifft(mst)
        dimension timelap(mst),cflap(mst)

c        dimension zlt1(20000),zlt2(20000),zlt3(20000)
        dimension zlt1(3*mst*mel),zlt2(3*mst*mel),zlt3(3*mst*mel)
c        dimension zld1(20000),zld2(20000)
        dimension zld1(3*mst*mel),zld2(3*mst*mel)
c        dimension slt1(200),slt2(200),slt3(200),sld1(200),sld2(200)
        dimension slt1(3*mst*mel)
c        dimension sfarr(400),sfarrzone(400),sfarrnonzone(400)
        dimension sfarr(4*mnd),sfarrzone(4*mnd),sfarrnonzone(4*mnd)
c        dimension zzlrar(160000)
        dimension zzlrar(16*mnd*mnd)
c          ,zzarr(310)
c        dimension zzloadz(400),zztempz(400)
        dimension zzloadz(4*mnd),zztempz(4*mnd)
c        dimension linv(400),minv(400)

        dimension zzarr(1+9*mel)
        dimension linv(2*mnd),minv(2*mnd)
c        dimension nplot(10)
        dimension nplot(mpl)
c        DIMENSION NHPOS(50)
c        DIMENSION ZLOAD1(300)
        DIMENSION NHPOS(mel)
        DIMENSION ZLOAD1(3*mnd)
c        dimension zloadf(200),ztempf(200)
        dimension zloadf(2*mnd),ztempf(2*mnd)
        dimension zloadf1(3*mnd)
c
c       for traping IEE floatingipoint exceptions
       iii=ieee_handler("set", "common", common_handler)

        OPEN (UNIT=1,FILE='BEMLTDAT',STATUS='UNKNOWN')
        OPEN (UNIT=2,FILE='BEMLTSTR',STATUS='UNKNOWN')
        OPEN (UNIT=3,FILE='BEMLTTRA',STATUS='UNKNOWN')
c        OPEN (UNIT=6,FILE='BEMLTOUT',STATUS='UNKNOWN')
cc       file to save the tractions of the free field, in case
cc       of subregion method of calculating them
c        OPEN (UNIT=4,FILE='/tmp/a131.7',STATUS='UNKNOWN')
cc       file to save the stresses of the free field, in case
cc       of subregion method of calculating them
c        OPEN (UNIT=5,FILE='/tmp/a132.7',STATUS='UNKNOWN')
c        OPEN (UNIT=7,FILE='/tmp/a133.7',STATUS='UNKNOWN')
c        OPEN (UNIT=8,FILE='/tmp/a134.7',STATUS='UNKNOWN')
c        OPEN (UNIT=9,FILE='/tmp/a135.7',STATUS='UNKNOWN')
c        OPEN (UNIT=10,FILE='/tmp/a136.7',STATUS='UNKNOWN')
c        OPEN (UNIT=11,FILE='/tmp/a137.7',STATUS='UNKNOWN')
c        OPEN (UNIT=12,FILE='/tmp/a138.7',STATUS='UNKNOWN')
cc       file to save displacements for postbem ploting
c        open (unit=11,file='BEMLTDES',status='new')
cc       file to save hoop stresses for postbem ploting
c        OPEN (UNIT=13,FILE='BEMLTHOO',STATUS='UNKNOWN')
cc       file to save displacements for postbem ploting
cc        open (unit=22,file='BEMLTDES',status='unknown')
cc       file to read plot informations saved in preprocessing
c        open (unit=23,file='POSTDAT',status='unknown')
C
c       of subregion method of calculating them
        OPEN (UNIT=4,STATUS='SCRATCH')
c       file to save the stresses of the free field, in case
c       of subregion method of calculating them
c        OPEN (UNIT=5,STATUS='SCRATCH')
c        OPEN (UNIT=17,STATUS='SCRATCH')
c       files 7,8: temporary files for saving the displacements
c       and traction kernels (not of inclusion)
        OPEN (UNIT=7,STATUS='SCRATCH')
        OPEN (UNIT=8,STATUS='SCRATCH')
        OPEN (UNIT=9,STATUS='SCRATCH')
c        OPEN (UNIT=10,STATUS='SCRATCH')
        OPEN (UNIT=11,STATUS='SCRATCH')
        OPEN (UNIT=12,STATUS='SCRATCH')
cc       file to save displacements for postbem ploting
c        open (unit=22,file='BEMLTDES',status='unknown')
c       file to read plot informations saved in preprocessing
        open (unit=23,file='POSTDAT',status='unknown')
c       file to save hoop stresses for postbem ploting
        OPEN (UNIT=13,FILE='BEMLTHOO',STATUS='UNKNOWN')
c       files 14,15: temporary files for saving the displacements
c       and traction kernels of inclusion
        OPEN (UNIT=14,STATUS='SCRATCH')
        OPEN (UNIT=15,STATUS='SCRATCH')
c       file for temporary saving of left hand side array
        OPEN (UNIT=16,STATUS='SCRATCH')
        OPEN (UNIT=18,STATUS='SCRATCH')

        OPEN (UNIT=19,STATUS='SCRATCH')
c        OPEN (UNIT=20,STATUS='SCRATCH')
c        OPEN (UNIT=21,STATUS='SCRATCH')
C
        REWIND 1
        REWIND 2
        REWIND 3
        REWIND 4
c        rewind 5
c        REWIND 6
        REWIND 7
        REWIND 8
        REWIND 9
c        REWIND 10
        REWIND 11
        REWIND 12
        REWIND 13
        rewind 14
        rewind 15
        rewind 16
c        rewind 17
        rewind 18
        rewind 19
c        rewind 20
c        rewind 21
c        rewind 22
        REWIND 23
c        WRITE(6,900)
c  900   FORMAT(2X,//2X,'2-D LAPLACE TRANSFORM BEM PROGRAM'//)
cc        IMAX=860
c        IMAX=1701
c        JMAX=210
c        LMAX=80500
c        mmax=201
c        mlmax=301
C
c         NP1 : number of nodes for calculating hoop stresses
c         nffs:
c           1: the incident field is reflected on the free surface, by
c              the laws of reflection of harmonic waves, and
c              the free field is calculated by adding the incident field
c              and the reflected field
c           2: the free field is calculated by deviding the region in 
c              subregions and taking the equilibrium conditions on the 
c              surface of the subsurface
c           3: simulation case

c        READ(1,*) NNODE,NELEM,NTINT,NP1
        READ(1,*) NNODE,NELEM,NTINT,NP1,ndoma,nffs
c       case that domain not Laplace, write a message and stop the program
        if(ndoma.ne.1) then
          write(6,142)ndoma
  142     FORMAT(1X,'THIS PROGRAM IS ONLY FOR LAPLACE DOMAIN NDOMA=',
     &         I3)
          go to 101
        endif
c        IF(NP1.LT.nnode) then 
        IF(NP1.LT.nelem) then 
          READ(1,*)(NHPOS(k),k=1,np1)
        endif
c       3 stresses per node
        NP=3*NP1
c       read the parameter and model of viscosity 
        read(1,*) fvisc,nmode
c
C       READ THE 3 NECESSERY NUMBERS
C       FOR FFT
c        READ(1,*)(jp1(k),k=1,3)
        READ(1,*)(jp1(k),k=1,5)
        N1=NNODE
        N2=2*NNODE
        N3=NELEM
c       quadratic elements
        N4=3
        N5=3*NNODE
        n6=4*nelem
        n7=np
        n8=5
        n9=ntint
C
C       DYNAMIC STORAGE ALLOCATION
C
c        I1=1
c        I2=I1+N3*N4
c        I3=I2+N1*N4
c        I4=I3+N3*N4
c        I5=I4+N1*n4
c        I6=I5+N1*n4
c        I7=I6+N1*n4
c        i8=i7+n2
c        i9=i8+n3
c        i10=i9+n1
        I1=1
        I2=I1+N3*N4
        I3=I2+N3*N4
        I4=I3+N1*n4
        I5=I4+N1*n4
        i6=i5+n2
        i7=i6+n3
        i8=i7+n1

c
        J1=1
        J2=J1+N1
        J3=J2+N1
        L1=1
c        L2=L1+N2*N2
c        L3=L2+N2*N2
c        L4=L3+N2
c        L5=L4+N2
        L2=L1+N2
        L3=L2+N2

c
        m1=1
        m2=m1+n6
c
        ml1=1
c        ml2=ml1+n1
c        ml3=ml2+n1
c        ml4=ml3+n1
        ml2=ml1+n4*n3
        ml3=ml2+n4*n3
        ml4=ml3+n4*n3
C
c        IF(I10.GT.IMAX) GO TO 10
c        IF(J3.GT.JMAX) GO TO 10
c        IF(L5.GT.LMAX) GO TO 10
c        if(m2.gt.mmax) go to 10
c        if(ml4.gt.mlmax) go to 10
c        CALL DATAS(RARR(J1),RARR(J2),IARR(I1),IARR(I2),iarr(i3),
c     &       iarr(i4),iarr(i5),iarr(i6),iarr(i7),iarr(i8),iarr(i9),
c     &       ntot,nplot,N1,N2,N3,N4,N9,ndzon,mpl)
        CALL DATAS(RARR(J1),RARR(J2),IARR(I1),IARR(I2),iarr(i3),
     &       iarr(i4),iarr(i5),iarr(i6),iarr(i7),
     &       ntot,nplot,N1,N2,N3,N4,N9,ndzon,mpl)
c                 iarr(i8),iarr(i9),
        CALL GAUSWE
        CALL TRIGO(wr,wi,nf,dlancz,n9)
        CALL LAPIN(RARR(J1),RARR(J2),IARR(I1),IARR(I2),iarr(i3),
     &    iarr(i4),iarr(i5),iarr(i6),iarr(i7),
     &    zzlrar,ZARR(L1),ZARR(L2),zzloadz,
     &    zztempz,ZLOAD1,zloadf,zloadf1,ztempf,ndzon,NHPOS,
     &    timelap,cflap,wr,wi,jp1,nf,dlancz,axfft,bxfft,arfft,aifft,
     &    zlt1,zlt2,zlt3,zld1,zld2,slt1,
     &    sfarr,sfarrzone,sfarrnonzone,
     &    zzarr(ml1),zzarr(ml2),zzarr(ml3),
     &    N1,N2,N3,N4,N5,n7,n8,n9,linv,minv,ntot,nplot)
C             sreallap,slt2,slt3,sld1,sld2,iarr(i8),iarr(i9),ZARR(L3),ZARR(L4),
        CLOSE (1)
        CLOSE (2)
        CLOSE (3)
        CLOSE (4)
c        close (5)
c
        CLOSE (6)
        CLOSE (7)
        CLOSE (8)
        CLOSE (9)
c        CLOSE (10)
        CLOSE (11)
        CLOSE (12)
        CLOSE (13)
        CLOSE (14)
        CLOSE (15)
        CLOSE (16)
c        close (17)
        close (18)
        close (19)
c        close (20)
c        close (21)
c        close (22)
        CLOSE (23)
        GO TO 101
c  10    CONTINUE
c        WRITE(6,100) I3,J3,L5,IMAX,JMAX,LMAX
c  100   FORMAT(1X,'INSUFFICIENT CORE SPACE ',6I5)
  101   STOP
        END

C -----------------------------------------------
c        SUBROUTINE DATAS(XCR,YCR,KONN,konzon,konnot,nonn,nonzon,
c     &    nonnot,MBCU,lbcu,ndtype,ntot,nplot,N1,N2,N3,N4,n9,
c     &    ndzon,mpl)
        SUBROUTINE DATAS(XCR,YCR,KONN,konzon,nonn,nonzon,
     &    MBCU,lbcu,ndtype,ntot,nplot,N1,N2,N3,N4,n9,
     &    ndzon,mpl)
C
C       READS IN DATA AND PRINTS IT OUT FOR CHECKING
C       LOADING WILL BE READ LATER
C       index for the existence of a free surface
c       nfre=0: no free surface exists
c       nfre=1: free surface exists (seminfinity)
C
        IMPLICIT REAL*8 (A-H,O-Y)
        DIMENSION XCR(N1),YCR(N1)
        DIMENSION KONN(N3,N4),MBCU(N2)
        DIMENSION konzon(N3,N4)
        dimension nonzon(n1,n4)
c          ,konnot(N3,N4),nonnot(n1,n4)
        dimension nonn(n1,n4)
        dimension lbcu(n3),ndtype(n1)
c        dimension nplot(10)
        dimension nplot(mpl)
        character fnamepro*14,fnamepos*14
c        COMMON/CONTR/NNODE,NELEM,NTINT,NP
        COMMON/CONTR/NNODE,NELEM,NP
c                        ##########################
c                        !!!!!!!!!!!!!!!!!!!!!!!!!!
c                 modification for viscoelesticity
        COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM,fvisc,nmode
        common/free/nffs,nfre
C
        READ(1,*) RLA,RMU,RNU,RHO,C1,C2,TIM,nfre
        READ(1,*) XCR
        READ(1,*) YCR
        READ(1,*) KONN
        READ(1,*) nonn
        READ(1,*) MBCU
        READ(1,*) lbcu
        read(1,*) ndtype
        READ(1,*) konzon
c        READ(1,*) konnot
        READ(1,*) nonzon
c        READ(1,*) nonnot
c
c        read(23,*)nplot
        read(23,*)ntot,nplot,fnamepro,fnamepos
        OPEN (UNIT=6,FILE=fnamepro,STATUS='UNKNOWN')
        REWIND 6
c       find the number of nodes for the subregion (simulation case)
c       ATTENTION: for future versions must take and the nodes of zones
        do 68 i=1,n1
          if(ndtype(i).eq.3) ndzon=ndzon+1
  68    continue
C
C       PRINTING OF INPUT
C
C        WRITE(6,100)NNODE,NELEM,NTINT,NP,RLA,RMU,RNU,RHO,C1,C2,TIM
        WRITE(6,100)NNODE,NELEM,n9,NP,RLA,RMU,RNU,RHO,C1,C2,TIM
  100   FORMAT(2X,' NNODE=',I3,2X,' NELEM=',I3,2X,' NTINT=',I3,2X,
     &      ' NP=',I3//2X,' RLA=',E12.6,2X,' RMU=',E12.6,2X,
     &      ' RNU=',E12.6,2X,' RHO=',E12.6//2X,
     &      ' C1=',E12.6,2X,' C2=',E12.6,2X,' TIM=',E12.6//)
        WRITE(6,94)
  94    FORMAT(2X,' X - COORDINATES'/)
        WRITE(6,101)(XCR(I),I=1,NNODE)
        WRITE(6,95)
  95    FORMAT(//2X,' Y - COORDINATES'/)
        WRITE(6,101,err=81)(YCR(I),I=1,NNODE)
  81    WRITE(6,96)
  96    FORMAT(//2X,' CONNECTIVITY'/)
        WRITE(6,97,err=82)((KONN(I,J),J=1,3),I=1,NELEM)
  97    FORMAT(2X,3I5)
        write(6,58)
  58    format(//2x,' BOUNDARY CONDITIONS')
  82    WRITE(6,98,err=83) (MBCU(I),I=1,n2)
c  98    FORMAT(//2X,' BOUNDARY CONDITIONS'/2X,15I4)
  98    FORMAT(25(i3))
        write(6,59)
  59    format(//2x,' POSITION OF THE ELEMENT ')
  83    WRITE(6,103,err=84) (lbcu(I),I=1,n3)
c  103   FORMAT(//2X,' POSITION OF THE ELEMENT '/2X,15I4)
  103   format(25(i3))
        write(6,60)
  60    format(//2X,' TYPE OF THE NODE ')
  84    WRITE(6,104,err=85) (ndtype(I),I=1,n1)
c  104   FORMAT(//2X,' TYPE OF THE NODE '/2X,15I4)
  104   format(25(i3))
C
  101   FORMAT(8(2X,E12.6))
  85    RETURN
        END
C -----------------------------------------------
C
        SUBROUTINE GAUSWE
C
C       GAUSSIAN QUDRATURE POINTS AND WEIGHTS
C
        IMPLICIT REAL*8 (A-H,O-Y)
        COMMON/GAUSS/GPN(10,10),WTN(10,10),GPS(10,10),WTS(10,10)
C
C       NORMAL - GAUSS-LEGENDRE
C
        CALL CLRARR(GPN,100)
        CALL CLRARR(WTN,100)
        CALL CLRARR(GPS,100)
        CALL CLRARR(WTS,100)
C
        GPN(2,1)=-0.5773502692D0
        WTN(2,1)= 1.0000000000D0
        GPN(2,2)=+0.5773502692D0
        WTN(2,2)= 1.0000000000D0
C
        GPN(3,1)=-0.7745966692D0
        WTN(3,1)= 0.5555555556D0
        GPN(3,2)=+0.0000000000D0
        WTN(3,2)= 0.8888888889D0
        GPN(3,3)=+0.7745966692D0
        WTN(3,3)= 0.5555555556D0
C
        GPN(4,1)=-0.8611363116D0
        WTN(4,1)= 0.3478548451D0
        GPN(4,2)=-0.3399810436D0
        WTN(4,2)= 0.6521451548D0
        GPN(4,3)= 0.3399810436D0
        WTN(4,3)= 0.6521451548D0
        GPN(4,4)= 0.8611363116D0
        WTN(4,4)= 0.3478548451D0
C
        GPN(5,1)=-0.9061798459D0
        WTN(5,1)= 0.2369268850D0
        GPN(5,2)=-0.5384693101D0
        WTN(5,2)= 0.4786286704D0
        GPN(5,3)= 0.0D0
        WTN(5,3)= 0.5688888889D0
        GPN(5,4)= 0.5384693101D0
        WTN(5,4)= 0.4786286704D0
        GPN(5,5)= 0.9061798459D0
        WTN(5,5)= 0.2369268850D0
C
        GPN(6,1)=-0.9324695142D0
        WTN(6,1)= 0.1713244924D0
        GPN(6,2)=-0.6612093865D0
        WTN(6,2)= 0.3607615730D0
        GPN(6,3)=-0.2386191861D0
        WTN(6,3)= 0.4679139346D0
        GPN(6,4)= 0.2386191861D0
        WTN(6,4)= 0.4679139346D0
        GPN(6,5)= 0.6612093865D0
        WTN(6,5)= 0.3607615730D0
        GPN(6,6)= 0.9324695142D0
        WTN(6,6)= 0.1713244924D0
C
        GPN(7,1)=-0.9491079123D0
        WTN(7,1)= 0.1294849662D0
        GPN(7,2)=-0.7415311856D0
        WTN(7,2)= 0.2797053915D0
        GPN(7,3)=-0.4058451514D0
        WTN(7,3)= 0.3818300505D0
        GPN(7,4)= 0.0D0
        WTN(7,4)= 0.4179591837D0
        GPN(7,5)= 0.4058451514D0
        WTN(7,5)= 0.3818300505D0
        GPN(7,6)= 0.7415311856D0
        WTN(7,6)= 0.2797053915D0
        GPN(7,7)= 0.9491079123D0
        WTN(7,7)= 0.1294849662D0
C
        GPN(8,1)=-0.9602898565D0
        WTN(8,1)= 0.1012285363D0
        GPN(8,2)=-0.7966664774D0
        WTN(8,2)= 0.2223810344D0
        GPN(8,3)=-0.5255324099D0
        WTN(8,3)= 0.3137066459D0
        GPN(8,4)=-0.1834346425D0
        WTN(8,4)= 0.3626837834D0
        GPN(8,5)= 0.1834346425D0
        WTN(8,5)= 0.3626837834D0
        GPN(8,6)= 0.5255324099D0
        WTN(8,6)= 0.3137066459D0
        GPN(8,7)= 0.7966664774D0
        WTN(8,7)= 0.2223810344D0
        GPN(8,8)= 0.9602898565D0
        WTN(8,8)= 0.1012285363D0
C
        GPN(9,1)=-0.9681602395D0
        WTN(9,1)= 0.0812743884D0
        GPN(9,2)=-0.8360311073D0
        WTN(9,2)= 0.1806481607D0
        GPN(9,3)=-0.6133714327D0
        WTN(9,3)= 0.2606106964D0
        GPN(9,4)=-0.3242534234D0
        WTN(9,4)= 0.3123470770D0
        GPN(9,5)= 0.0D0
        WTN(9,5)= 0.3302393550D0
        GPN(9,6)= 0.3242534234D0
        WTN(9,6)= 0.3123470770D0
        GPN(9,7)= 0.6133714327D0
        WTN(9,7)= 0.2606106964D0
        GPN(9,8)= 0.8360311073D0
        WTN(9,8)= 0.1806481607D0
        GPN(9,9)= 0.9681602395D0
        WTN(9,9)= 0.0812743884D0
C
        GPN(10,1)=-0.9739065285D0
        WTN(10,1)= 0.0666713443D0
        GPN(10,2)=-0.8650633667D0
        WTN(10,2)= 0.1494513492D0
        GPN(10,3)=-0.6794095683D0
        WTN(10,3)= 0.2190863625D0
        GPN(10,4)=-0.4333953941D0
        WTN(10,4)= 0.2692667193D0
        GPN(10,5)=-0.1488743390D0
        WTN(10,5)= 0.2955242247D0
        GPN(10,6)= 0.1488743390D0
        WTN(10,6)= 0.2955242247D0
        GPN(10,7)= 0.4333953941D0
        WTN(10,7)= 0.2692667193D0
        GPN(10,8)= 0.6794095683D0
        WTN(10,8)= 0.2190863625D0
        GPN(10,9)= 0.8650633667D0
        WTN(10,9)= 0.1494513492D0
        GPN(10,10)= 0.9739065285D0
        WTN(10,10)= 0.0666713443D0
C
C       SINGULAR - LOG WEIGHTED
C
C
        GPS(2,1)= 0.1120088061D0
        GPS(2,2)= 0.6022769081D0
        WTS(2,1)= 0.7185393190D0
        WTS(2,2)= 0.2814606809D0
C
        GPS(3,1)= 0.06389079308D0
        GPS(3,2)= 0.3689970637D0
        GPS(3,3)= 0.7668803039D0
        WTS(3,1)= 0.5134045522D0
        WTS(3,2)= 0.3919800412D0
        WTS(3,3)= 0.09461540656D0
C
        GPS(4,1)= 0.04144848019D0
        GPS(4,2)= 0.2452749143D0
        GPS(4,3)= 0.5561654535D0
        GPS(4,4)= 0.8489823945D0
        WTS(4,1)= 0.3834640681D0
        WTS(4,2)= 0.3868753177D0
        WTS(4,3)= 0.1904351269D0
        WTS(4,4)= 0.03922548712D0
C
        GPS(5,1)= 0.02913447215D0
        GPS(5,2)= 0.1739772133D0
        GPS(5,3)= 0.4117025202D0
        GPS(5,4)= 0.6773141745D0
        GPS(5,5)= 0.8947713610D0
        WTS(5,1)= 0.2978934717D0
        WTS(5,2)= 0.3497762265D0
        WTS(5,3)= 0.2344882900D0
        WTS(5,4)= 0.09893045951D0
        WTS(5,5)= 0.01891155214D0
C
        GPS(6,1)= 0.02163400584D0
        GPS(6,2)= 0.1295833911D0
        GPS(6,3)= 0.3140204499D0
        GPS(6,4)= 0.5386572173D0
        GPS(6,5)= 0.7569153373D0
        GPS(6,6)= 0.9226688513D0
        WTS(6,1)= 0.2387636625D0
        WTS(6,2)= 0.3082865732D0
        WTS(6,3)= 0.2453174265D0
        WTS(6,4)= 0.1420087565D0
        WTS(6,5)= 0.05545462232D0
        WTS(6,6)= 0.01016895869D0
C
        GPS(7,1)= 0.01671935540D0
        GPS(7,2)= 0.1001856779D0
        GPS(7,3)= 0.2462942462D0
        GPS(7,4)= 0.4334634932D0
        GPS(7,5)= 0.6323509880D0
        GPS(7,6)= 0.8111186267D0
        GPS(7,7)= 0.9408481667D0
        WTS(7,1)= 0.1961693894D0
        WTS(7,2)= 0.2703026442D0
        WTS(7,3)= 0.2396818730D0
        WTS(7,4)= 0.1657757748D0
        WTS(7,5)= 0.08894322713D0
        WTS(7,6)= 0.03319430435D0
        WTS(7,7)= 0.005932787015D0
C
        GPS(8,1)= 0.01332024416D0
        GPS(8,2)= 0.07975042901D0
        GPS(8,3)= 0.1978710293D0
        GPS(8,4)= 0.3541539943D0
        GPS(8,5)= 0.5294585752D0
        GPS(8,6)= 0.7018145299D0
        GPS(8,7)= 0.8493793204D0
        GPS(8,8)= 0.9533264500D0
        WTS(8,1)= 0.1644166047D0
        WTS(8,2)= 0.2375256100D0
        WTS(8,3)= 0.2268419844D0
        WTS(8,4)= 0.1757540790D0
        WTS(8,5)= 0.1129240302D0
        WTS(8,6)= 0.05787221071D0
        WTS(8,7)= 0.02097907374D0
        WTS(8,8)= 0.003686407104D0
        RETURN
        END
C -----------------------------------------------
        SUBROUTINE CLRZRR(ZA,NA)
C
C       BLANKS COMPLEX ARRAY
C
        IMPLICIT COMPLEX*16 (Z)
c        DIMENSION ZA(1)
        DIMENSION ZA(na)
        DO 1 I=1,NA
        ZA(I)=DCMPLX(0.0D0,0.0D0)
  1     CONTINUE
        RETURN
        END
C -----------------------------------------------
        SUBROUTINE CLRARR(AA,NA)
C
C       BLANKS REAL ARRAY
C
        IMPLICIT REAL*8 (A)
c        DIMENSION AA(1)
        dimension aa(na)
        DO 1 I=1,NA
        AA(I)=0.0D0
  1     CONTINUE
        RETURN
        END
cC -----------------------------------------------
c        SUBROUTINE PLOTGR(A,B,NRDM,NR)
cC
cC       PLOTS B VERSUS A AND TABULATES
cC
c        IMPLICIT REAL*8 (A-H,O-Y)
c        CHARACTER*1 BLANK,DOT,EX,LINE(61)
c        DIMENSION A(NRDM),B(NRDM)
c        DATA BLANK,DOT,EX/ ' ','.','X' /
cC
c        WRITE(6,100)
c        DO 1 I=1,61
c  1     LINE(I)=DOT
c        WRITE(6,101) LINE
c        BMAX=DABS(B(1))
c        DO 2 I=1,NR
c        BMIN=DABS(B(I))
c        IF (BMIN.GT.BMAX) BMAX=BMIN
c  2     CONTINUE
cC               MODIFICATION  ERROR OF THE PROGRAM   !!!!!!!
cC       IF (BMAX.EQ.0.) BAX=1.
c        IF (BMAX .LT. 1.0D-15) BMAX=1.0D0
cC               END MODIF   !!!!!!!!!
c        DO 3 I=1,NR
c        AA=A(I)
c        BB=B(I)
c        DO 4 J=1,61
c  4     LINE(J)=BLANK
c        LINE(31)=DOT
cc        JUG=30.0D0*(BB/BMAX+1.0D0)+1.5D0
c        JUG=30.0*(BB/BMAX+1.0)+1.5
c        LINE(JUG)=EX
c        WRITE(6,102) LINE,AA,BB
c  3     CONTINUE
c  100   FORMAT(//,74X,'VARIABLE A',10X,'VARIABLE B')
c  101   FORMAT(4X,61A1)
c  102   FORMAT(4X,61A1,2(5X,E15.6))
c        RETURN
c        END
C------------------------------------------------
C
        SUBROUTINE TRIGO(wr,wi,nf,dlancz,n9)
C
C       NTINT ROOTS OF UNITY, LANCZOS CONVERGENCE FACTORS
C
C       NTINT IS NUMBER OF TIME AND/OR LAPLACE SPECTRUM POINTS
C       (IT CAN BE SET TO EITHER 20 OR 50)
C
        IMPLICIT REAL*8 (A-H,O-Y)
c                         ##############################
c                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c       modification for 200 time steps
c        COMMON/TRIG/WR(50),WI(50),JP1(5),JC1(5),NF
c        COMMON/LANC/DLANCZ(50)
c        COMMON/TRIG/WR(200),WI(200),JP1(5),JC1(5),NF
        dimension wr(n9),wi(n9)
c         COMMON/LANC/DLANCZ(200)
        dimension dlancz(n9)
C
C             MODIF !!!!!!!!!!!!!!
c        CALL CLRARR(WR,50)
c        CALL CLRARR(WI,50)
c        CALL CLRARR(DLANCZ,50)
C             AND MODIF  !!!!!!
c        N=NTINT
        N=n9
c           attention       ##############################
c                           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c        IF (N.EQ.20) GO TO 5
c        JP1(1)=2
c        JP1(2)=5
c        JP1(3)=5
c        GO TO 15
c  5     JP1(1)=2
c        JP1(2)=2
c        JP1(3)=5
  15    NF=3
c
c
C
C       N ROOTS OF UNITY
C
        PI=3.141592653598979323D0
        WR(1)=1.0D0
        WI(1)=0.0D0
        COEF=2.0D0*PI
        MN=N/2+1
        DO 25 I=2,MN
        REX=DFLOAT(I-1)
        ANG=COEF*REX/DFLOAT(N)
        WR(I)=DCOS(ANG)
        WI(I)=DSIN(ANG)
  25    CONTINUE
        L=MN+1
        DO 45 I=L,N
        J=N-I+2
        WR(I)=WR(J)
        WI(I)=-WI(J)
  45    CONTINUE
C
C       LANCZOS FACTORS
C
        DN=DFLOAT(N)
        DPIN=3.14159265358979323D+00/DN
        DLANCZ(1)=1.0D+00
        DO 65 K=2,N
        DS=DPIN*DFLOAT(K-1)
        L=-1
        DTERM=1.0D+00
        DSOM=1.0D+00
  35    L=L+1
        DNOM=DFLOAT((2*L+2)*(2*L+3))
        DTERM=-DTERM*DS*DS/DNOM
        DSOM=DSOM+DTERM
        IF (DABS(DTERM) .GE. 1.0D-16) GO TO 35
        DLANCZ(K)=DSOM
  65    CONTINUE
        RETURN
        END
C------------------------------------------------
        SUBROUTINE LAPIN(XCR,YCR,KONN,konzon,nonn,nonzon,
     &    MBCU,lbcu,ndtype,zlrar,ZLOAD,ZTEMP,
     &    zloadz,ztempz,ZLOAD1,zloadf,zloadf1,ztempf,ndzon,NHPOS,
     &    time,cf,wr,wi,jp1,nf,dlancz,ax,bx,arfft,aifft,
     &    zlt1,zlt2,zlt3,zld1,zld2,slt1,
     &    sfdiag,sfdiagzone,SFDIAGnonzone,
     &    zz1,zz2,zz3,N1,N2,N3,N4,N5,n7,
     &    n8,n9,linv,minv,ntot,nplot)

c            sreal,slt2,slt3,sld1,sld2,konnot,nonnot,ZGM,ZFM,
c
C       PERFORMS NUMERICAL INVERSION OF LAPLACE TRANSFORM DOMAIN
C
C       ZLT CONTAINS NTINT POINTS OF NP FUNCTIONS IN LAPLACE
C       TRANSFORM DOMAIN THAT WILL BE INVERTED
C       SR CONTAINS THE REAL PART OF THESE NP FUNCTIONS AFTER
C       THEY HAVE BEEN INVERTED TO THE TIME DOMAIN, ALSO,
C       SR AGAIN CONTAINS NTINT POINTS
C
        IMPLICIT COMPLEX*16 (Z)
        IMPLICIT REAL*8 (A-H,O-Y)
c
        DIMENSION XCR(N1),YCR(N1)
c        DIMENSION ZGM(N2,N2),ZFM(N2,N2)
        dimension ZLOAD(N2),ZTEMP(N2)
        DIMENSION KONN(N3,N4),MBCU(N2)
        DIMENSION konzon(N3,N4)
c            ,konnot(N3,N4)
        dimension nonn(n1,n4)
        dimension nonzon(n1,n4)
c            ,nonnot(n1,n4)
        dimension lbcu(n3),ndtype(n1)
        dimension time(n9),cf(n9)
c         sreal(n9),
        dimension wr(n9),wi(n9),jp1(n8)
        dimension dlancz(n9)
        dimension arfft(n9),aifft(n9)
        dimension ax(n9),bx(n9)
c        dimension zlt1(n9,n1),zlt2(n9,n1),zlt3(n9,n1)
        dimension zlt1(n9,n4,n3),zlt2(n9,n4,n3),zlt3(n9,n4,n3)
c        dimension zld1(n9,n1),zld2(n9,n1)
        dimension zld1(n9,n4,n3),zld2(n9,n4,n3)
c        dimension slt1(n9),slt2(n9),slt3(n9),sld1(n9),sld2(n9)
        dimension slt1(n9,n4,n3)
        dimension sfdiag(n1,2,2)
        dimension zlrar(n2+2*ndzon,n2+2*ndzon)
        dimension zloadz(n2+2*ndzon),ztempz(n2+2*ndzon)
        dimension sfdiagnonzone(n1,2,2)
c        dimension sfdiagzone(ndzon,2,2)
        dimension sfdiagzone(n1,2,2)

        dimension linv(n2+2*ndzon),minv(n2+2*ndzon)
        dimension zloadf(n2),ztempf(n2)
        dimension zloadf1(n5)

c        dimension zz1(n1),zz2(n1),zz3(n1)
        dimension zz1(n4,n3),zz2(n4,n3),zz3(n4,n3)
        DIMENSION NHPOS(N3)
        DIMENSION ZLOAD1(N5)
c        dimension nplot(10)
        dimension nplot(ntot)
        COMMON/CONTR/NNODE,NELEM,NP
c                 modification for viscoelesticity
        COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM,fvisc,nmode
c
        common/free/nffs,nfre
c        N=NTINT
        N=n9
C
c       in case that there is free surface, in the subregion case (nffs=2)
c       loop for the 2 steps a,b
c       in the case that there is simulation surface (nffs=3), loop
c       for the 3 steps a,b,c
c         a. ic=1, calculation of the free field of the 
c             i)  non free surface boundary, or
c             ii) simulation boundary
c         b. ic=2, with this free field, calculate the 
c             i)  scattered field on non free surface boundary
c             ii) non simulation boundary (free field)
cc         c. ic=3, for the case only without free surface, with loads
cc              the free fields of the previous 2 steps, calculate the
cc              final scattered field

        ic=1
        if(nffs.eq.2) ic=2
        if(nffs.eq.3) ic=2
c
c          call clrzrr(zlt1,n9*n1)
c          call clrzrr(zlt2,n9*n1)
c          call clrzrr(zlt3,n9*n1)
c          call clrzrr(zld1,n9*n1)
c          call clrzrr(zld2,n9*n1)
          call clrzrr(zlt1,n9*n4*n3)
          call clrzrr(zlt2,n9*n4*n3)
          call clrzrr(zlt3,n9*n4*n3)
          call clrzrr(zld1,n9*n4*n3)
          call clrzrr(zld2,n9*n4*n3)

          PI=3.14159265359D0
          OMEGA=2.0D0*PI/TIM
          AK=6.0D0/TIM
          DT=TIM/DFLOAT(N)
          DO 2 K=1,N
          TIME1=DT*DFLOAT(K-1)
          CF(K)=2.0D0*DEXP(AK*TIME1)/TIM
          TIME(K)=TIME1
          AKMG=OMEGA*DFLOAT(K-1)
          IF (k .EQ. 1) CALL CLRARR(SFDIAG,n1*4)
          IF (k .EQ. 1) CALL CLRARR(SFDIAGzone,ndzon*4)
          IF (k .EQ. 1) CALL CLRARR(SFDIAGnonzone,n1*4)
c          CALL CLRARR(SFDIAG,n1*4)
c          CALL CLRARR(SFDIAGzone,ndzon*4)
c          CALL CLRARR(SFDIAGnonzone,n1*4)
        do 41 icffs=1,ic
c          rewind 2
c          rewind 3
          rewind 4
c          rewind 5
          rewind 17
          rewind 7
          rewind 8
c          rewind 10
          rewind 11
          rewind 12
          rewind 14
          rewind 15
c          call clrzrr(zz1,n1)
c          call clrzrr(zz2,n1)
c          call clrzrr(zz3,n1)
          call clrzrr(zz1,n4*n3)
          call clrzrr(zz2,n4*n3)
          call clrzrr(zz3,n4*n3)
          call clrzrr(zload1,n5)
          call clrzrr(zloadf1,n5)
          call clrzrr(zload,n2)
          call clrzrr(zloadf,n2)
          call clrzrr(ztemp,n2)
          call clrzrr(ztempf,n2)

          CALL BEM(AK,AKMG,XCR,YCR,KONN,konzon,nonn,nonzon,
     &      MBCU,lbcu,ndtype,zlrar,ZLOAD,ZTEMP,zloadz,
     &      ztempz,ZLOAD1,zloadf,zloadf1,ztempf,ndzon,K,NHPOS,zlt1,zlt2,
     &      zlt3,zld1,zld2,sfdiag,sfdiagzone,SFDIAGnonzone,
     &      zz1,zz2,zz3,N1,N2,N3,N4,N5,n7,n9,icffs,linv,minv)
c                 konnot,nonnot,ZGM,ZFM,
  41      continue
  2     CONTINUE
c ************** restriction 10 plots ***************
c       calculate and save in file 13, the hoop stresses and x,y global 
c       displacements for ploting with postbem
        np1=np/3
        kdiak1=0
        kdiak2=0
        kdiak3=0
        kdiak4=0
        kdiak5=0
c        do 20 l=1,10
        do 20 l=1,ntot
c         case of hoop stresses plots  
c         normal stress in radial direction system (nplot=1)
c         or in Cartesian system (nplot=11)
          if(((nplot(l).eq.1).or.(nplot(l).eq.11)).and.(kdiak1.ne.1))
     &    then
            kdiak1=1
c           Invertion of the Laplace transform to obtain the hoop stresses
            DO 1 I=1,NP1
c             loop for the 3 nodes of the quadratic element
              do 1 j=1,n4
                DO 3 K=1,N
c                  ZZ=ZLT1(K,i)
                  ZZ=ZLT1(K,j,i)
                  ZZ=ZZ*DLANCZ(K)
                  AX(K)=DREAL(ZZ)
                  BX(K)=DIMAG(ZZ)
  3             CONTINUE
                AX(1)=0.5D0*AX(1)
                CALL FFCT(AX,BX,N,wr,wi,jp1,nf,arfft,aifft,n8)
                AX(1)=2.0D0*AX(1)
                DO 6 K=1,N
c                  slt1(K,I)=AX(K)
c                  slt1(K,I)=CF(K)*slt1(K,I)
                  slt1(K,j,I)=AX(K)
                  slt1(K,j,I)=CF(K)*slt1(K,j,I)
c                  slt1(k)=AX(K)
c     4             CONTINUE
c                   DO 6 K=1,N
c                  slt1(k)=CF(K)*slt1(k)
  6             CONTINUE
c                WRITE(13,*) slt1
  1         CONTINUE
            WRITE(13,*) slt1
          endif
c         tangential stress in radial direction system (nplot=2)
c         or in Cartesian system (nplot=12)
          if(((nplot(l).eq.2).or.(nplot(l).eq.12)).and.(kdiak2.ne.1))
     &    then
            kdiak2=1
c           Invertion of the Laplace transform to obtain the hoop stresses
            DO 11 I=1,NP1
c             loop for the 3 nodes of the quadratic element
              do 11 j=1,n4
              DO 31 K=1,N
c                ZZ=ZLT2(K,i)
                ZZ=ZLT2(K,j,i)
                ZZ=ZZ*DLANCZ(K)
                AX(K)=DREAL(ZZ)
                BX(K)=DIMAG(ZZ)
  31          CONTINUE
              AX(1)=0.5D0*AX(1)
              CALL FFCT(AX,BX,N,wr,wi,jp1,nf,arfft,aifft,n8)
              AX(1)=2.0D0*AX(1)
              DO 61 K=1,N
c                slt1(K,I)=AX(K)
c                slt1(K,I)=CF(K)*slt1(K,I)
                slt1(K,j,I)=AX(K)
                slt1(K,j,I)=CF(K)*slt1(K,j,I)
c                slt2(k)=AX(K)
c     41          CONTINUE
c                 DO 61 K=1,N
c                slt2(k)=CF(K)*slt2(k)
  61          CONTINUE
c                WRITE(13,*) slt2
  11        CONTINUE
            WRITE(13,*) slt1
          endif
c         shear stress in radial direction system (nplot=3)
c         or in Cartesian system (nplot=13)
          if(((nplot(l).eq.3).or.(nplot(l).eq.13)).and.(kdiak3.ne.1))
     &    then
            kdiak3=1
c           Invertion of the Laplace transform to obtain the hoop stresses
            DO 12 I=1,NP1
c             loop for the 3 nodes of the quadratic element
              do 12 j=1,n4
              DO 32 K=1,N
c                ZZ=ZLT3(K,i)
                ZZ=ZLT3(K,j,i)
                ZZ=ZZ*DLANCZ(K)
                AX(K)=DREAL(ZZ)
                BX(K)=DIMAG(ZZ)
  32          CONTINUE
              AX(1)=0.5D0*AX(1)
              CALL FFCT(AX,BX,N,wr,wi,jp1,nf,arfft,aifft,n8)
              AX(1)=2.0D0*AX(1)
              DO 62 K=1,N
c                slt1(K,I)=AX(K)
c                slt1(K,I)=CF(K)*slt1(K,I)
                slt1(K,j,I)=AX(K)
                slt1(K,j,I)=CF(K)*slt1(K,j,I)
c                slt3(k)=AX(K)
c      42          CONTINUE
c                  DO 62 K=1,N
c                slt3(k)=CF(K)*slt3(k)
  62          CONTINUE
c              WRITE(13,*) slt3
  12        CONTINUE
            WRITE(13,*) slt1
          endif
c         case of deformed shape plot(nplot=5) in radial direction,
c         or x displacements plot in Cartesian system(nplot=21)
          if(((nplot(l).eq.5).or.(nplot(l).eq.21)).and.(kdiak4.ne.1))
     &    then
            kdiak4=1
c           Invertion of the Laplace transform to obtain the displacements
c           x displacements
c            DO 13 I=1,n1
            DO 13 I=1,n3


c             loop for the 3 nodes of the quadratic element
              do 13 j=1,n4
              DO 33 K=1,N
c                ZZ=ZLd1(K,I)
                ZZ=ZLd1(K,j,I)
                ZZ=ZZ*DLANCZ(K)
                AX(K)=DREAL(ZZ)
                BX(K)=DIMAG(ZZ)
  33          CONTINUE
              AX(1)=0.5D0*AX(1)
              CALL FFCT(AX,BX,N,wr,wi,jp1,nf,arfft,aifft,n8)
              AX(1)=2.0D0*AX(1)
              DO 63 K=1,N
c                slt1(K,I)=AX(K)
c                slt1(K,I)=CF(K)*slt1(K,I)
                slt1(K,j,I)=AX(K)
                slt1(K,j,I)=CF(K)*slt1(K,j,I)
c                sld1(k)=AX(K)
c    43          CONTINUE
c                DO 63 K=1,N
c                sld1(k)=CF(K)*sld1(k)
  63          CONTINUE
c              WRITE(13,*) sld1
  13        CONTINUE
            WRITE(13,*) slt1
          endif
c         case of deformed shape plot(nplot=5) in radial direction,
c         or y displacements plot in Cartesian system(nplot=22)
          if(((nplot(l).eq.5).or.(nplot(l).eq.22)).and.(kdiak5.ne.1))
     &    then
            kdiak5=1
c           Invertion of the Laplace transform to obtain the displacements
c           y displacements
c            DO 14 I=1,n1
            DO 14 I=1,n3

c             loop for the 3 nodes of the quadratic element
              do 14 j=1,n4
              DO 34 K=1,N
c                ZZ=ZLd2(K,I)
                ZZ=ZLd2(K,j,I)
                ZZ=ZZ*DLANCZ(K)
                AX(K)=DREAL(ZZ)
                BX(K)=DIMAG(ZZ)
  34          CONTINUE
              AX(1)=0.5D0*AX(1)
              CALL FFCT(AX,BX,N,wr,wi,jp1,nf,arfft,aifft,n8)
              AX(1)=2.0D0*AX(1)
              DO 64 K=1,N
c                slt1(K,I)=AX(K)
c                slt1(K,I)=CF(K)*slt1(K,I)
                slt1(K,j,I)=AX(K)
                slt1(K,j,I)=CF(K)*slt1(K,j,I)
c                sld2(k)=AX(K)
c    44          CONTINUE
c                DO 64 K=1,N
c                sld2(k)=CF(K)*sld2(k)
  64          CONTINUE
c              WRITE(13,*) sld2
  14        CONTINUE
            WRITE(13,*) slt1
          endif
  20    continue
cC
cC       PRINTING
cC
c        DO 11 I=1,NP
c        DO 10 K=1,N
c        SREAL(K)=SR(K,I)
cc        SREAL(K)=SR(K)
c  10    CONTINUE
c        WRITE(6,100) I
c  100   FORMAT(///10X,'TRANSFORM INVERSION NO. ',I3,///)
c        CALL PLOTGR(TIME,SREAL,N,N)
c  11    CONTINUE
        RETURN
        END
c
C------------------------------------------------
        SUBROUTINE FFCT(VR,VI,NT,wr,wi,jp1,nf,ar,ai,n8)
C
C       FAST FOURIER ALGORITHM FOR INVERSE TRANSFORMATION
C
C       VR, VI ARE REAL AND IMAGINARY PARTS OF FUNCTION TO BE INVERTED
C       NT IS NUMBER OF LAPLACE SPECTRUM POINTS
C
        IMPLICIT REAL*8 (A-H,O-Y)
c        DIMENSION VR(1),VI(1)
        dimension vr(nt),vi(nt)
c        DIMENSION AR(50),AI(50),AUXR(5),AUXI(5)
c        DIMENSION AR(200),AI(200),
        dimension AUXR(5),AUXI(5)
c        COMMON/TRIG/WR(50),WI(50),JP1(5),JC1(5),NF
c        COMMON/TRIG/WR(200),WI(200),JP1(5),JC1(5),NF
        dimension ar(nt),ai(nt)
        dimension wr(nt),wi(nt),jp1(n8),jc1(5)
C
        N=NT
        K=1
        J=0
        GO TO 8
  4     LP=1
        LS=0
        NM=NF-1
        DO 9 KK=1,NM
        LC=JC1(KK+1)-1
        LP=LP*JP1(KK)
        LS=LS+LC*LP
  9     CONTINUE
        LS=LS+JC1(1)
        AR(K)=VR(LS)
        AI(K)=VI(LS)
        K=K+1
  7     JC1(J)=JC1(J)+1
        NST=JC1(J)-JP1(J)
        IF (NST) 8,8,6
  8     JM=J+1
        IF (JM-NF) 2,2,1
  2     DO 3 L=JM,NF
        JC1(L)=1
  3     CONTINUE
  1     J=NF
        GO TO 4
  6     J=J-1
        IF (J) 5,5,7
  5     DO 11 I=1,N
        VR(I)=AR(I)
        VI(I)=AI(I)
  11    CONTINUE
        NUR=N
        DO 25 I=1,NF
        K=NF-I+1
        NU=JP1(K)
        NUR=NUR/NU
        NUD=N/(NUR*NU)
        LRD=NUR*NUD
        LDU=NUD*NU
        LS=(NUR-1)*LDU+1
        DO 24 KRE=1,LS,LDU
        KRE1=KRE-1
        DO 22 LE=1,NUD
        LE1=LE-1
        JAC=NUR*LE1
        DO 21 NS=1,NU
        AUXR(NS)=0.0D0
        AUXI(NS)=0.0D0
        LA=LE+KRE1
        DO 20 NUT=1,NU
        M=JAC*(NUT-1)
        MQ=M/N
        M=M+1-MQ*N
        Q=WI(M)
        AUXR(NS)=AUXR(NS)+VR(LA)*WR(M)-VI(LA)*Q
        AUXI(NS)=AUXI(NS)+VR(LA)*Q+VI(LA)*WR(M)
        LA=LA+NUD
  20    CONTINUE
        JAC=JAC+LRD
  21    CONTINUE
        LA=LE+KRE1
        DO 23 NS=1,NU
        VR(LA)=AUXR(NS)
        VI(LA)=AUXI(NS)
        LA=LA+NUD
  23    CONTINUE
  22    CONTINUE
  24    CONTINUE
  25    CONTINUE
        RETURN
        END
C
C--------------------------------------------------------------
        SUBROUTINE BEM(AK,AKMG,XCR,YCR,KONN,konzon,nonn,
     &    nonzon,MBCU,lbcu,ndtype,zlrar,ZLOAD,ZTEMP,
     &    zloadz,ztempz,ZLOAD1,zloadf,zloadf1,ztempf,ndzon,
     &    NCOUNT,NHPOS,
     &    zlt1,zlt2,zlt3,zld1,zld2,sfdiag,sfdiagzone,SFDIAGnonzone,
     &    zz1,zz2,zz3,N1,N2,N3,N4,N5,n7,n9,icffs,linv,minv)
c            konnot,nonnot,ZGM,ZFM,
C
C       THE BEM SOLUTION BEGINS HERE
C
c          AK,AKMG,XCR,YCR,KONN,konzon,konnot,nonn,
c     &    nonzon,nonnot,MBCU,lbcu,ndtype,ZGM,ZFM,zlrar,ZLOAD,ZTEMP,
C
C       THE BEM SOLUTION BEGINS HERE
C
        IMPLICIT COMPLEX*16 (Z)
        IMPLICIT REAL*8 (A-H,O-Y)
        DIMENSION XCR(N1),YCR(N1),KONN(N3,N4),MBCU(N2)
        dimension konzon(n3,n4)
c              ,konnot(n3,n4)
        dimension nonn(n1,n4)
        dimension nonzon(n1,n4)
c           ,nonnot(n1,n4)
        dimension lbcu(n3),ndtype(n1)
        DIMENSION ZLOAD(N2),ZTEMP(N2)
c            ZGM(N2,N2),ZFM(N2,N2),
        DIMENSION NHPOS(N3)
        DIMENSION ZLOAD1(N5)
        dimension zloadf(n2),ztempf(n2)
c            ,zloadf1(n5)
        dimension sfdiag(n1,2,2),sfdiagnonzone(n1,2,2)
        dimension zlrar(n2+2*ndzon,n2+2*ndzon)
        dimension zloadz(n2+2*ndzon),ztempz(n2+2*ndzon)
c        dimension sfdiagzone(ndzon,2,2)
        dimension sfdiagzone(n1,2,2)







c        dimension zz1(n1),zz2(n1),zz3(n1)
        dimension zz1(n4,n3),zz2(n4,n3),zz3(n4,n3)
c        dimension zlt1(n9,n1),zlt2(n9,n1),zlt3(n9,n1)
        dimension zlt1(n9,n4,n3),zlt2(n9,n4,n3),zlt3(n9,n4,n3)
c        dimension zld1(n9,n1),zld2(n9,n1)
        dimension zld1(n9,n4,n3),zld2(n9,n4,n3)
        dimension linv(n2+2*ndzon),minv(n2+2*ndzon)
c
        COMMON/CONTR/NNODE,NELEM,NP
        common/free/nffs,nfre
c
        ZK=DCMPLX(AK,AKMG)
c       in case that there is free surface, in the subregion case (nffs=2)
c       loop for the 2 steps a,b
c       in the case that there is simulation surface (nffs=3), loop
c       for the 3 steps a,b,c
c         a. ic=1, calculation of the free field of the 
c             i)  non free surface boundary, or
c             ii) simulation boundary
c         b. ic=2, with this free field, calculate the 
c             i)  scattered field on non free surface boundary
c             ii) non simulation boundary (free field)
c
c      in subregion case of calculating the free field(2nd loop),
c      proceed in numerical integration inside the zone
c      ATTENTION: for future versions must take(check) and the nodes of zones


c     for subregion or simulation case, and for the 1st loop
      if(((nffs.eq.2).or.(nffs.eq.3)).and.(icffs.eq.1)) then





c        if((nffs.eq.3).and.
c        if(icffs.eq.1)then




c          CALL nozone(NCOUNT,zk,XCR,YCR,KONN,N1,n2,N3,N4,
c     &             SFDIAGnonzone,lbcu,ndtype,icffs)
          CALL nozone(NCOUNT,zk,XCR,YCR,KONN,N1,n2,N3,N4,
     &             SFDIAGzone,lbcu,ndtype,icffs)
c          CALL nozone(NCOUNT,zk,XCR,YCR,KONN,nonn,N1,n2,N3,N4,
c     &             SFDIAGzone,lbcu,ndtype,icffs)




          CALL zone(NCOUNT,zk,XCR,YCR,KONN,konzon,nonzon,n1,n2,n3,n4,
     &          ndzon,SFDIAGnonzone,lbcu,ndtype,icffs)
c          CALL GENRleft(NCOUNT,KONN,konzon,konnot,MBCU,zlrar,
c     &         N1,N2,n3,N4,ndzon,SFDIAGnonzone,lbcu,ndtype,icffs)
c          CALL GENRleft(NCOUNT,KONN,konzon,konnot,nonzon,nonnot,MBCU,
c     &       zlrar,N1,N2,n3,N4,ndzon,SFDIAGnonzone,lbcu,ndtype,icffs)
c          CALL GENRleft(NCOUNT,KONN,konzon,konnot,nonzon,nonnot,MBCU,
c     &       zlrar,N1,N2,n3,N4,ndzon,SFDIAGzone,lbcu,ndtype,icffs)
c          CALL GENRleft(NCOUNT,KONN,konzon,konnot,nonzon,MBCU,
c     &       zlrar,N1,N2,n3,N4,ndzon,SFDIAGzone,lbcu,ndtype,icffs)
          CALL GENRleft(NCOUNT,KONN,konzon,nonzon,MBCU,
     &       zlrar,N1,N2,n3,N4,ndzon,SFDIAGzone,lbcu,ndtype,icffs)





c         in subregion case of calculating the free field(1st loop),
c         proceed in numerical integration inside the zone

           CALL GENRzone(NCOUNT,KONN,konzon,nonzon,MBCU,zlrar,
     &         N1,N2,n3,N4,ndzon,SFDIAGnonzone,lbcu,ndtype,icffs)
c           CALL GENRright(NCOUNT,KONN,konzon,konnot,nonn,nonzon,
c     &              MBCU,zlrar,ZLOAD,ZTEMP,zloadz,ztempz,N1,N2,
c     &              n3,N4,ndzon,lbcu,ndtype,icffs,linv,minv)
           CALL GENRright(NCOUNT,KONN,konzon,nonn,nonzon,
     &              MBCU,zlrar,ZLOAD,ZTEMP,zloadz,ztempz,N1,N2,
     &              n3,N4,ndzon,lbcu,ndtype,icffs,linv,minv)



c         save the displacements of the 1st step
          rewind 18
          write(18,*) ztemp
c         read the scattered field tractions on the zone nodes, from file 19
c         These tractions are with opposite sign (we suppose that t zone =- t,
c         and so we found t)
c         read the incident tractions
          call clrzrr(zload,n2)
          call clrzrr(zloadf,n2)





          rewind 9
          read(9,*)zload
c         take the incident tractions on the non-zone nodes
          do 29 ii=1,n1
              ln2=nonn(ii,2)
              jref=2*(ln2-1)
            if(ndtype(ii).eq.3) then
              zload(jref+1)=dcmplx(0.0d0,0.0d0)
              zload(jref+2)=dcmplx(0.0d0,0.0d0)
cc             zero displacements
c              ztemp(jref+1)=dcmplx(0.0d0,0.0d0)
c              ztemp(jref+2)=dcmplx(0.0d0,0.0d0)



            endif
cc           zero displacements
c            ztemp(jref+1)=dcmplx(0.0d0,0.0d0)
c            ztemp(jref+2)=dcmplx(0.0d0,0.0d0)


  29      continue

          rewind 19
c          read(19,*)zload
          read(19,*)zloadf
c         add the two fields
          do 19 ii=1,n2
c            zload(ii)=zload(ii)-zloadf(ii)
            zload(ii)=zload(ii)+zloadf(ii)
c            zload(ii)=zload(ii)-zloadf(ii)
c            zload(ii)=-zloadf(ii)
c            zload(ii)=zloadf(ii)



c            zload(ii)=zload(ii)-zloadf(ii)
  19      continue




cc         take the opposite sign
c          do 69 ii=1,n2
c            zload(ii)=-zload(ii)
c  69      continue




c         call the HOOP subroutine, to calculate the radial, tangential and
c         shear stresses on the simulation boundary.
c          CALL HOOP(KONN,nonn,ZLOAD,ZTEMP,ZLOAD1,XCR,YCR,
c     &          n1,n2,n3,n4,n5,zz1,zz2,zz3,zk,icffs,ndtype)
          CALL HOOP(KONN,nonn,ZLOAD,ZTEMP,ZLOAD1,zloadf1,XCR,YCR,
     &          n1,n2,n3,n4,n5,zz1,zz2,zz3,zk,icffs,ndtype)


          go to 41
        endif
c       for subregion or simulation case and 2nd loop
        if(((nffs.eq.2).or.(nffs.eq.3)).and.(icffs.eq.2)) then
c        if((nffs.eq.3).and.(icffs.eq.2)) then


c         in the 2nd step, zload and ztemp
c         contain the free field loads, on the subregion nodes, plus 
c         the incident field on all nodes.
c         read the incident tractions
c          rewind 3
c          read(3,*)zload
          call clrzrr(zload,n2)
          call clrzrr(zloadf,n2)
          call clrzrr(ztemp,n2)
          rewind 9
          read(9,*)zload

cc         for the non-subregion nodes take opposite loads
c          do 29 ii=1,n1
c            if(ndtype(ii).ne.3) then
c              ln2=nonn(ii,2)
c              jref=2*(ln2-1)
c              zload(jref+1)=-zload(jref+1)
c              zload(jref+2)=-zload(jref+2)
c            endif
c  29      continue


c         read the scattered field
          rewind 19
          read(19,*)zloadf
c         read the displacements of the 1st time step
          rewind 18
          read(18,*)ztemp





c         add the two fields (opposite sign for the scattered)
          do 39 ii=1,n1
            ln2=nonn(ii,2)
            jref=2*(ln2-1)
c           accept the displacements of the 1st loop as boundary condition
c           for the free surface in the 2nd loop
            if(ndtype(ii).eq.2) then
              mbcu(jref+1)=1
              mbcu(jref+2)=1
              zload(jref+1)=ztemp(jref+1)
              zload(jref+2)=ztemp(jref+2)
            endif




            if(ndtype(ii).eq.3) then
              zload(jref+1)=zload(jref+1)-zloadf(jref+1)
              zload(jref+2)=zload(jref+2)-zloadf(jref+2)
c              zload(jref+1)=zload(jref+1)+zloadf(jref+1)
c              zload(jref+2)=zload(jref+2)+zloadf(jref+2)



c              zload(jref+1)=-zloadf(jref+1)
c              zload(jref+2)=-zloadf(jref+2)
            endif
c            if(ndtype(ii).eq.2) then
cccc              zload(jref+1)=-zload(jref+1)
cccc              zload(jref+2)=-zload(jref+2)
c              zload(jref+1)=dcmplx(0.0d0,0.0d0)
c              zload(jref+2)=dcmplx(0.0d0,0.0d0)
c            endif




c              zload(ii)=zload(ii)-zloadf(ii)
c              zload(ii)=zload(ii)+zloadf(ii)
c              zload(ii)=-zloadf(ii)


c            zload(ii)=zload(ii)+zloadf(ii)


  39      continue




          CALL NMAIN(NCOUNT,zk,XCR,YCR,KONN,N1,n2,N3,N4,
     &           SFDIAG,lbcu,ndtype,icffs)
c          CALL GENR(NCOUNT,KONN,MBCU,ZGM,ZFM,ZLOAD,ZTEMP,N1,N2,
c     &       n3,N4,SFDIAG,lbcu,ndtype,icffs,linv,minv)
          CALL GENR(NCOUNT,KONN,MBCU,zlrar,ZLOAD,ZTEMP,N1,N2,
     &       n3,N4,SFDIAG,lbcu,ndtype,icffs,linv,minv)
c       read the displacements of the 1st step
        rewind 18
        read(18,*) ztempf
cc       add the displacements of the 2 steps
c        do 49 ii=1,n2
c          ztemp(ii)=ztemp(ii)+ztempf(ii)
cc          ztemp(ii)=ztemp(ii)-ztempf(ii)
cc          ztemp(ii)=ztempf(ii)
c  49    continue










          go to 91
c        endif

      endif
      CALL NMAIN(NCOUNT,ZK,XCR,YCR,KONN,N1,n2,N3,N4,
     &           SFDIAG,lbcu,ndtype,icffs)
c      CALL GENR(NCOUNT,KONN,MBCU,ZGM,ZFM,ZLOAD,ZTEMP,N1,N2,
c     &            n3,N4,SFDIAG,lbcu,ndtype,icffs,linv,minv)
      CALL GENR(NCOUNT,KONN,MBCU,zlrar,ZLOAD,ZTEMP,N1,N2,
     &            n3,N4,SFDIAG,lbcu,ndtype,icffs,linv,minv)
  91  continue
c
      CALL HOOP(KONN,nonn,ZLOAD,ZTEMP,ZLOAD1,zloadf1,XCR,YCR,
     &          n1,n2,n3,n4,n5,zz1,zz2,zz3,zk,icffs,ndtype)



















c     for the 1st loop
c     for subregion or simulation case and 1st loop
      if(((nffs.eq.2).or.(nffs.eq.3)).and.(icffs.eq.1)) go to 41
c      if((nffs.eq.3).and.(icffs.eq.1)) go to 41



cc       read the displacements of the 1st step
c        rewind 18
c        read(18,*) ztempf
cc       add the displacements of the 2 steps
c        do 49 ii=1,n2
c          ztemp(ii)=ztemp(ii)+ztempf(ii)
cc          ztemp(ii)=ztemp(ii)-ztempf(ii)
cc          ztemp(ii)=ztempf(ii)
c  49    continue     







c         tractions
          NP1=NP/3
          DO 1 J=1,NP1
c            k=3*j-2
c           case np=nnode, results of hoop stresses are calculated
c           for all the nodes
c            IF(NP1.EQ.nnode) NJ=J
            IF(NP1.EQ.nelem) NJ=J
c           case np<nnode, results of hoop stresses for some nodes
c           case np<nelem, results of hoop stresses for some elements
            IF(NP1.LT.NELEM) NJ=NHPOS(J)

cc            ZLT(NCOUNT,k)=ZZ1(NJ)
c            ZLT1(NCOUNT,j)=ZZ1(NJ)
cc            k=k+1
cc            ZLT(NCOUNT,k)=ZZ2(NJ)
c            ZLT2(NCOUNT,j)=ZZ2(NJ)
cc            k=k+1
cc            ZLT(NCOUNT,k)=ZZ3(NJ)
c            ZLT3(NCOUNT,j)=ZZ3(NJ)
cc            k1=2*j-1
cc            jref=2*(ln2-1)
cc            zld(ncount,k1)=ztemp(jref)
cc            k1=k1+1
cc            jref=jref+1
cc            zld(ncount,k1)=ztemp(jref)

c           loop for the 3 nodes of the quadratic element
            do 11 l=1,n4
c              ZLT1(NCOUNT,j)=ZZ1(l,NJ)
c              ZLT2(NCOUNT,j)=ZZ2(l,NJ)
c              ZLT3(NCOUNT,j)=ZZ3(l,NJ)
              ZLT1(ncount,l,j)=ZZ1(l,NJ)
              ZLT2(ncount,l,j)=ZZ2(l,NJ)
              ZLT3(ncount,l,j)=ZZ3(l,NJ)
  11        continue
  1       CONTINUE
c
c         displacements
c          do 2 jcl=1,nnode
          do 2 jcl=1,nelem
c            ln2=nonn(jcl,2)
            ln2=konn(jcl,i)
            jref=2*(ln2-1)
c            zld(ncount,jref+1)=ztemp(jref+1)
c            zld(ncount,jref+2)=ztemp(jref+2)
c            zld1(ncount,ln2)=ztemp(jref+1)
c            zld2(ncount,ln2)=ztemp(jref+2)
              zld1(ncount,i,jcl)=ztemp(jref+1)
              zld2(ncount,i,jcl)=ztemp(jref+2)
  2       continue
  41    continue
        RETURN
        END
C----------------------------------------------------------------------
C
        SUBROUTINE NMAIN(NCOUNT,ZK,XCR,YCR,KONN,N1,n2,N3,N4,SFDIAG,
     &                   lbcu,ndtype,icffs)
C
C       START NODAL COLLOCATION SCHEME AND SUBDIVISION INTO SUBELEMENT
C
        IMPLICIT COMPLEX*16 (Z)
        IMPLICIT REAL*8 (A-H,O-Y)
        DIMENSION XCR(N1),YCR(N1),KONN(N3,N4)
        dimension sfdiag(n1,2,2)
c        DIMENSION mbcu(n2)
        dimension lbcu(n3),ndtype(n1)
        DIMENSION XY(3,2),LDIR(3),RLNA(3),SHP(3),DSHP(3),RN(2)
c        COMMON/CONTR/NNODE,NELEM,NTINT,NP
        COMMON/CONTR/NNODE,NELEM,NP
        COMMON/GAUSS/GPN(10,10),WTN(10,10),GPS(10,10),WTS(10,10)
        common/free/nffs,nfre
c
        NOD=3
        TOL1=1.0D-02
        TOL2=1.0D-04
C
        DO 2 ICL=1,NNODE
        XPT=XCR(ICL)
        YPT=YCR(ICL)
        DO 3 JCL=1,NELEM
        DO 4 I=1,NOD
        LN=KONN(JCL,I)
        LDIR(I)=LN
        XY(I,1)=XCR(LN)
        XY(I,2)=YCR(LN)
  4     CONTINUE
C
C       CHECK FOR SINGULAR ELEMENTS
C
        IN=0
        NSUB=1
        DO 6 I=1,NOD
        RLNA(I)=DSQRT((XY(I,1)-XPT)**2+(XY(I,2)-YPT)**2)
        IF (RLNA(I) .LT. TOL1) IN=I
  6     CONTINUE
        IF (IN .EQ. 0) GO TO 10
C
C       SINGULAR CASE SUBDIVISION
C
        NSUB=4
        NPS=6
        CALL INTEGR(NCOUNT,ZK,NPS,ICL,jcl,IN,NSUB,NOD,XPT,
     &         YPT,XY,SFDIAG,n1,n2,n3,lbcu,ndtype,icffs)
        GO TO 3
C
C       NONSINGULAR CASE SUBDIVISION
C
  10    CONTINUE
        RMIN=RLNA(1)
        DO 7 I=2,NOD
        IF (RMIN .GT. RLNA(I)) RMIN=RLNA(I)
  7     CONTINUE
C
        NGP=6
        ELEG=0.0D0
        DO 8 IX=1,NGP
        EK=GPN(NGP,IX)
        WT=WTN(NGP,IX)
        CALL CSHAPE(EK,NOD,XY,SHP,DSHP,RN,DETJ)
        ELEG=ELEG+DETJ*WT
  8     CONTINUE
C
        RATIO=RMIN/ELEG
        IF (RATIO .LT. 1.0D0) NSUB=3
        IF (RATIO .GE. 1.0D0 .AND. RATIO .LE. 2.0D0) NSUB=2
C
        ELEG=ELEG/NSUB
        RAT1=TOL2*ELEG
        RAT2=ELEG/(4.0D0*RMIN)
        RPN=0.5D0*DLOG(RAT1)/DLOG(RAT2)+1.0D0
        NPN=IDINT(RPN)
        IF (NPN .LT. 2) NPN=2
        IF (NPN .GT. 10) NPN=10
        CALL INTEGR(NCOUNT,ZK,NPN,ICL,jcl,IN,NSUB,NOD,XPT,
     &       YPT,XY,SFDIAG,n1,n2,n3,lbcu,ndtype,icffs)
  3     CONTINUE
  2     CONTINUE
        RETURN
        END
C---------------------------------------------------------------
        SUBROUTINE INTEGR(NCOUNT,ZK,NGP,ICL,jcl,IN,NSUB,NOD,
     &      XPT,YPT,XY,SFDIAG,n1,n2,n3,lbcu,ndtype,icffs)
C
C       FORMATION OF 2*2 SUBMATRICES VIA NUMERICAL INTEGRATION
C
        IMPLICIT COMPLEX*16 (Z)
        IMPLICIT REAL*8 (A-H,O-Y)
        DIMENSION SKT(3,2,2),SU(2,2),ST(2,2)
        dimension sfdiag(n1,2,2)
        DIMENSION XY(3,2),YX(3,2),SHP(3),DSHP(3),RN(2)
        DIMENSION ZKU(3,2,2),ZKT(3,2,2),ZU(2,2),ZT(2,2),X(2),Y(2)
c        DIMENSION mbcu(n2)
        dimension lbcu(n3),ndtype(n1)
c
        COMMON/GAUSS/GPN(10,10),WTN(10,10),GPS(10,10),WTS(10,10)
        common/free/nffs,nfre
C
        CALL CLRARR(YX,6)
        CALL CLRZRR(ZKU,12)
        CALL CLRZRR(ZKT,12)
        CALL CLRARR(SKT,12)
        DO 10 IS=1,NSUB
        DO 10 IX=1,NGP
        IF (IN .GT. 0) GO TO 40
C
C       NONSINGULAR CASE
C
        EK=GPN(NGP,IX)
        WT=WTN(NGP,IX)
        GO TO 41
C
C       SINGULAR CASE
C
  40    EK=GPN(NGP,IX)
        WT=WTN(NGP,IX)
  41    CONTINUE
C
C       START AT LINEAR SUBELEMENT LEVEL
C
C                  MODIFICATION       !!!!!!!!!!!
        RINC=2.0D0/DFLOAT(NSUB)
        YX(1,1)=-1.0D0+RINC*DFLOAT(IS-1)
        YX(2,1)=YX(1,1)+RINC
C                 END MODIF    !!!!!!!!!
        CALL CSHAPE(EK,2,YX,SHP,DSHP,RN,DETJ1)
C
C       CONTINUE WITH LOCAL QUADRATIC ELEMENT LEVEL
C
        EKS=0.0D0
        DO 11 I=1,2
        EKS=EKS+YX(I,1)*SHP(I)
  11    CONTINUE
C
        CALL CSHAPE(EKS,NOD,XY,SHP,DSHP,RN,DETJ2)
        DETJ=DETJ1*DETJ2
C
C       NOW AT GLOBAL QUADRATIC ELEMENT LEVEL
C
        XG=0.0D0
        YG=0.0D0
        DO 12 I=1,NOD
        XG=XG+XY(I,1)*SHP(I)
        YG=YG+XY(I,2)*SHP(I)
  12    CONTINUE
        X(1)=XG
        X(2)=YG
        Y(1)=XPT
        Y(2)=YPT
        CALL KERNEL(ZU,ZT,X,Y,RN,ZK)
        CALL ZERNEL(XPT,YPT,XG,YG,RN,SU,ST)

cccccccccccccccccccccccccccccccccccccccccccccccccccc
cc             ATTENTION    !!!!!!!!!!!!!!!!
cc       for free surface nodes in the 2nd loop the displacement
cc       kernel is 0, because it is multiplyed by the zero tractions(with
cc       opposite sign) on the free surface.
cc       free surface node or element
cc        if((ndtype(icl).eq.2).or.(lbcu(jcl).eq.1)) then


cc       case of free surface
c        if((lbcu(jcl).eq.2)) then
c          do 52 ii=1,2
c            do 52 jj=1,2
c              zt(ii,jj)=dcmplx(0.0d0,0.0d0)
c              st(ii,jj)=0.0d0
c  52      continue
c        endif

cc        if (lbcu(jcl).eq.2) then
c         if((ndtype(icl).eq.2).or.(lbcu(jcl).eq.2)) then
c          do 50 ii=1,2
c            do 50 jj=1,2 
c              zu(ii,jj)=dcmplx(0.0d0,0.0d0)
c              su(ii,jj)=0.0d0
c  50      continue
c        endif



c **************************
c case of free surface
C this gives better results
         if((lbcu(jcl).eq.2)) then
          do 52 ii=1,2
            do 52 jj=1,2
              zt(ii,jj)=dcmplx(0.0d0,0.0d0)
              st(ii,jj)=0.0d0
  52      continue
        endif
c **************************



        IF (IN .NE. 0) GO TO 18
C
C       NUMERICAL INTEGRATION - NONSINGULAR CASE
C
        DO 16 K1=1,NOD
        AMUL=SHP(K1)*DETJ*WT
        DO 17 K2=1,2
        DO 17 K3=1,2
        ZKU(K1,K2,K3)=ZKU(K1,K2,K3)+ZU(K2,K3)*AMUL
        ZKT(K1,K2,K3)=ZKT(K1,K2,K3)+ZT(K2,K3)*AMUL
        IF (NCOUNT .EQ. 1)
     &    SKT(K1,K2,K3)=SKT(K1,K2,K3)+ST(K2,K3)*AMUL
  17    CONTINUE
  16    CONTINUE
        GO TO 10
  18    CONTINUE
C
C       NUMERICAL INTEGRATION - SINGULAR CASE
C
        DO 19 K1=1,NOD
        AMUL=SHP(K1)*DETJ*WT
        DO 20 K2=1,2
        DO 20 K3=1,2
        ZKU(K1,K2,K3)=ZKU(K1,K2,K3)+ZU(K2,K3)*AMUL
        IF (K1 .EQ. IN) GO TO 21
        SKT(K1,K2,K3)=SKT(K1,K2,K3)+ST(K2,K3)*AMUL
        ZKT(K1,K2,K3)=ZKT(K1,K2,K3)+ZT(K2,K3)*AMUL
        GO TO 20
  21    ZKT(K1,K2,K3)=ZKT(K1,K2,K3)+
     &          (ZT(K2,K3)-ST(K2,K3))*AMUL

ccccccccccccccccccccccccccccccccccccccccccccccccccc
cc #################################################
cc  modification in accordance to the same subroutine of
cc  the program BEMTD (time domain)
c        IF (K1 .EQ. IN.and.ncount.eq.1) GO TO 21
c        if(ncount.ne.1) go to 99
c        SKT(K1,K2,K3)=SKT(K1,K2,K3)+ST(K2,K3)*AMUL
cc        ZKT(K1,K2,K3)=ZKT(K1,K2,K3)+ZT(K2,K3)*AMUL
c   99   ZKT(K1,K2,K3)=ZKT(K1,K2,K3)+ZT(K2,K3)*AMUL
c        GO TO 20
c  21    ZKT(K1,K2,K3)=ZKT(K1,K2,K3)+
c     &          (ZT(K2,K3)-ST(K2,K3))*AMUL
  20    CONTINUE
  19    CONTINUE
  10    CONTINUE
C
C       SUBMATRICES ARE STORED ON TAPES FOR LATER ASSEMBLY
C
        DO 28 K1=1,NOD
        WRITE(7,*) ((ZKU(K1,K2,K3),K2=1,2),K3=1,2)
  28    CONTINUE
        DO 29 K1=1,NOD
        WRITE(8,*) ((ZKT(K1,K2,K3),K2=1,2),K3=1,2)
  29    CONTINUE
C
C       DIAGONAL TERMS
C
        IF (NCOUNT .GT. 1) RETURN
        DO 30 K1=1,NOD
        IF (K1 .EQ. IN) GO TO 30
        DO 31 K2=1,2
        DO 31 K3=1,2
        SFDIAG(ICL,K2,K3)=SFDIAG(ICL,K2,K3)
     &                    -SKT(K1,K2,K3)
  31    CONTINUE
  30    CONTINUE
        RETURN
        END

C----------------------------------------------------------------------
C
        SUBROUTINE nozone(NCOUNT,zk,XCR,YCR,KONN,N1,n2,N3,N4,SFDIAG,
     &                   lbcu,ndtype,icffs)
c        SUBROUTINE nozone(NCOUNT,zk,XCR,YCR,KONN,nonn,N1,n2,N3,N4,
c     &                   SFDIAG,lbcu,ndtype,icffs)



C
C       START NODAL COLLOCATION SCHEME AND SUBDIVISION INTO SUBELEMENT
C
        IMPLICIT COMPLEX*16 (Z)
        IMPLICIT real*8 (a-h,o-y)
        DIMENSION XCR(N1),YCR(N1),KONN(N3,N4)
        dimension sfdiag(n1,2,2)
c        DIMENSION mbcu(n2)
        dimension lbcu(n3),ndtype(n1)
c        dimension nonn(n1,n4)


        DIMENSION XY(3,2),LDIR(3),RLNA(3),SHP(3),DSHP(3),RN(2)
c        COMMON/CONTR/NNODE,NELEM,NTINT,NP
        COMMON/CONTR/NNODE,NELEM,NP
        COMMON/GAUSS/GPN(10,10),WTN(10,10),GPS(10,10),WTS(10,10)
ccccccc
ccccccc
ccccccc
        common/free/nffs,nfre
ccccccc
ccccccc
ccccccc
        NOD=3
        TOL1=1.0D-02
        TOL2=1.0D-04
C
        DO 2 ICL=1,NNODE
        XPT=XCR(ICL)
        YPT=YCR(ICL)
        DO 3 JCL=1,NELEM
        DO 4 I=1,NOD
cc       change the coordinates of the connected nodes, in order
cc       to obtain inward normal unit vector (zone)
c        if(ndtype(jcl).eq.3)  then
c          LN=nONN(JCL,nod+1-I)
c          go to 83  
c        endif  
        LN=KONN(JCL,I)
c  83    continue
        LDIR(I)=LN
        XY(I,1)=XCR(LN)
        XY(I,2)=YCR(LN)
  4     CONTINUE
C
C       CHECK FOR SINGULAR ELEMENTS
C
        IN=0
        NSUB=1
        DO 6 I=1,NOD
        RLNA(I)=DSQRT((XY(I,1)-XPT)**2+(XY(I,2)-YPT)**2)
        IF (RLNA(I) .LT. TOL1) IN=I
  6     CONTINUE
        IF (IN .EQ. 0) GO TO 10
C
C       SINGULAR CASE SUBDIVISION
C
        NSUB=4
        NPS=6
        CALL intnozon(NCOUNT,zk,NPS,ICL,jcl,IN,NSUB,NOD,XPT,
     &              YPT,XY,SFDIAG,n1,n2,n3,ndtype,lbcu,icffs)
        GO TO 3
C
C       NONSINGULAR CASE SUBDIVISION
C
  10    CONTINUE
        RMIN=RLNA(1)
        DO 7 I=2,NOD
        IF (RMIN .GT. RLNA(I)) RMIN=RLNA(I)
  7     CONTINUE
C
        NGP=6
        ELEG=0.0D0
        DO 8 IX=1,NGP
        EK=GPN(NGP,IX)
        WT=WTN(NGP,IX)
        CALL CSHAPE(EK,NOD,XY,SHP,DSHP,RN,DETJ)
        ELEG=ELEG+DETJ*WT
  8     CONTINUE
C
        RATIO=RMIN/ELEG
        IF (RATIO .LT. 1.0D0) NSUB=3
        IF (RATIO .GE. 1.0D0 .AND. RATIO .LE. 2.0D0) NSUB=2
C
        ELEG=ELEG/NSUB
        RAT1=TOL2*ELEG
        RAT2=ELEG/(4.0D0*RMIN)
        RPN=0.5D0*DLOG(RAT1)/DLOG(RAT2)+1.0D0
        NPN=IDINT(RPN)
        IF (NPN .LT. 2) NPN=2
        IF (NPN .GT. 10) NPN=10
        CALL intnozon(NCOUNT,zk,NPN,ICL,jcl,IN,NSUB,NOD,XPT,
     &            YPT,XY,SFDIAG,n1,n2,n3,ndtype,lbcu,icffs)
  3     CONTINUE
  2     CONTINUE
        RETURN
        END
c
c
C---------------------------------------------------------------
        SUBROUTINE intnozon(NCOUNT,zk,NGP,ICL,jcl,IN,NSUB,NOD,
     &          XPT,YPT,XY,SFDIAG,n1,n2,n3,ndtype,lbcu,icffs)
C
C       FORMATION OF 2*2 SUBMATRICES VIA NUMERICAL INTEGRATION
C
        IMPLICIT COMPLEX*16 (Z)
        IMPLICIT real*8 (a-h,o-y)
        DIMENSION SKT(3,2,2),SU(2,2),ST(2,2)
        dimension sfdiag(n1,2,2)
        DIMENSION XY(3,2),YX(3,2),SHP(3),DSHP(3),RN(2)
        DIMENSION ZKU(3,2,2),ZKT(3,2,2),ZU(2,2),ZT(2,2),X(2),Y(2)
c        DIMENSION mbcu(n2)
        dimension lbcu(n3),ndtype(n1)

c
        COMMON/GAUSS/GPN(10,10),WTN(10,10),GPS(10,10),WTS(10,10)
        common/free/nffs,nfre
C
        CALL CLRARR(YX,6)
        CALL CLRZRR(ZKU,12)
        CALL CLRZRR(ZKT,12)
        CALL CLRARR(SKT,12)

        DO 10 IS=1,NSUB
        DO 10 IX=1,NGP
        IF (IN .GT. 0) GO TO 40
C
C       NONSINGULAR CASE
C
        EK=GPN(NGP,IX)
        WT=WTN(NGP,IX)
        GO TO 41
C
C       SINGULAR CASE
C
  40    EK=GPN(NGP,IX)
        WT=WTN(NGP,IX)
  41    CONTINUE
C
C       START AT LINEAR SUBELEMENT LEVEL
C
        RINC=2.0D0/DFLOAT(NSUB)
        YX(1,1)=-1.0D0+RINC*DFLOAT(IS-1)
        YX(2,1)=YX(1,1)+RINC
        CALL CSHAPE(EK,2,YX,SHP,DSHP,RN,DETJ1)
C       CONTINUE WITH LOCAL QUADRATIC ELEMENT LEVEL
C
        EKS=0.0D0
        DO 11 I=1,2
        EKS=EKS+YX(I,1)*SHP(I)
  11    CONTINUE
C
        CALL CSHAPE(EKS,NOD,XY,SHP,DSHP,RN,DETJ2)
        DETJ=DETJ1*DETJ2
C
C       NOW AT GLOBAL QUADRATIC ELEMENT LEVEL
C
        XG=0.0D0
        YG=0.0D0
        DO 12 I=1,NOD
        XG=XG+XY(I,1)*SHP(I)
        YG=YG+XY(I,2)*SHP(I)
  12    CONTINUE
        X(1)=XG
        X(2)=YG
        Y(1)=XPT
        Y(2)=YPT
        CALL KERNEL(ZU,ZT,X,Y,RN,ZK)

cccccccc
cccccccc
cccccccc
c
        CALL ZERNEL(XPT,YPT,XG,YG,RN,SU,ST)

ccccccccccccccccccccccccccccccccccccccccccccccccccc
c             ATTENTION    !!!!!!!!!!!!!!!!
c       for free surface nodes in this 1st loop the displacement
c       kernel is not 0, because in the right hand side of the
c       equation is multiplyed by the non zero incident tractions(with
c      opposite sign) on the free surface.
cc      
ccc       for free surface node, and in the 2nd loop the displacement
cc       kernel is  0, because in the right hand side of the
cc       equation is multiplyed by the zero incident tractions
c        if((icffs.eq.2).and. ((ndtype(icl).eq.2).or.
c     &              (lbcu(jcl).eq.2))) then
c          do 50 ii=1,2
c            do 50 jj=1,2
c              zu(ii,jj)=dcmplx(0.0d0,0.0d0)
c              su(ii,jj)=0.0d0
c  50      continue
c        endif













        IF (IN .NE. 0) GO TO 18
C
C       NUMERICAL INTEGRATION - NONSINGULAR CASE
C
        DO 16 K1=1,NOD
        AMUL=SHP(K1)*DETJ*WT
        DO 17 K2=1,2
        DO 17 K3=1,2
        ZKU(K1,K2,K3)=ZKU(K1,K2,K3)+ZU(K2,K3)*AMUL
        ZKT(K1,K2,K3)=ZKT(K1,K2,K3)+ZT(K2,K3)*AMUL
        IF (NCOUNT .EQ. 1)
     &    SKT(K1,K2,K3)=SKT(K1,K2,K3)+ST(K2,K3)*AMUL
  17    CONTINUE
  16    CONTINUE
        GO TO 10
  18    CONTINUE
C
C       NUMERICAL INTEGRATION - SINGULAR CASE
C
        DO 19 K1=1,NOD
        AMUL=SHP(K1)*DETJ*WT
        DO 20 K2=1,2
        DO 20 K3=1,2
        ZKU(K1,K2,K3)=ZKU(K1,K2,K3)+ZU(K2,K3)*AMUL
        IF (K1 .EQ. IN) GO TO 21
        SKT(K1,K2,K3)=SKT(K1,K2,K3)+ST(K2,K3)*AMUL
        ZKT(K1,K2,K3)=ZKT(K1,K2,K3)+ZT(K2,K3)*AMUL
        GO TO 20
  21    ZKT(K1,K2,K3)=ZKT(K1,K2,K3)+
     &          (ZT(K2,K3)-ST(K2,K3))*AMUL

  20    CONTINUE
  19    CONTINUE
  10    CONTINUE
C
C       SUBMATRICES ARE STORED ON TAPES FOR LATER ASSEMBLY
C
        DO 28 K1=1,NOD
cc         in subregion case, put the opposite sign for the subregion
cc         boundary only (unit normal opposite)
c          if(lbcu(jcl).eq.3) then
c            do 64 k2=1,2
c              do 64 k3=1,2
c                sku(k1,k2,k3)=-sku(k1,k2,k3)
c  64        continue
c          endif
          WRITE(7,*) ((ZKU(K1,K2,K3),K2=1,2),K3=1,2)
  28    CONTINUE
        DO 29 K1=1,NOD
cc         in subregion case, put the opposite sign for the subregion
cc         boundary only (unit normal opposite)
c          if(lbcu(jcl).eq.3) then 
c            do 66 k2=1,2
c              do 66 k3=1,2
c                skt(k1,k2,k3)=-skt(k1,k2,k3)
c  66        continue
c          endif
        WRITE(8,*) ((ZKT(K1,K2,K3),K2=1,2),K3=1,2)
  29    CONTINUE
C
C       DIAGONAL TERMS
C
        IF (NCOUNT .GT. 1) RETURN
        DO 30 K1=1,NOD
        IF (K1 .EQ. IN) GO TO 30
        DO 31 K2=1,2
        DO 31 K3=1,2
        SFDIAG(ICL,K2,K3)=SFDIAG(ICL,K2,K3)
     &                    -SKT(K1,K2,K3)
  31    CONTINUE
  30    CONTINUE
        RETURN
        END
c
C------------------------------------------------
      SUBROUTINE zone(NCOUNT,zk,XCR,YCR,KONN,konzon,nonzon,
     &     n1,n2,n3,n4,ndzon,SFDIAG,lbcu,ndtype,icffs)
C
C     START NODAL COLLOCATION SCHEME AND SUBDIVISION INTO SUBELEMENT
C 
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT real*8 (a-h,o-y)
      DIMENSION XCR(N1),YCR(N1),KONN(N3,N4)
c      dimension sfdiag(n1,2,2)
c      dimension sfdiagzone(ndzon,2,2)
c      dimension sfdiag(ndzon,2,2)
      dimension sfdiag(n1,2,2)














c      DIMENSION mbcu(n2)
      dimension lbcu(n3),ndtype(n1)
      DIMENSION XY(3,2),LDIR(3),RLNA(3),SHP(3),DSHP(3),RN(2)
c      COMMON/CONTR/NNODE,NELEM,NTINT,NP
      COMMON/CONTR/NNODE,NELEM,NP
      COMMON/GAUSS/GPN(10,10),WTN(10,10),GPS(10,10),WTS(10,10)
      dimension konzon(n3,n4)
      dimension nonzon(n1,n4)
      common/free/nffs,nfre
c
      NOD=3
      TOL1=1.0D-02
      TOL2=1.0D-04
C
      DO 2 ICL=1,NNODE
c       take only the nodes with ndtype=1(inclusion case),
c                             or ndtype=3(simulation case)
c        if((ndtype(icl).ne.1).or.(ndtype(icl).ne.3)) go to 2
        if((ndtype(icl).eq.1).or.(ndtype(icl).eq.3)) then
          XPT=XCR(ICL)
          YPT=YCR(ICL)
          DO 3 JCL=1,NELEM
c           take only the elements with lbcu=1(inclusion case),
c                                or lbcu=3(simulation case)
c            if((lbcu(jcl).ne.1).or.(lbcu(jcl).ne.3)) go to 3
            if((lbcu(jcl).eq.1).or.(lbcu(jcl).eq.3)) then
              DO 4 I=1,NOD
c                LN=KONN(JCL,I)
c             change the coordinates of the connected nodes, in order
c             to obtain inward normal unit vector (zone)
                LN=KONN(JCL,nod+1-I)



                LDIR(I)=LN
                XY(I,1)=XCR(LN)
                XY(I,2)=YCR(LN)
  4           CONTINUE



C
C             CHECK FOR SINGULAR ELEMENTS
C
              IN=0
              NSUB=1
              DO 6 I=1,NOD
                RLNA(I)=DSQRT((XY(I,1)-XPT)**2+(XY(I,2)-YPT)**2)
                IF (RLNA(I) .LT. TOL1) IN=I
  6           CONTINUE
              IF (IN .EQ. 0) GO TO 10
C
C             SINGULAR CASE SUBDIVISION
C
              NSUB=4
              NPS=6
              CALL intzon(NCOUNT,zk,NPS,ICL,jcl,IN,NSUB,NOD,XPT,
     &          YPT,XY,SFDIAG,n1,n2,n3,n4,ndzon,ndtype,lbcu,
     &          icffs,konzon,nonzon)
              GO TO 3
C
C             NONSINGULAR CASE SUBDIVISION
C
  10          CONTINUE
              RMIN=RLNA(1)
              DO 7 I=2,NOD
                IF (RMIN .GT. RLNA(I)) RMIN=RLNA(I)
  7           CONTINUE
C
              NGP=6
              ELEG=0.0D0
              DO 8 IX=1,NGP
                EK=GPN(NGP,IX)
                WT=WTN(NGP,IX)
                CALL CSHAPE(EK,NOD,XY,SHP,DSHP,RN,DETJ)
                ELEG=ELEG+DETJ*WT
  8           CONTINUE
C
              RATIO=RMIN/ELEG
              IF (RATIO .LT. 1.0D0) NSUB=3
              IF (RATIO .GE. 1.0D0 .AND. RATIO .LE. 2.0D0) NSUB=2
C
              ELEG=ELEG/NSUB
              RAT1=TOL2*ELEG
              RAT2=ELEG/(4.0D0*RMIN)
              RPN=0.5D0*DLOG(RAT1)/DLOG(RAT2)+1.0D0
              NPN=IDINT(RPN)
              IF (NPN .LT. 2) NPN=2
              IF (NPN .GT. 10) NPN=10
              CALL intzon(NCOUNT,zk,NPN,ICL,jcl,IN,NSUB,NOD,XPT,
     &          YPT,XY,SFDIAG,n1,n2,n3,n4,ndzon,ndtype,lbcu,
     &          icffs,konzon,nonzon)
            endif
  3       CONTINUE
        endif
  2   CONTINUE
      RETURN
      END
c
C---------------------------------------------------------------
        SUBROUTINE intzon(NCOUNT,zk,NGP,ICL,jcl,IN,NSUB,NOD,
     &    XPT,YPT,XY,SFDIAG,n1,n2,n3,n4,ndzon,ndtype,lbcu,
     &    icffs,konzon,nonzon)
C
C       FORMATION OF 2*2 SUBMATRICES VIA NUMERICAL INTEGRATION
C
        IMPLICIT COMPLEX*16 (Z)
        IMPLICIT real*8 (a-h,o-y)
        DIMENSION SKT(3,2,2),SU(2,2),ST(2,2)
c        dimension sfdiagzone(ndzon,2,2)
c        dimension sfdiag(ndzon,2,2)
        dimension sfdiag(n1,2,2)









        dimension konzon(n3,n4)
        dimension nonzon(n1,n4)
        DIMENSION XY(3,2),YX(3,2),SHP(3),DSHP(3),RN(2)
        DIMENSION ZKU(3,2,2),ZKT(3,2,2),ZU(2,2),ZT(2,2),X(2),Y(2)
c        DIMENSION mbcu(n2)
        dimension lbcu(n3),ndtype(n1)
c
        COMMON/GAUSS/GPN(10,10),WTN(10,10),GPS(10,10),WTS(10,10)
        common/free/nffs,nfre
C
        CALL CLRARR(YX,6)
        CALL CLRZRR(ZKU,12)
        CALL CLRZRR(ZKT,12)
        CALL CLRARR(SKT,12)

        DO 10 IS=1,NSUB
        DO 10 IX=1,NGP
        IF (IN .GT. 0) GO TO 40
C
C       NONSINGULAR CASE
C
        EK=GPN(NGP,IX)
        WT=WTN(NGP,IX)
        GO TO 41
C
C       SINGULAR CASE
C
  40    EK=GPN(NGP,IX)
        WT=WTN(NGP,IX)
  41    CONTINUE
C
C       START AT LINEAR SUBELEMENT LEVEL
C
C                  MODIFICATION       !!!!!!!!!!!
        RINC=2.0D0/DFLOAT(NSUB)
        YX(1,1)=-1.0D0+RINC*DFLOAT(IS-1)
        YX(2,1)=YX(1,1)+RINC
C                 END MODIF    !!!!!!!!!
        CALL CSHAPE(EK,2,YX,SHP,DSHP,RN,DETJ1)
C
C       CONTINUE WITH LOCAL QUADRATIC ELEMENT LEVEL
C
        EKS=0.0D0
        DO 11 I=1,2
        EKS=EKS+YX(I,1)*SHP(I)
  11    CONTINUE
C
        CALL CSHAPE(EKS,NOD,XY,SHP,DSHP,RN,DETJ2)
        DETJ=DETJ1*DETJ2
C
C       NOW AT GLOBAL QUADRATIC ELEMENT LEVEL
C
        XG=0.0D0
        YG=0.0D0
        DO 12 I=1,NOD
        XG=XG+XY(I,1)*SHP(I)
        YG=YG+XY(I,2)*SHP(I)
  12    CONTINUE
        X(1)=XG
        X(2)=YG
        Y(1)=XPT
        Y(2)=YPT
        CALL KERNEL(ZU,ZT,X,Y,RN,zk)
        CALL ZERNEL(XPT,YPT,XG,YG,RN,SU,ST)

ccccccccccccccccccccccccccccccccccccccccccccccccccc
c         ATTENTION    !!!!!!!!!!!
c       the displacement kernel isn't 0 because in this 1st step
c       the incident traction on the free surface isn't 0
c       for free surface node
c        if((ndtype(icl).eq.2).or.(lbcu(jcl).eq.2)) then
c          do 50 ii=1,2
c            do 50 jj=1,2
c              zu(ii,jj)=dcmplx(0.0d0,0.0d0)
c              su(ii,jj)=0.0d0
c  50      continue
c        endif
        IF (IN .NE. 0) GO TO 18
C
C       NUMERICAL INTEGRATION - NONSINGULAR CASE
C
        DO 16 K1=1,NOD
        AMUL=SHP(K1)*DETJ*WT
        DO 17 K2=1,2
        DO 17 K3=1,2
        ZKU(K1,K2,K3)=ZKU(K1,K2,K3)+ZU(K2,K3)*AMUL
        ZKT(K1,K2,K3)=ZKT(K1,K2,K3)+ZT(K2,K3)*AMUL
        IF (NCOUNT .EQ. 1)
     &    SKT(K1,K2,K3)=SKT(K1,K2,K3)+ST(K2,K3)*AMUL
  17    CONTINUE
  16    CONTINUE
        GO TO 10
  18    CONTINUE
C
C       NUMERICAL INTEGRATION - SINGULAR CASE
C
        DO 19 K1=1,NOD
        AMUL=SHP(K1)*DETJ*WT
        DO 20 K2=1,2
        DO 20 K3=1,2
        ZKU(K1,K2,K3)=ZKU(K1,K2,K3)+ZU(K2,K3)*AMUL
        IF (K1 .EQ. IN) GO TO 21
        SKT(K1,K2,K3)=SKT(K1,K2,K3)+ST(K2,K3)*AMUL
        ZKT(K1,K2,K3)=ZKT(K1,K2,K3)+ZT(K2,K3)*AMUL
        GO TO 20
  21    ZKT(K1,K2,K3)=ZKT(K1,K2,K3)+
     &          (ZT(K2,K3)-ST(K2,K3))*AMUL

  20    CONTINUE
  19    CONTINUE
  10    CONTINUE
C
C       SUBMATRICES ARE STORED ON TAPES FOR LATER ASSEMBLY
C
        DO 28 K1=1,NOD
cc         put the opposite sign for the subregion
cc         boundary (unit normal opposite)
cc          if(lbcu(jcl).eq.3) then
c          do 64 k2=1,2
c            do 64 k3=1,2
c              zku(k1,k2,k3)=-zku(k1,k2,k3)
c  64      continue
        WRITE(14,*) ((zKU(K1,K2,K3),K2=1,2),K3=1,2)
  28    CONTINUE
        DO 29 K1=1,NOD
c          do 66 k2=1,2
c            do 66 k3=1,2
c              zkt(k1,k2,k3)=-zkt(k1,k2,k3)
c              skt(k1,k2,k3)=-skt(k1,k2,k3)
c  66        continue
        WRITE(15,*) ((zKT(K1,K2,K3),K2=1,2),K3=1,2)
  29    CONTINUE


C
C       DIAGONAL TERMS
C
        IF (NCOUNT .GT. 1) RETURN
c        ln2=nonzon(icl,2)


















        DO 30 K1=1,NOD
        IF (K1 .EQ. IN) GO TO 30
        DO 31 K2=1,2
        DO 31 K3=1,2
        SFDIAG(ICL,K2,K3)=SFDIAG(ICL,K2,K3)
     &                    -SKT(K1,K2,K3)
c        SFDIAG(ln2,K2,K3)=SFDIAG(ln2,K2,K3)
c     &                    -SKT(K1,K2,K3)
















  31    CONTINUE
  30    CONTINUE
        RETURN
        END
c
C------------------------------------------------
        SUBROUTINE CSHAPE(EK,NOD,XY,SHP,DSHP,RN,DETJ)
C
C       SHAPE FUNCTIONS, OUTWARD NORMAL, AND DET. OF JACOBIAN FOR LINE
C       ELEMENTS LINEAR AND QUADRATIC TYPES
C
        IMPLICIT REAL*8 (A-H,O-Y)
        DIMENSION XY(3,2),SHP(3),DSHP(3),RN(2),AJAC(2)
        tol=1.0e-11
c        pi=3.14159265359d0
C
        IF (NOD .EQ. 3) GO TO 1
C
C       LINEAR
C
        SHP(1)=0.5d0*(1.0d0-EK)
        SHP(2)=0.5d0*(1.0d0+EK)
C
        DSHP(1)=-0.5d0
        DSHP(2)=0.5d0
        GO TO 2
C
C       QUADRATIC
C
  1     CONTINUE
        SHP(1)=-0.5d0*EK+0.5d0*EK**2
        SHP(2)=1.0d0-EK**2
        SHP(3)=0.5d0*EK+0.5d0*EK**2
C
        DSHP(1)=-0.5d0+EK
        DSHP(2)=-2.0d0*EK
        DSHP(3)=0.5d0+EK
  2     CONTINUE
C
C       JACOBIAN MATRIX, DETERMINANT, NORMAL VECTOR
C
        DO 10 I=1,2
        AJAC(I)=0.0D0
        DO 10 J=1,NOD
        AJAC(I)=AJAC(I)+DSHP(J)*XY(J,I)
  10    CONTINUE
C
        DETJ=DSQRT(AJAC(1)**2+AJAC(2)**2)

c       calculate the rn(1)= cos of normal to the element unit outward vector
c                            with the global x-axis
c                            (direction cosine of unit outward normal vector)
c                     rn(2)= sin of normal to the element unit outward vector
c                            with the global x-axis
        aj1=ajac(1)
        aj2=ajac(2)
c       detj always possitive or 0
        if(detj.lt.tol)then
c          if (detj.lt.0.0d0) detj=-tol
c          if (detj.gt.0.0d0) detj=tol
          if (detj.eq.0.0d0) then
            if (aj2.gt.0.0d0) rn(1)=0.0d0
            if (aj2.lt.0.0d0) rn(1)=0.0d0
            if (aj1.lt.0.0d0) rn(2)=1.0d0
            if (aj2.gt.0.0d0) rn(2)=-1.0d0
            go to 34
          endif
        endif
c        RN(1)=AJAC(2)/DETJ
c        RN(2)=-AJAC(1)/DETJ
        RN(1)=aj2/detj
        RN(2)=-aj1/detj
  34    continue
        RETURN
        END
cC------------------------------------------------
c        SUBROUTINE CSHAPE1(EK,NOD,XY,SHP,DSHP,RN,DETJ,ajac)
cC
cC       SHAPE FUNCTIONS, OUTWARD NORMAL, AND DET. OF JACOBIAN FOR LINE
cC       ELEMENTS LINEAR AND QUADRATIC TYPES
cC
c        IMPLICIT REAL*8 (A-H,O-Y)
c        DIMENSION XY(3,2),SHP(3),DSHP(3),RN(2),AJAC(2)
c        tol=1.0e-11
cc        pi=3.14159265359d0
cC
c        IF (NOD .EQ. 3) GO TO 1
cC
cC       LINEAR
cC
c        SHP(1)=0.5d0*(1.0d0-EK)
c        SHP(2)=0.5d0*(1.0d0+EK)
cC
c        DSHP(1)=-0.5d0
c        DSHP(2)=0.5d0
c        GO TO 2
cC
cC       QUADRATIC
cC
c  1     CONTINUE
c        SHP(1)=-0.5d0*EK+0.5d0*EK**2
c        SHP(2)=1.0d0-EK**2
c        SHP(3)=0.5d0*EK+0.5d0*EK**2
cC
c        DSHP(1)=-0.5d0+EK
c        DSHP(2)=-2.0d0*EK
c        DSHP(3)=0.5d0+EK
c  2     CONTINUE
cC
cC       JACOBIAN MATRIX, DETERMINANT, NORMAL VECTOR
cC
c        DO 10 I=1,2
c        AJAC(I)=0.0D0
c        DO 10 J=1,NOD
c        AJAC(I)=AJAC(I)+DSHP(J)*XY(J,I)
c  10    CONTINUE
cC
c        DETJ=DSQRT(AJAC(1)**2+AJAC(2)**2)
c
cc       calculate the rn(1)= cos of normal to the element unit outward vector
cc                            with the global x-axis
cc                            (direction cosine of unit outward normal vector)
cc                     rn(2)= sin of normal to the element unit outward vector
cc                            with the global x-axis
c        aj1=ajac(1)
c        aj2=ajac(2)
cc       detj always possitive or 0
c        if(detj.lt.tol)then
cc          if (detj.lt.0.0d0) detj=-tol
cc          if (detj.gt.0.0d0) detj=tol
c          if (detj.eq.0.0d0) then
c            if (aj2.gt.0.0d0) rn(1)=0.0d0
c            if (aj2.lt.0.0d0) rn(1)=0.0d0
c            if (aj1.lt.0.0d0) rn(2)=1.0d0
c            if (aj2.gt.0.0d0) rn(2)=-1.0d0
c            go to 34
c          endif
c        endif
cc        RN(1)=AJAC(2)/DETJ
cc        RN(2)=-AJAC(1)/DETJ
c        RN(1)=aj2/detj
c        RN(2)=-aj1/detj
c  34    continue
c        RETURN
c        END
C------------------------------------------------
        SUBROUTINE KERNEL(ZU,ZT,X,Y,RN,Z)
C
C       EVALUAT. OF 2-D FUNDAMENTAL SOLUTIONS IN LAPLACE TRANSF. DOMAIN
C
C       ZU IS DISPLACEMENT FUNDAMENTAL SOLUTION
C       ZT IS TRACTION FUNDAMENTAL SOLUTION
C       X,Y ARE FIELD AND SOURCE POINTS, RN IS NORMAL VECTOR
C       Z IS LAPLACE TRANSFORM PARAMETER
C
        IMPLICIT REAL*8 (A-H,O-Y)
        IMPLICIT COMPLEX*16 (Z)
c                        ##########################
c                        !!!!!!!!!!!!!!!!!!!!!!!!!!
c                 modification for viscoelesticity
        COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM,fvisc,nmode
        DIMENSION ZU(2,2),ZT(2,2),X(2),Y(2),RN(2),RC(2)
        DIMENSION D(2,2)
c        DATA D/1.0D0,0.0D0,0.0D0,1.0D0/
c       Kronecker delta
        d(1,1)=1.0d0
        d(2,1)=0.0d0
        d(1,2)=0.0d0
        d(2,2)=1.0d0
C
c                        ##########################
c                        !!!!!!!!!!!!!!!!!!!!!!!!!!
cc                 modification for viscoelesticity
        tol=1.0d-9
cc       Kelvin solid
c        if(dabs(fvisc).lt.tol) then
c          zmu=rmu
c          go to 5
c        endif
cc        if(nmode.eq.1)zmu=rmu+fvisc*rmu*z
c        if(nmode.eq.1)zmu=rmu+fvisc*z
cc       Maxwell fluid
cc        if(nmode.eq.2)zmu=rmu*z/(z+1.0d0/fvisc)
c        if(nmode.eq.2)zmu=rmu*z/(z+rmu/fvisc)
cc
c   5    continue

c       find the complex values of Lame's constants
        call viscoel(z,zmu,zla,znu,zzc1,zzc2) 
        PI=3.14159265359D0
        R=0.0D0
        DO 10 I=1,2
        R=R+(X(I)-Y(I))**2
  10    CONTINUE
        R=DSQRT(R)
C
        DO 20 I=1,2
        RC(I)=(X(I)-Y(I))/R
  20    CONTINUE
C
        DRN=0.0D0
        DO 21 I=1,2
        DRN=DRN+RC(I)*RN(I)
  21    CONTINUE
C
        ZC1=Z*R/zzC1
        ZC2=Z*R/zzC2
        CALL BESSEL(ZC1,ZK01,ZK11,ZK21)
        CALL BESSEL(ZC2,ZK02,ZK12,ZK22)
        ZPSI=ZK02+(1.0D0/ZC2)*(ZK12-(zzC2/zzC1)*ZK11)
        ZCHI=ZK22-(zzC2**2/zzC1**2)*ZK21
C
        ZDPSI=-(1.0D0/R)*ZK02-(Z/zzC2+2.0D0/(ZC2*R))*ZK12+
     &         (zzC2**2/(zzC1**2*R))*
     &          ZK01+((2.0D0*zzC2**2)/(Z*zzC1*R**2))*ZK11
        ZDCHI=-(2.0D0/R)*ZK02-(Z/zzC2+4.0D0/(ZC2*R))*ZK12+
     &         ((2.0D0*zzC2**2)/
     &         (zzC1**2*R))*ZK01+((zzC2**2*Z)/zzC1**3+(4.0D0*zzC2**2)/
     &         (Z*zzC1*R**2))*ZK11
C
        ZE1=ZDPSI-ZCHI/R
        ZE2=2.0D0*ZCHI/R
        ZE3=2.0D0*ZDCHI*DRN
        ZE4=(zzC1**2/zzC2**2-2.0D0)*(ZDPSI-ZDCHI-ZCHI/R)
C
        DO 50 I=1,2
        DO 50 J=1,2
c                        ##########################
c                        !!!!!!!!!!!!!!!!!!!!!!!!!!
c                 modification for viscoelesticity
c        ZU(I,J)=(.5D0/(PI*RMU))*(ZPSI*D(I,J)-ZCHI*RC(I)*RC(J))
        ZU(I,J)=(.5D0/(PI*ZMU))*(ZPSI*D(I,J)-ZCHI*RC(I)*RC(J))
        ZT(I,J)=(.5D0/PI)*(ZE1*(DRN*D(I,J)+RC(J)*RN(I))-ZE2*
     &        (RN(J)*RC(I)-2.D0*RC(I)*RC(J)*DRN)-ZE3*RC(I)*
     &        RC(J)+ZE4*RC(I)*RN(J))
  50    CONTINUE
        RETURN
        END
c
C------------------------------------------------
        SUBROUTINE BESSEL(Z,ZK0,ZK1,ZK2)
C
C       OTHER EXPANSIONS ARE POSSIBLE AS WELL
C
C       Z IS THE BESSEL FUNCTION ARGUMENT
C
        IMPLICIT COMPLEX*16 (Z)
        IMPLICIT REAL*8 (A-H,O-Y)
C
        ZS=8.0D0*Z
        ZY=CDSQRT(Z)*CDEXP(Z)
        SIMAG=DIMAG(Z)
        SREAL=DREAL(Z)
        TOL=1.0D-09
        G=0.5772156649D0
        H=1.253314137D0
C
        IF (CDABS(Z) .LE. TOL) GO TO 10
        IF (CDABS(Z) .LE. 3.5D0) GO TO 20
        IF (CDABS(Z) .GT. 100.0D0) GO TO 40
        IF (SIMAG .GT. 0.0D0 .AND. DABS(SREAL) .LT. TOL) GO TO 40
        GO TO 30
C
C       BESSEL FUNCTION K0
C
  10    ZK0=-CDLOG(Z)
        GO TO 11
  20    ZI0=1.0D0+Z**2/4.0D0+Z**4/64.0D0+Z**6/2304.0D0+
     &      Z**8/147456.0D0+Z**10/14745600.0D0+Z**12/2123366400.0D0
        ZK0=-(CDLOG(Z/2.0D0)+G)*ZI0+Z**2/4.0D0+(1.5D0*Z**4)/
     &      64.0D0+(1.8333333333D0*Z**6)/2304.0D0+(2.0833333333D0*
     &      Z**8)/147456.0D0+(2.2833333333D0*Z**10)/14745600.0D0+
     &      (2.45D0*Z**12)/2123366400.0D0
        GO TO 21
  30    ZK0=H*(1.0D0-1.0D0/ZS+9.0D0/(2.0D0*ZS**2)-
     &     225.0D0/(6.0D0*ZS**3)+11025.0D0/(24.0D0*ZS**4)-
     &     893025.0D0/(120.0D0*ZS**5)+108056025.0D0/(720.0D0*
     &     ZS**6))/ZY
        GO TO 31
  40    ZK0=DCMPLX(0.0D0,0.0D0)
        GO TO 41
C       BESSEL FUNCTION K1
C
  11    ZK1=1.0D0/Z
        GO TO 50
  21    ZI1=(Z/2.0D0)*(1.0D0+Z**2/8.0D0+Z**4/192.0D0+Z**6/9216.0D0+
     &       Z**8/737280.0D0+Z**10/88473600.0D0+
     &       Z**12/14863564800.0D0)
        ZK1=1.0D0/Z+CDLOG(Z/2.0D0)*ZI1-(Z/4.0D0)*(-.1544313298D0+
     &      (1.345568670D0*Z**2)/8.0D0+(2.178902003D0*Z**4)/
     &      192.0D0+(2.762235336*Z**6)/9216.0D0+
     &      (3.212235336D0*Z**8)/
     &      737280.0D0+(3.578902003*Z**10)/88473600.0D0+
     &      (3.888425813D0*Z**12)/14863564800.0D0)
        GO TO 50
  31    ZK1=H*(1.0D0+3.0D0/ZS-15.0D0/(2.0D0*ZS**2)+315.0D0/
     &       (6.0D0*ZS**3)-14175.0D0/(24.0D0*ZS**4)+
     &       1091475.0D0/(120.0D0*ZS**5)-127702575.0D0/
     &       (720.0D0*ZS**6))/ZY
        GO TO 50
  41    ZK1=DCMPLX(0.0D0,0.0D0)
C
  50    CONTINUE
C
C       BESSEL FUNCTION K2
C
        ZK2=(2.0D0/Z)*ZK1+ZK0
        RETURN
        END
C------------------------------------------------
        SUBROUTINE ZERNEL(Y1,Y2,X1,X2,RN,SU,ST)
C
C       2-D ELASTOSTATIC KERNELS
C
        IMPLICIT REAL*8 (A-H,O-Y)
        DIMENSION RN(2),SU(2,2),ST(2,2),D(2,2),RC(2)
c                 modification for viscoelesticity
        COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM,fvisc,nmode
c        DATA D/1.0D0,0.0D0,0.0D0,1.0D0/
c       Kronecker delta
        d(1,1)=1.0d0
        d(2,1)=0.0d0
        d(1,2)=0.0d0
        d(2,2)=1.0d0

C
c                        ##########################
c                        !!!!!!!!!!!!!!!!!!!!!!!!!!
c        ATTENTION:  REAL NUMBERS, CALCULATION WITH RMU INSTEAD OF ZMU
c                 modification for viscoelesticity
c       Kelvin solid
c        if(nmode.eq.1)zmu=rmu+rmu*fvisc*z
c       Maxwell fluid
c        if(nmode.eq.2)zmu=rmu*z/(z+1.0d0/fvisc)
c
        CALL CLRARR(SU,4)
        CALL CLRARR(ST,4)
C
        PI=3.14159265359D0
        R=DSQRT((X1-Y1)**2+(X2-Y2)**2)
        RC(1)=(X1-Y1)/R
        RC(2)=(X2-Y2)/R
        DRN=RC(1)*RN(1)+RC(2)*RN(2)
C
        DO 1 I=1,2
        DO 1 J=1,2
        SU(I,J)=(-(3.0D0-4.0D0*RNU)*DLOG(R)*D(I,J)+RC(I)*RC(J))/
     &              (8.0D0*PI*RMU*(1.0D0-RNU))
c        SU(I,J)=(-(3.0D0-4.0D0*RNU)*DLOG(R)*D(I,J)+RC(I)*RC(J))/
c     &              (8.0D0*PI*ZMU*(1.0D0-RNU))
        ST(I,J)=-(DRN*((1.0D0-2.0D0*RNU)*D(I,J)+2.0D0*RC(I)*RC(J))-
     &          (1.0D0-2.0D0*RNU)*(RC(I)*RN(J)-RC(J)*RN(I)))/
     &          (4.0D0*PI*(1.0D0-RNU)*R)
  1     CONTINUE
        RETURN
        END
C---------------------------------------------------------------------
C
c        SUBROUTINE GENR(NCOUNT,KONN,MBCU,ZGM,ZFM,ZLOAD,ZTEMP,
c     &              N1,N2,n3,n4,SFDIAG,lbcu,ndtype,icffs,linv,minv)
        SUBROUTINE GENR(NCOUNT,KONN,MBCU,zlrar,ZLOAD,ZTEMP,
     &              N1,N2,n3,n4,SFDIAG,lbcu,ndtype,icffs,linv,minv)
c
c
C       READING OF LOADING AND DEFORMATION OF SYSTEM MATRICES
C
        IMPLICIT COMPLEX*16 (Z)
        IMPLICIT REAL*8 (A-H,O-Y)
        DIMENSION KONN(N3,N4),MBCU(N2)
c        DIMENSION ZGM(N2,N2),ZFM(N2,N2)
        dimension zlrar(n2,n2)


        dimension ZLOAD(N2),ZTEMP(N2)
        DIMENSION ZKU(3,2,2),ZKT(3,2,2),DJUMP(2,2)
        dimension sfdiag(n1,2,2)
        dimension lbcu(n3),ndtype(n1)
        dimension linv(n2),minv(n2)
c        dimension linv(n2+2*ndzon),minv(n2+2*ndzon)
        COMMON/CONTR/NNODE,NELEM,NP
c                 modification for viscoelesticity
        COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM,fvisc,nmode
        common/free/nffs,nfre
C
cc       jump term for infinity
c        DATA DJUMP1/1.0D0,0.0D0,0.0D0,1.0D0/
c       jump term for infinity
        djump(1,1)=1.0d0
        djump(2,1)=0.0d0
        djump(1,2)=0.0d0
        djump(2,2)=1.0d0  
cc       jump term for semi-infinity     
c        if(nfre.eq.1) then
c          do 91 i=1,2
c            do 91 j=1,2
c              djump(i,j)=0.5d0*djump(i,j)
c  91      continue
c        endif




        N51=N2*N2
c        CALL CLRZRR(ZGM,N51)
c        CALL CLRZRR(ZFM,N51)
        CALL CLRZRR(zlrar,N51)





C       DATA DJUMP/1.0D0,0.0D0,0.0D0,1.0D0/
C
C       READ AND WRITE LOAD IN LAPLACE TRANSFORM DOMAIN
C
c        if(nffs.eq.1) then
cc          rewind(9,*)
cc          read(9,*) zload
c          read(3,*)zload
c        endif






c        READ(1,*) ZLOAD
c        READ(3,*) ZLOAD
  43    continue
cc       take the loads on the actual free surface nodes zero
c        do 68 ii=1,n1
cc          ln2=nonn(ii,2)
cc          jref=2*(ln2-1)
c          jref=2*(ii-1)
c          if(ndtype(ii).eq.2) then
c            zload(jref+1)=dcmplx(0.0d0,0.0d0)
c            zload(jref+2)=dcmplx(0.0d0,0.0d0)
c          endif
c  68    continue            









        rewind 9
        WRITE(9,*) (ZLOAD(I),I=1,N2)
        WRITE(6,100) NCOUNT
        WRITE(6,103) (ZLOAD(I),I=1,N2)
C
        NOD=3
        GSCALE=RMU
c        REWIND 7
        REWIND 8
        DO 21 ICL=1,NNODE
        IROW=2*(ICL-1)
        DO 31 JCL=1,NELEM
        DO 41 K1=1,NOD
        LN=KONN(JCL,K1)
        JCOL=2*(LN-1)
c        READ(7,*) ((ZKU(K1,K2,K3),K2=1,2),K3=1,2)
        READ(8,*) ((ZKT(K1,K2,K3),K2=1,2),K3=1,2)
        DO 61 I=1,2
        DO 61 J=1,2

c        ZFM(IROW+I,JCOL+J)=ZFM(IROW+I,JCOL+J)+
c     &          ZKT(K1,I,J)
        zlrar(IROW+I,JCOL+J)=zlrar(IROW+I,JCOL+J)+
     &          ZKT(K1,I,J)
  61     CONTINUE
  41     CONTINUE
  31     CONTINUE
C
C       ADJUST DIAGONAL TERMS USING CONCEPT OF SUMMAT. OF OFF-DIAGONALS
C       DJUMP ARE JUMP TERMS FOR EXTERIOR TYPE PROBLEM
C
        DO 92 I=1,2
        DO 92 J=1,2
c        ZFM(IROW+I,IROW+J)=ZFM(IROW+I,IROW+J)+
c     &     SFDIAG(ICL,I,J)+DJUMP(I,J)
c        if(ndtype(icl).eq.2) then
c          zlrar(IROW+I,IROW+J)=zlrar(IROW+I,IROW+J)+
c     &       SFDIAG(ICL,I,J)+0.5d0*DJUMP(I,J)
cc          zlrar(IROW+I,IROW+J)=zlrar(IROW+I,IROW+J)+
ccc     &       SFDIAG(ICL,I,J)+0.0d0*DJUMP(I,J)
ccc
ccc
c          go to 92
c        endif




















        zlrar(IROW+I,IROW+J)=zlrar(IROW+I,IROW+J)+
     &     SFDIAG(ICL,I,J)+DJUMP(I,J)
  92     CONTINUE



  21     CONTINUE
        rewind 16
        write(16,*) zlrar


        CALL CLRZRR(zlrar,N51)
        REWIND 7
c        REWIND 8
        DO 2 ICL=1,NNODE
        IROW=2*(ICL-1)
        DO 3 JCL=1,NELEM
        DO 4 K1=1,NOD
        LN=KONN(JCL,K1)
        JCOL=2*(LN-1)
        READ(7,*) ((ZKU(K1,K2,K3),K2=1,2),K3=1,2)
c        READ(8,*) ((ZKT(K1,K2,K3),K2=1,2),K3=1,2)
        DO 6 I=1,2
        DO 6 J=1,2
ccccccccccccccccccccccccccccccccccccccccccccccccccc
cc       case of free surface: ignore the displacement kernel G (ZU) for
cc       the nodes and elements of the free surface, during the collocation
c        if((mbcu(2*icl).eq.2).or.(lbcu(jcl).eq.1)) go to 50
ccc

c       for free surface elements in the 2nd loop the displacement
c       kernel is 0, because it is multiplyed by the zero tractions(with
c       opposite sign) on the free surface.
c       free surface node or element
c        if((ndtype(icl).eq.2).or.(lbcu(jcl).eq.1)) then
c        if ((icffs.eq.2).and.(lbcu(jcl).eq.2)) go to 50






c          ZGM(IROW+I,JCOL+J)=ZGM(IROW+I,JCOL+J)+
c     &            ZKU(K1,I,J)*GSCALE
          zlrar(IROW+I,JCOL+J)=zlrar(IROW+I,JCOL+J)+
     &            ZKU(K1,I,J)*GSCALE

c  50    continue
c  60    continue
  6     CONTINUE
  4     CONTINUE
  3     CONTINUE
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
c        ZLOAD(I)=ZLOAD(I)
  10    CONTINUE
c        CALL BOUND(ZGM,ZFM,ZLOAD,ZTEMP,MBCU,linv,minv,N2)
c        CALL BOUND(ZGM,ZFM,ZLOAD,ZTEMP,MBCU,N2,linv,minv)
        CALL BOUND(zlrar,ZLOAD,ZTEMP,MBCU,N2,linv,minv)



C
C       SORT OUT DISPLACEMENTS AND TRACTIONS-STORE IN ZTEMP,ZLOAD RESP.
C
        DO 11 I=1,N2
        ZTEMP(I)=ZLOAD(I)
  11    CONTINUE
        REWIND 9
        READ(9,*) (ZLOAD(I),I=1,N2)
        REWIND 9
C
        DO 12 I=1,N2
        LC=MBCU(I)
        IF (LC .EQ. 0) GO TO 12
cccccccccccccccccccccccccccccccccccccccccccccccccc
c        IF (lc .EQ. 2) GO TO 12
c        IF (lc .EQ. 3) GO TO 12
c
        ZSTUMP=ZTEMP(I)
        ZTEMP(I)=ZLOAD(I)
        ZLOAD(I)=ZSTUMP*GSCALE
c        ZLOAD(I)=ZSTUMP





  12    CONTINUE
C
C       PRINTING OF TRACTIONS AND DISPLACEM. IN LAPLACE TRANSFORM DOMAIN
C
cc       save the displacements for use in ploting in postbem
c        write(11,*) ztemp
c        in subregion case of calculating the free field write
c         if((nffs.eq.2).and.(icffs.eq.1)) then
         if((nffs.gt.1).and.(icffs.eq.1)) then
           write(6,110)
           go to 49
         endif
         if((nffs.eq.3).and.(icffs.eq.2)) then
           write(6,115)
           go to 49
         endif
         write(6,120)
  49     continue
        WRITE(6,101) NCOUNT
        WRITE(6,103) (ZLOAD(I),I=1,N2)
        WRITE(6,102) NCOUNT
        WRITE(6,103) (ZTEMP(I),I=1,N2)
  100   FORMAT(//2X,' LOADING AT STEP=',I3)
  101   FORMAT(//2X,' TRACTIONS AT STEP=',I3)
  102   FORMAT(//2X,' DISPLACEMENTS AT STEP=',I3)
  103   FORMAT(8(2X,E12.6))
  110   format(//2x,' FREE FIELD CALCULATIONS (1st STEP)')
  115   format(//2x,' FREE FIELD CALCULATIONS (2nd STEP)')
  120   format(//2x,' SCATTERED FIELD CALCULATIONS ')
        RETURN
        END
C---------------------------------------------------------------------
c        SUBROUTINE GENRleft(NCOUNT,KONN,konzon,konnot,MBCU,zlrar,
c     &                N1,N2,n3,n4,ndzon,SFDIAG,lbcu,ndtype,icffs)
c        SUBROUTINE GENRleft(NCOUNT,KONN,konzon,konnot,nonzon,nonnot,
c     &        MBCU,zlrar,N1,N2,n3,n4,ndzon,SFDIAG,lbcu,ndtype,icffs)
        SUBROUTINE GENRleft(NCOUNT,KONN,konzon,nonzon,
     &        MBCU,zlrar,N1,N2,n3,n4,ndzon,SFDIAG,lbcu,ndtype,icffs)


c
c       forms the left array zlrar
c
C       READING OF LOADING AND DEFORMATION OF SYSTEM MATRICES
C
        IMPLICIT COMPLEX*16 (Z)
        IMPLICIT real*8 (a-h,o-y)
        DIMENSION KONN(N3,N4),MBCU(N2)
        dimension konzon(n3,n4)
c              ,konnot(n3,n4)
        dimension nonzon(n1,n4)
c              ,nonnot(n1,n4)


        dimension zlrar(n2+2*ndzon,n2+2*ndzon)
        DIMENSION ZKU(3,2,2),ZKT(3,2,2),DJUMP(2,2)
        dimension sfdiag(n1,2,2)
        dimension lbcu(n3),ndtype(n1)



        COMMON/CONTR/NNODE,NELEM,NP
        COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM,fvisc,nmode
        common/free/nffs,nfre
C
cc       jump term for infinity
c        DATA DJUMP1/1.0D0,0.0D0,0.0D0,1.0D0/
c       jump term for infinity
        djump(1,1)=1.0d0
        djump(2,1)=0.0d0
        djump(1,2)=0.0d0
        djump(2,2)=1.0d0

c        djump1(1,1)=0.5d0
c        djump1(2,1)=0.0d0
c        djump1(1,2)=0.0d0
c        djump1(2,2)=0.5d0

cc       jump term for semi-infinity
c        if(nfre.eq.1) then
c          do 91 i=1,2
c            do 91 j=1,2
c              djump(i,j)=0.5d0*djump(i,j)
c  91      continue
c         endif









        N51=N2*N2
        nt=n2+2*ndzon
        nt2=nt*nt
        call clrzrr(zlrar,nt2)
        NOD=3
        GSCALE=RMU
c        REWIND 7
        REWIND 8
        DO 2 ICL=1,NNODE
        IROW=2*(ICL-1)
        DO 3 JCL=1,NELEM
        DO 4 K1=1,NOD
        LN=KONN(JCL,K1)
c        LN=KONN(JCL,nod+1-K1)


        JCOL=2*(LN-1)
c        READ(7,*) ((ZKU(K1,K2,K3),K2=1,2),K3=1,2)
        READ(8,*) ((ZKT(K1,K2,K3),K2=1,2),K3=1,2)
        DO 6 I=1,2
        DO 6 J=1,2
ccccccccccccccccccccccccccccccccccccccccccccccccccc
cc       case of free surface: ignore the displacement kernel G (ZU) for
cc       the nodes and elements of the free surface, during the collocation
         zlrar(IROW+I,JCOL+J)=zlrar(IROW+I,JCOL+J)+
     &          ZKT(K1,I,J)
c  60    continue
  6     CONTINUE
  4     CONTINUE
  3     CONTINUE
C
C       ADJUST DIAGONAL TERMS USING CONCEPT OF SUMMAT. OF OFF-DIAGONALS
C       DJUMP ARE JUMP TERMS FOR EXTERIOR TYPE PROBLEM
C


c        if(ndtype(icl).eq.3) go to 2





        DO 9 I=1,2
        DO 9 J=1,2
c        if(ndtype(icl).eq.2) then
c          zlrar(IROW+I,IROW+J)=zlrar(IROW+I,IROW+J)+
c     &       SFDIAG(ICL,I,J)+0.5d0*DJUMP(I,J)
ccc          zlrar(IROW+I,IROW+J)=zlrar(IROW+I,IROW+J)+
ccc     &       SFDIAG(ICL,I,J)+0.0d0*DJUMP(I,J)
ccccc
ccccc
ccccc
c           go to 9
c        endif

















c
c
c
cc        if(ndtype(icl).eq.2) go to 27
          zlrar(IROW+I,IROW+J)=zlrar(IROW+I,IROW+J)+
     &       SFDIAG(ICL,I,J)+DJUMP(I,J)




  9     CONTINUE
  2     CONTINUE
c
c       part of displacement kernel matrix, that multiplys the
c       tractions on the subregion only 
c       (-G12, -G22)
c        isim1=0
        REWIND 7
        DO 23 ICL=1,NNODE
c          if(ndtype(icl).eq.3) then


            IROW=2*(ICL-1)
c            kdiak=0
            DO 33 JCL=1,NELEM
              DO 47 K1=1,NOD
                READ(7,*) ((ZKU(K1,K2,K3),K2=1,2),K3=1,2)
c               exclude the terms of the submatrix that maltiply tractions on
c               non-sudregion boundary
                if((lbcu(jcl).eq.1).or.(lbcu(jcl).eq.3)) then
                  ln=konzon(jcl,k1)
c                  ln=konzon(jcl,nod+1-k1)
                  jcol=n2+2*(ln-1)


                  DO 63 I=1,2
                    DO 63 J=1,2
                     zlrar(IROW+I,JCOL+J)=zlrar(IROW+I,JCOL+J)-
     &                    ZKU(K1,I,J)*GSCALE


  63              CONTINUE
                endif
  47          CONTINUE 
  33        CONTINUE

  23    continue

        return
        end
c
C---------------------------------------------------------------------
C
        SUBROUTINE GENRzone(NCOUNT,KONN,konzon,nonzon,MBCU,zlrar,
     &            N1,N2,n3,n4,ndzon,SFDIAG,lbcu,ndtype,icffs)
c
c
C       READING OF LOADING AND DEFORMATION OF SYSTEM MATRICES
C
        IMPLICIT COMPLEX*16 (Z)
        IMPLICIT real*8 (a-h,o-y)
        DIMENSION KONN(N3,N4),MBCU(N2)
        dimension konzon(n3,n4)
c        dimension sfdiagzone(ndzon,2,2)
c        dimension sfdiag(ndzon,2,2)
        dimension sfdiag(n1,2,2)



        dimension nonzon(n1,n4)
        DIMENSION ZKU(3,2,2),ZKT(3,2,2),DJUMP(2,2)
c        dimension sfdiag(n1,2,2)
        dimension lbcu(n3),ndtype(n1)
c        DIMENSION ZGMz(ndzon,ndzon),ZFMz(ndzon,ndzon)
c        DIMENSION ZGMz(n2,n2),ZFMz(n2,n2)
c        dimension zgm2(n2,n2)
        dimension zlrar(n2+2*ndzon,n2+2*ndzon)
c                         ##############################
c                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c        dimension linv(n2),minv(n2)
C                END MOPDIF  %%%%%%%%%%%%%%%%%
c        COMMON/CONTR/NNODE,NELEM,NTINT,NP
        COMMON/CONTR/NNODE,NELEM,NP
c                        ##########################
        COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM,fvisc,nmode
        common/free/nffs,nfre
C
c       jump term for infinity


        djump(1,1)=0.0d0
        djump(2,1)=0.0d0
        djump(1,2)=0.0d0
        djump(2,2)=0.0d0




cc       jump term for semi-infinity
c        if(nfre.eq.1) then
c          do 91 i=1,2
c            do 91 j=1,2
c              djump(i,j)=0.5d0*djump(i,j)
c  91      continue
c        endif


C
c        DATA DJUMP/1.0D0,0.0D0,0.0D0,1.0D0/
c        DATA DJUMP/0.5D0,0.0D0,0.0D0,0.5D0/
c        DATA DJUMP/0.0D0,0.0D0,0.0D0,0.0D0/
        N51=N2*N2
c        nz=n2-2*ndzon


        NOD=3
        GSCALE=RMU
        REWIND 15
c       the traction kernels first
c        calculate the matrix F
c        isim1=0
        DO 2 ICL=1,NNODE
          kdiak=0
c


c          if((nffs.eq.3).and.((icffs.eq.1).or.(icffs.eq.2)).and.
          if((ndtype(icl).eq.1).or.(ndtype(icl).eq.3)) then

c            isim1=isim1+1
c            irow=n2+2*(icl-1)
c            irow=n2+2*(isim1-1)
            ln1=nonzon(icl,2)
            irow=n2+2*(ln1-1)
            DO 3 JCL=1,NELEM
c             collect the elements of the zone or subregion only
              if((lbcu(jcl).eq.1).or.(lbcu(jcl).eq.3)) then
c                kdiak=1



                DO 4 K1=1,NOD
c                  READ(14,*) ((ZKU(K1,K2,K3),K2=1,2),K3=1,2)
                  READ(15,*) ((ZKT(K1,K2,K3),K2=1,2),K3=1,2)
c                 change the connectivity of the nodes of the element, in order
c                 to obtain inward normal unit vector (zone)
                   LN=KONN(JCL,nod+1-K1)

c                    ln=konzon(jcl,k1)
c                    jcol=n2-2*ndzon+2*(ln-1)

                    jcol=2*(ln-1)

                    DO 6 I=1,2
                      DO 6 J=1,2
                        zlrar(IROW+I,JCOL+J)=zlrar(IROW+I,JCOL+J)+
     &                        ZKT(K1,I,J)
  6                 CONTINUE
  4             CONTINUE
              endif
  3         CONTINUE
C
C       ADJUST DIAGONAL TERMS USING CONCEPT OF SUMMAT. OF OFF-DIAGONALS
C       DJUMP ARE JUMP TERMS FOR inTERIOR TYPE PROBLEM (=0)
C
c             ln2=nonzon(icl,2)
c             ln2=nonzon(isim1,2)
c             jcol=n2-2*ndzon+2*(ln2-1)
c             jcol=2*(ln1-1)
c             jcol=nz+2*(ln1-1)
             jcol=2*(icl-1)


c            if(kdiak.eq.1) then
               DO 9 I=1,2
                 DO 9 J=1,2
c                   zlrar(IROW+I,jcol+J)=zlrar(IROW+I,jcol+J)+
c     &               SFDIAG(icl,I,J)+DJUMP(I,J)
                   zlrar(IROW+I,jcol+J)=zlrar(IROW+I,jcol+J)+
     &               SFDIAG(icl,I,J)-DJUMP(I,J)
c                   zlrar(IROW+I,jcol+J)=zlrar(IROW+I,jcol+J)+
c     &               SFDIAG(icl,I,J)-DJUMP(I,J)














c                   zlrar(IROW+I,jcol+J)=zlrar(IROW+I,jcol+J)-
c     &               SFDIAG(icl,I,J)+DJUMP(I,J)







  9           CONTINUE
c            endif
























          endif
  2     continue
c
c       the displacement kernels next
c       calculate the matrix G
        REWIND 14
c        isim1=0
        DO 21 ICL=1,NNODE
c          if((nffs.eq.3).and.((icffs.eq.1).or.(icffs.eq.2)).and.
          if((ndtype(icl).eq.1).or.(ndtype(icl).eq.3)) then


c            isim1=isim1+1
c            irow=n2+2*(icl-1)
c            irow=n2+2*(isim1-1)
            ln1=nonzon(icl,2)
            irow=n2+2*(ln1-1)
            DO 31 JCL=1,NELEM
c             collect the elements of the zone or subregion only
              if((lbcu(jcl).eq.1).or.(lbcu(jcl).eq.3)) then

                DO 41 K1=1,NOD
                  READ(14,*) ((ZKU(K1,K2,K3),K2=1,2),K3=1,2)
c                   READ(15,*) ((ZKT(K1,K2,K3),K2=1,2),K3=1,2)
c                 take only the elements with lbcu=1(inclusion case),
c                                    or lbcu=3(simulation case)
c                 collect the elements of the subregion only 
c                  if(lbcu(jcl).eq.3) then 
c                  ln=konzon(jcl,k1)
c                  jcol=n2+2*(ln-1)
c                  ln=nonzon(icl,2)
c                  ln=konzon(jcl,k1)

c                 change the connectivity of the nodes of the element, in order
c                 to obtain inward normal unit vector (zone)
                  ln=konzon(jcl,nod+1-k1)   


                  jcol=n2+2*(ln-1)
c                  JCOL=2*(LN-1)
                  DO 61 I=1,2
                    DO 61 J=1,2
                      zlrar(IROW+I,JCOL+J)=zlrar(IROW+I,JCOL+J)+
     &                    ZKU(K1,I,J)*gscale
c                      zlrar(IROW+I,JCOL+J)=zlrar(IROW+I,JCOL+J)-
c     &                    ZKU(K1,I,J)*gscale


  61              CONTINUE
c                  endif
  41            CONTINUE
              endif
  31        CONTINUE
          endif
  21    CONTINUE
c
c       save the array zlrar in file 16
        rewind 16
        write(16,*) zlrar
        RETURN
        END
C---------------------------------------------------------------------
C
c        SUBROUTINE GENRright(NCOUNT,KONN,konzon,konnot,nonn,nonzon,
c     &       MBCU,zlrar,ZLOAD,ZTEMP,zloadz,ztempz,N1,N2,n3,n4,ndzon,
c     &       lbcu,ndtype,icffs,linv,minv)
        SUBROUTINE GENRright(NCOUNT,KONN,konzon,nonn,nonzon,
     &       MBCU,zlrar,ZLOAD,ZTEMP,zloadz,ztempz,N1,N2,n3,n4,ndzon,
     &       lbcu,ndtype,icffs,linv,minv)


c
c
C       READING OF LOADING AND DEFORMATION OF SYSTEM MATRICES
C
        IMPLICIT COMPLEX*16 (Z)
        IMPLICIT real*8 (a-h,o-y)
        DIMENSION KONN(N3,N4),MBCU(N2)
        DIMENSION ZLOAD(N2),ZTEMP(N2) 
        DIMENSION ZKU(3,2,2)
c             ,DJUMP(2,2)
        dimension konzon(n3,n4)
c             ,konnot(n3,n4)
c        dimension sfdiag(n1,2,2)
        dimension lbcu(n3),ndtype(n1)
        dimension zlrar(n2+2*ndzon,n2+2*ndzon)
        dimension zloadz(n2+2*ndzon),ztempz(n2+2*ndzon)
        dimension linv(n2+2*ndzon),minv(n2+2*ndzon)
        dimension nonn(n1,n4)
        dimension nonzon(n1,n4)
        COMMON/CONTR/NNODE,NELEM,NP
        COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM,fvisc,nmode
        common/free/nffs,nfre
C
c        DATA DJUMP/0.0D0,0.0D0,0.0D0,0.0D0/
        N51=N2*N2
        nt=n2+2*ndzon
        nt2=nt*nt
        call clrzrr(zlrar,nt2)
        call clrzrr(zloadz,nt)
        NOD=3 
        GSCALE=RMU
        REWIND 7
c
c       the displacement kernels 
c        isim1=0
        DO 2 ICL=1,NNODE
c        in subregion case of calculating the free field(2nd loop),take the
c        parts of the displacement kernel, for the formation of the right 
c        hand matrix of the equation
c              G11,G21 : parts of the non inclusion boundary
           irow=2*(icl-1)
           DO 3 JCL=1,NELEM
               DO 4 K1=1,NOD 
                 READ(7,*) ((ZKU(K1,K2,K3),K2=1,2),K3=1,2)
c                 READ(8,*) ((ZKT(K1,K2,K3),K2=1,2),K3=1,2)
c                elements with lbcu=1(inclusion case),
c                           or lbcu=3(simulation case)
c                exclude the terms of the submatrix that maltiply tractions on
c                sudregion boundary
c                 if(lbcu(jcl).eq.3)then
                 if(lbcu(jcl).ne.3)then

                   LN=KONN(JCL,K1)
                   jcol=2*(ln-1)
                   DO 6 I=1,2
                     DO 6 J=1,2

                       zlrar(IROW+I,JCOL+J)=zlrar(IROW+I,JCOL+J)+
     &                     ZKU(K1,I,J)*gscale


  6                CONTINUE
                 endif

  4         CONTINUE
  3     CONTINUE

  2    CONTINUE
c
c       traction array (only for the non-zone nodes)
c       read the incident field
c        rewind 3
        READ(3,*) ZLOAD
        rewind 9
        write(9,*) zload
c       loads on the nonzone nodes
c        z1=dcmplx(0.0d0,-1.0d0)
        do 89 i=1,n1
          ln2=nonn(i,2)
          jref=2*(ln2-1)
c          zloadz(jref+1)=z1*zload(jref+1)
c          zloadz(jref+2)=z1*zload(jref+2)
          zloadz(jref+1)=zload(jref+1)
          zloadz(jref+2)=zload(jref+2)
c          zloadz(jref+1)=-zload(jref+1)
c          zloadz(jref+2)=-zload(jref+2)






































          if(ndtype(i).eq.3) then
            zload(jref+1)=dcmplx(0.0d0,0.0d0)
            zload(jref+2)=dcmplx(0.0d0,0.0d0)
            zloadz(jref+1)=dcmplx(0.0d0,0.0d0)
            zloadz(jref+2)=dcmplx(0.0d0,0.0d0)
          endif
  89    continue
c        write(4,*)zload
        do 67 ii=1,n2 

C
C         SOLVE FOR UNKNOWN DISPLACEM. AND/OR LOADS-WATCH SCALING FACTOR
C
          LC=MBCU(ii)
          IF (LC .EQ. 1) GO TO 10
          ZLOADz(ii)=ZLOADz(ii)/GSCALE
c          ZLOADz(ii)=ZLOADz(ii)


  10      CONTINUE
  67    continue
c       form the right hand side of the equation
        nka=n2+2*ndzon
        call clrzrr(ztempz,nka)
        CALL ZVPROD(zlrar,ZLOADz,NKA,ztempz)
c
c       solve the equation
        call clrzrr(zlrar,nt2)
c       read the left hand side array zlrar from file 16
        rewind 16
        read(16,*) zlrar
c       invert this
C        CALL ZMINV(zlrar,nt2)
        CALL ZMINV(zlrar,nt,linv,minv)
c        ndz=n2+2*ndzon
c        CALL ZMINV(zlrar,nt,linv,minv,ndz)
c       finally multiply the inverted zlrar with the right hand side ztempz
        call clrzrr(zloadz,nt)
        CALL ZVPROD(zlrar,ztempz,NKA,zloadz)
c       take as loads in matrix zload, the incident loads on the
c       non-subregion boundary nodes as 0, plus the free field loads
c       (with opposite sign) on the subregion boundary.

c       find the free field tractions (only on zone nodes)
        call clrzrr(zload,n2)
        do 75 i=1,n1

          if(ndtype(i).eq.3) then
            ln2=nonzon(i,2)
            jref1=2*(ln2-1)
            ln2=nonn(i,2)
            jref2=2*(ln2-1)

            lc=mbcu(jref2+1)
            if(lc.eq.1) go to 71
            zload(jref2+1)=zloadz(n2+jref1+1)*gscale
            go to 72
  71        continue
            zload(jref2+1)=zloadz(n2+jref1+1)
  72        continue
            lc=mbcu(jref2+2)
            if(lc.eq.1) go to 73
            zload(jref2+2)=zloadz(n2+jref1+2)*gscale
            go to 74
  73        continue
            zload(jref2+2)=zloadz(n2+jref1+2)
  74        continue

          endif
  75    continue

        rewind 19
        write(19,*) zload
c       extract the dispalcements from the array zloadz
        call clrzrr(ztemp,n2)
        do 78 i=1,n1
          ln2=nonn(i,2)
          jref=2*(ln2-1)

          ztemp(jref+1)=zloadz(jref+1)
          ztemp(jref+2)=zloadz(jref+2)
cc         for the subregion nodes, put negative sign
c          if(ndtype(i).eq.3) then
c            ztemp(jref+1)=-zloadz(jref+1)
c            ztemp(jref+2)=-zloadz(jref+2)
c          endif

  78    continue

        do 92 i=1,n2
        lc=mbcu(i)
        if(lc.eq.0) go to 92
        ZSTUMP=ZTEMP(I)
        ZTEMP(I)=ZLOAD(I)
        ZLOAD(I)=ZSTUMP*GSCALE
  92  continue

c        rewind 20
c        write(20,*) zload
c        rewind 18
c        write(18,*) ztemp
        RETURN
        END
C------------------------------------------------
c        SUBROUTINE BOUND(ZGM,ZFM,ZLOAD,ZTEMP,MBCU,N2,linv,minv)
        SUBROUTINE BOUND(zlrar,ZLOAD,ZTEMP,MBCU,N2,linv,minv)
C
C       MATRIX RE-ARRANGEM. ACCORDING TO ROUNDARY CONDIT. AND SOLUTION
C
        IMPLICIT COMPLEX*16 (Z)
        IMPLICIT REAL*8 (A-H,O-Y)
        DIMENSION ZLOAD(N2),ZTEMP(N2),MBCU(N2)
c          ZGM(N2,N2),ZFM(N2,N2)



c                         ##############################
c                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        dimension linv(n2),minv(n2)
c        dimension linv(n2+2*ndzon),minv(n2+2*ndzon)
        dimension zlrar(n2,n2)
C
C       INTERCHANGE COLUMNS OF ZGM, ZFM MATRICES
C
        NKA=N2

c          ?????????????????????
c        there is not capability for interchange of columns
c        DO 10 J=1,NKA
c        IF (MBCU(J) .EQ. 0) GO TO 10
c        DO 11 I=1,NKA
c        ZTEMP(I)=ZFM(I,J)
c        ZFM(I,J)=-ZGM(I,J)
c        ZGM(I,J)=-ZTEMP(I)
c  11    CONTINUE
c  10    CONTINUE









C
C       SOLUTION FOR RESPONSE
C
c        CALL ZVPROD(ZGM,ZLOAD,NKA,ZTEMP)
        CALL ZVPROD(zlrar,ZLOAD,NKA,ZTEMP)
c       read ZFM
        rewind 16
        read(16,*) zlrar
c        CALL ZMINV(ZFM,NKA,linv,minv)
        CALL ZMINV(zlrar,NKA,linv,minv)
c        CALL ZVPROD(ZFM,ZTEMP,NKA,ZLOAD)
        CALL ZVPROD(zlrar,ZTEMP,NKA,ZLOAD)
        RETURN
        END
C------------------------------------------------
c        SUBROUTINE ZMINV(A,N)
        SUBROUTINE ZMINV(A,N,L,M)
C
C       COMPLEX MATRIX INVERSION --- no LIMIT, square MATRIX
C
        IMPLICIT COMPLEX*16 (A-H)
        IMPLICIT REAL*8 (R)
C        DIMENSION A(1),L(100),M(100)
c        DIMENSION A(1),L(100),M(100)
c        DIMENSION A(n*n),L(300),M(300)
        DIMENSION A(n*n),L(n),M(n)
c
        do 5 i=1,n
          l(n)=0
          m(n)=0
  5     continue
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
        RA=CDABS(CA)
        RB=CDABS(BIGA)
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
  80    IF (CDABS(BIGA)-1.0D-15) 90,100,100
  90    CONTINUE
        WRITE (NLIN,95) K
  95    FORMAT(20X,'DURING A MATRIX INVERSION OPERATION A ZERO',/
     &      20X,' DETERMINANT HAS BEEN DETECTED. THE RUN HAS
     &      BEEN',/' TERMINATED. ROW OR COLUMN NO. WAS ',I4//)
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
        SUBROUTINE ZVPROD(ZA,ZB,N,ZC)
C
C       PRODUCT OF COMPLEX MATRIX TIMES COMPLEX VECTOR
C
        IMPLICIT COMPLEX*16 (Z)
        DIMENSION ZA(N,N),ZB(N),ZC(N)
C
        DO 20 I=1,N
        ZC(I)=dcmplx(0.0d0,0.0d0)
        DO 10 J=1,N
        ZC(I)=ZC(I)+ZA(I,J)*ZB(J)
  10    CONTINUE
  20    CONTINUE
        RETURN
        END
C------------------------------------------------
c        SUBROUTINE HOOP(KONN,nonn,ZLOAD,ZTEMP,ZLOAD1,
c     &          XCR,YCR,N1,N2,N3,N4,N5,zz1,zz2,zz3,z,icffs,ndtype)
        SUBROUTINE HOOP(KONN,nonn,ZLOAD,ZTEMP,ZLOAD1,zloadf1,
     &          XCR,YCR,N1,N2,N3,N4,N5,zz1,zz2,zz3,z,icffs,ndtype)








C
C       HOOP STRESSES IN THE LAPLACE DOMAIN
C
        IMPLICIT COMPLEX*16 (Z)
        IMPLICIT REAL*8 (A-H,O-Y)
        DIMENSION KONN(N3,N4),ZLOAD(N2),ZTEMP(N2),XCR(N1),YCR(N1)
        DIMENSION ZMA(7,7),ZV(7),ZW(7)
        DIMENSION ZTR(2),ZDI(2)
        DIMENSION LDIR(3),XY(3,2),SPOINT(3),SHP(3),DSHP(3),RN(2),RQ(2)
c        dimension zz1(n1),zz2(n1),zz3(n1)
        dimension zz1(n4,n3),zz2(n4,n3),zz3(n4,n3)
        DIMENSION ZLOAD1(N5)
        DIMENSION ZLOADf1(N5)

        dimension nonn(n1,n4)
        dimension ndtype(n1)
        dimension linv(7),minv(7)
        COMMON/CONTR/NNODE,NELEM,NP
c                 modification for viscoelesticity
        COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM,fvisc,nmode
        COMMON/GAUSS/GPN(10,10),WTN(10,10),GPS(10,10),WTS(10,10)
        common/free/nffs,nfre
        DATA SPOINT/-1.0d0,0.0d0,1.0d0/
 
c
        NOD=3
cc                 modification for viscoelesticity
c        tol=1.0d-9
cc       Kelvin solid
c        if(dabs(fvisc).lt.tol) then
c          zmu=rmu
c          zla=rla
c          go to 5
c        endif
cc       Kelvin solid
c        if(nmode.eq.1) then 
cc          zmu=rmu+rmu*fvisc*z
cc          zla=rla-(2.0d0/3.0d0)*rmu*fvisc*z
c          zmu=rmu+fvisc*z
c          zla=rla-(2.0d0/3.0d0)*fvisc*z
c        endif
cc       Maxwell fluid
c        if(nmode.eq.2) then
cc          zmu=rmu*z/(z+1.0d0/fvisc)
cc          zla=(rla*z+(1.0d0/fvisc)*(rla+(2.0d0/3.0d0)*rmu))/
cc     &          (z+1.0d0/fvisc)
c          zmu=rmu*z/(z+rmu/fvisc)
c          zla=(rla*z+(rmu/fvisc)*(rla+(2.0d0/3.0d0)*rmu))/
c     &          (z+rmu/fvisc)
c        endif
c  5     continue
c

C              READ THE INCIDENT LOADS IN L.T.
c       in usual case
c        if(((nffs.eq.1).or.(nffs.eq.2)).or.
c     &     ((nffs.eq.3).and.(icffs.eq.1)))then




cccc      in the 1st loop
ccc       if(icffs.eq.1) then
cc      in the 2nd loop
c       if(icffs.eq.2) then

c       find the complex values of Lame's constants
        call viscoel(z,zmu,zla,znu,zc1,zc2)

        if(icffs.eq.1)then
          READ(2,*) ZLOAD1
c         save the incident hoop stresses for using them in the second loop
          rewind 4
          write(4,*) zload1
        endif


c
c
c
c
c          go to 43
c
c        endif


c
c
c       for the 2nd loop only (subregion or simulation case)
        if(((nffs.eq.2).or.(nffs.eq.3)).and.(icffs.eq.2)) then
          rewind 11
          read(11,*) zload1
c          read(11,*) zloadf1


          rewind 4
          read(4,*) zloadf1


c
c
        endif

 43    continue
c      DO 3 JCL=1,nnode
      DO 3 JCL=1,nelem
c        ndiak=0
c        ZSRR=dcmplx(0.0d0,0.0d0)
c        ZSTT=dcmplx(0.0d0,0.0d0)
c        ZSRT=dcmplx(0.0d0,0.0d0)
c       calculate the coordinates of the previus,current(jcl) and next node
c       index for quadratic interpolation
        nod=3
c
c       in the case that the number of the first (or last) node, is the
c       same with the number of the current node (e.g. edge nodes in 
c       open surfaces),
c        ln1=nonn(jcl,1)
c        ln2=nonn(jcl,2)
c        ln3=nonn(jcl,3)
c       in the 1st loop take only the nodes on the zone boundary
c        if((nffs.eq.3).and.(icffs.eq.1).and.(ndtype(jcl).ne.3))
c     &                    go to 79


cc       for the 1st loop and for the simulation nodes, or the free surface nodes,
cc       take the incident hoop stresses 0
c        if(((nffs.eq.3).and.(icffs.eq.1)).and.
c     &         ((ndtype(jcl).eq.3).or.(ndtype(jcl).eq.2)) then

ccc       for the 1st loop and for the free surface nodes,
ccc       take the incident hoop stresses 0
cc        if((nffs.eq.3).and.(icffs.eq.1).and.(ndtype(jcl).eq.2)) then
ccc          ndiak=1
c        if((nffs.eq.3).and.(icffs.eq.1).and.(ndtype(jcl).eq.3)) then
cc          ndiak=1
c          jref=3*(ln2-1)
c          ZLOAD1(JREF+1)=dcmplx(0.0d0,0.0d0)
c          ZLOAD1(JREF+2)=dcmplx(0.0d0,0.0d0)
c          ZLOAD1(JREF+3)=dcmplx(0.0d0,0.0d0)
c        endif
c        if((nffs.eq.3).and.(icffs.eq.1).and.(ndtype(jcl).eq.3)) then
cc          ndiak=1
cc          jref=3*(ln2-1)
c          jref=3*(ldir(i)-1)
c          ZLOAD1(JREF+1)=dcmplx(0.0d0,0.0d0)
c          ZLOAD1(JREF+2)=dcmplx(0.0d0,0.0d0)
c          ZLOAD1(JREF+3)=dcmplx(0.0d0,0.0d0)
c        endif

c        if((nffs.eq.3).and.(icffs.eq.2)) then
cc              .and.(ndtype(jcl).eq.3)) then
cc          ndiak=1
c          jref=3*(ln2-1)
cc          read(4,*) zl1,zl2,zl3
c          if(ndtype(jcl).eq.3) then
c            ZLOAD1(JREF+1)=zloadf1(jref+1)-zload1(jref+1)
c            ZLOAD1(JREF+2)=zloadf1(jref+2)-zload1(jref+2)
c            ZLOAD1(JREF+3)=zloadf1(jref+3)-zload1(jref+3)
c          endif
c          if(ndtype(jcl).eq.2) then
c            ZLOAD1(JREF+1)=-zloadf1(jref+1)
c            ZLOAD1(JREF+2)=-zloadf1(jref+2)
c            ZLOAD1(JREF+3)=-zloadf1(jref+3)
c          endif

c          if((ndtype(jcl).eq.2).and.(icffs.eq.2)) then
c            jref=3*(ln2-1)
c            ZLOADf1(JREF+1)=dcmplx(0.0d0,0.0d0)
c            ZLOADf1(JREF+2)=dcmplx(0.0d0,0.0d0)
c            ZLOADf1(JREF+3)=dcmplx(0.0d0,0.0d0)
c          endif

c
c
c        endif


cc
c        if((nffs.eq.3).and.(icffs.eq.1).and.(ndtype(jcl).eq.3))
c     &          ndiak=1
cc
cc       for simulation node,change the connectivity of nodes, 
cc       in order to obtain inward normal unit vector (zone)
c        if(ndiak.eq.1) then
c          lnt=ln1
c          ln1=ln3
c          ln3=lnt
c        endif
cc
        DO 4 I=1,NOD
cc         for simulation node,change the coordinates of the connected nodes, 
cc         in order to obtain inward normal unit vector (zone)
c          if(ndiak.eq.1) then
c            LN=nONN(JCL,nod+1-I)
c            go to 83
c          endif
c          LN=nONN(JCL,I)
          LN=KONN(JCL,I)
c  83      continue
c
          LDIR(I)=LN
          XY(I,1)=XCR(LN)
          XY(I,2)=YCR(LN)
  4     CONTINUE
c



cc       interpolation in the middle node only for quadratic elements (eks=0)
c        eks=spoint(2)
c
cc       case that first node is equal with the current
cc       index for linear interpolation
c        if(ln1.eq.ln2) then
c          nod=2
cc         interpolation in the first node for linear element
cc          eks=0.0d0
c          eks=spoint(1)
cc          ldir(2)=ln1
c
c        endif
cc       case that the last node is equal with the current
cc       index for linear interpolation
c        if(ln2.eq.ln3) then
c          nod=2
cc         interpolation in the 2nd node for linear element
cc          eks=1.0d0
c          eks=spoint(3)
cc          ldir(2)=ln3  
c        endif
cc       if nod=3: quadratic interpolation
cc          nod=2: linear       -//-
cc       in case of linear interpolation, take the coordinates of the 
cc       two edge nodes (previus and next, not middle)
c        if(nod.eq.2) then 
cc          ln=nonn(jcl,3)
c          ldir(2)=ln3
c          xy(2,1)=xy(3,1)
c          xy(2,2)=xy(3,2)
c        endif
cc        CALL CSHAPE(EKS,NOD,XY,SHP,DSHP,RN,DETJ)

c       find the length of the element with nodes previus,current and next
        NGP=6
        call elemlen(xy,ngp,nod,eleg)
c       only for the current node 
        DO 6 I=1,NOD
        CALL CLRZRR(ZMA,49)
        CALL CLRZRR(ZV,7)
        CALL CLRZRR(ZW,7)
        ZSRR=dcmplx(0.0d0,0.0d0)
        ZSTT=dcmplx(0.0d0,0.0d0)
        ZSRT=dcmplx(0.0d0,0.0d0)
        if((nffs.eq.3).and.(icffs.eq.1).and.(ndtype(jcl).eq.3)) then
c          ndiak=1
c          jref=3*(ln2-1)
          jref=3*(ldir(i)-1)
          ZLOAD1(JREF+1)=dcmplx(0.0d0,0.0d0)
          ZLOAD1(JREF+2)=dcmplx(0.0d0,0.0d0)
          ZLOAD1(JREF+3)=dcmplx(0.0d0,0.0d0)
        endif

          if((ndtype(jcl).eq.2).and.(icffs.eq.2)) then
c            jref=3*(ln2-1)
            jref=3*(ldir(i)-1)
            ZLOADf1(JREF+1)=dcmplx(0.0d0,0.0d0)
            ZLOADf1(JREF+2)=dcmplx(0.0d0,0.0d0)
            ZLOADf1(JREF+3)=dcmplx(0.0d0,0.0d0)
          endif


        EKS=SPOINT(I)

        CALL CSHAPE(EKS,NOD,XY,SHP,DSHP,RN,DETJ)
c        rn2=RN(2)
c        rn1=RN(1)
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
c
c       x,y coordinates of tangential strain
        ZV(5)=ZV(5)/(0.5d0*ELEG)
        ZV(7)=ZV(7)/(0.5d0*ELEG)
c       #############
cc       another technic of calculating the Hoop stresses
cc       convert these strains in normal and tangential directions
c        zun=zv(5)*rn1+zv(7)*rn2
c        zus=-zv(5)*rn2+zv(7)*rn1
c        ztx=zv(2)
c        zty=zv(3)
c       #############
        ZMA(1,1)=dcmplx(-1.0d0,0.0d0)
c        ZMA(1,4)=RLA+2.0D0*RMU
        ZMA(1,4)=ZLA+2.0D0*ZMU
c        ZMA(1,7)=RLA
        ZMA(1,7)=ZLA
        ZMA(2,1)=RN(1)
        ZMA(2,3)=RN(2)
        ZMA(3,2)=RN(2)
        ZMA(3,3)=RN(1)
        ZMA(4,3)=dcmplx(-1.0d0,0.0d0)
        ZMA(4,5)=ZMU
c        ZMA(4,5)=RMU
        ZMA(4,6)=ZMU
c        ZMA(4,6)=RMU
        ZMA(5,4)=RQ(1)
        ZMA(5,6)=RQ(2)
        ZMA(6,2)=dcmplx(-1.0d0,0.0d0)
        ZMA(6,4)=ZLA
c        ZMA(6,4)=RLA
        ZMA(6,7)=ZLA+2.0D0*ZMU
c        ZMA(6,7)=RLA+2.0D0*RMU
        ZMA(7,5)=RQ(1)
        ZMA(7,7)=RQ(2)
        CALL ZMINV(ZMA,7,linv,minv)
        CALL ZVPROD(ZMA,ZV,7,ZW)
        Z11=ZW(1)
        Z22=ZW(2)
        Z12=ZW(3)
c        ZZ1(I,JCL)=Z11
c        ZZ2(I,JCL)=Z22
c        ZZ3(I,JCL)=Z12

c        ZZ1(JCL)=Z11
c        ZZ2(JCL)=Z22
c        ZZ3(JCL)=Z12
C              MODIFICATION   !!!!!!!!!!!!!!!!??????????
C              TRANSFORMATION IN RADIAL AND TANGENTIAL COORD.
        ZSRR=z11*RN(1)**2+z22*RN(2)**2+
     &             2*z12*RN(1)*RN(2)
        ZSTT=z11*RN(2)**2+z22*RN(1)**2-
     &             2*z12*RN(1)*RN(2)
        ZSRT=(z22-z11)*RN(1)*RN(2)+
     &             z12*(RN(1)**2-RN(2)**2)
  79    continue
c               ADD THE INCIDENT STRESSES (STRESSES WITHOUT THE EXISTENSE OF
C               THE CAVITY)
c       #############
cc       another technic of calculating the Hoop stresses
cc       convert these strains in normal and tangential directions
cc       transformation of the tractions from the global system to local
cc       traction normal to the surface of the element
c        jref=2*(ln2-1)
cc        ztx=zload(jref+1)
cc        zty=zload(jref+2)
c        ztn=ztx*rn1+zty*rn2
c        zts=-(ztx*rn2-zty*rn1)
cc       stresses in the normal(zsn) and tangential(zss) direction, 
cc       zsns is the shear stress
c        zsn=ztn
c        zsns=zts
c        zss=(rnu*zsn+2*rmu*zus)/(1.0d0-rnu)
c       #############

c       for the 1st loop
c        jref=3*(ln2-1)
        jref=3*(ldir(i)-1)

        if((nffs.eq.3).and.(icffs.eq.1)) then
c         the scattered stresses with possitive sign, because
c         the tractions in the 2nd loop will reverse (-), but the equevelant
c         stresses will be opposite (+)
c          ZLOAD1(JREF+1)=ZSRR+ZLOAD1(JREF+1)
c          ZLOAD1(JREF+2)=ZSTT+ZLOAD1(JREF+2)
c          ZLOAD1(JREF+3)=ZSRT+ZLOAD1(JREF+3) 

           ZLOAD1(JREF+1)=ZSRR+ZLOAD1(JREF+1)
           ZLOAD1(JREF+2)=ZSTT+ZLOAD1(JREF+2)
           ZLOAD1(JREF+3)=ZSRT+ZLOAD1(JREF+3) 

c          ZLOAD1(JREF+1)=ZSRR
c          ZLOAD1(JREF+2)=ZSTT 
c          ZLOAD1(JREF+3)=ZSRT

c          write(11,*) zz1p,zz2p,zz3p
          go to 99
        endif
c        if((nffs.eq.3).and.(icffs.eq.2)) then
c        ZZ1(JCL)=ZSRR+ZLOAD1(JREF+1)
c        ZZ2(JCL)=ZSTT+ZLOAD1(JREF+2)
c        ZZ3(JCL)=ZSRT+ZLOAD1(JREF+3)


c       case of free surface node
        if(ndtype(jcl).eq.2) then
c          ZZ1(JCL)=ZSRR+zloadf1(jref+1)
c          ZZ2(JCL)=ZSTT+zloadf1(jref+2)
c          ZZ3(JCL)=ZSRT+zloadf1(jref+3)
          ZZ1(I,JCL)=ZSRR+zloadf1(jref+1)
          ZZ2(I,JCL)=ZSTT+zloadf1(jref+2)
          ZZ3(I,JCL)=ZSRT+zloadf1(jref+3)
          go to 85
        endif
c        ZZ1(JCL)=ZSRR+ZLOAD1(JREF+1)+zloadf1(jref+1)
c        ZZ2(JCL)=ZSTT+ZLOAD1(JREF+2)+zloadf1(jref+2)
c        ZZ3(JCL)=ZSRT+ZLOAD1(JREF+3)+zloadf1(jref+3)
        ZZ1(I,JCL)=ZSRR+ZLOAD1(JREF+1)+zloadf1(jref+1)
        ZZ2(I,JCL)=ZSTT+ZLOAD1(JREF+2)+zloadf1(jref+2)
        ZZ3(I,JCL)=ZSRT+ZLOAD1(JREF+3)+zloadf1(jref+3)
c        ZZ1(JCL)=ZSRR-ZLOAD1(JREF+1)+zloadf1(jref+1)
c        ZZ2(JCL)=ZSTT-ZLOAD1(JREF+2)+zloadf1(jref+2)
c        ZZ3(JCL)=ZSRT-ZLOAD1(JREF+3)+zloadf1(jref+3)
c        ZZ1(JCL)=ZSRR+ZLOAD1(JREF+1)
c        ZZ2(JCL)=ZSTT+ZLOAD1(JREF+2)
c        ZZ3(JCL)=ZSRT+ZLOAD1(JREF+3)
c        ZZ1(JCL)=ZSRR+ZLOAD1(JREF+1)
c        ZZ2(JCL)=ZSTT+ZLOAD1(JREF+2)
c        ZZ3(JCL)=ZSRT+ZLOAD1(JREF+3)
c        ZZ1(JCL)=ZSRR+ZLOAD1(JREF+1)-zloadf1(jref+1)
c        ZZ2(JCL)=ZSTT+ZLOAD1(JREF+2)-zloadf1(jref+2)
c        ZZ3(JCL)=ZSRT+ZLOAD1(JREF+3)-zloadf1(jref+3)
  85    continue

c        go to 99
c        endif
  99    continue
  6     continue
  3     CONTINUE
        if((nffs.eq.3).and.(icffs.eq.1)) then
           rewind 11
           write(11,*) zload1
        endif
cc       add the hoop stresses of the 2 loops
c        if((nffs.eq.3).and.(icffs.eq.2)) then
c           rewind 11
c           read(11,*) zload1
c           do  71 JCL=1,nnode
c             ln2=nonn(jcl,2)
c             jref=3*(ln2-1)
c             ZZ1(JCL)=ZZ1(JCL)+ZLOAD1(JREF+1)
c             ZZ2(JCL)=ZZ2(JCL)+ZLOAD1(JREF+2)
c             ZZ3(JCL)=ZZ3(JCL)+ZLOAD1(JREF+3)
c  71       continue
c        endif


        RETURN
        END

c -------------------------------------------------------------------------
      subroutine elemlen(xy,ngp,nod,eleg)
c      
c     calculates the length of the element, by gaussian interpolation
c     nod: 3 for quadratic interpolation
c          2 for linear interpolation
c     xy : coordinates of the nodes of the element
c     ngp: number of points for gaussian interpolation
c     eleg: return value of the length of the element
c
      DIMENSION XY(3,2),SHP(3),DSHP(3),RN(2)
      COMMON/GAUSS/GPN(10,10),WTN(10,10),GPS(10,10),WTS(10,10)
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT REAL*8 (A-H,O-Y)
c
      ELEG=0.0d0
      DO 8 IX=1,NGP
        EK=GPN(NGP,IX)
        WT=WTN(NGP,IX)
        CALL CSHAPE(EK,NOD,XY,SHP,DSHP,RN,DETJ)
        ELEG=ELEG+DETJ*WT
 8    CONTINUE
      return
      end
c
c -------------------------------------------------------------------------
      subroutine viscoel(z,zmu,zla,znu,zc1,zc2)
c
c     calculates the complex 
c     z: Laplace transform parameter
c     zmu,zla: return value for the complex value of Lame's constants
      IMPLICIT COMPLEX*16 (Z)
      IMPLICIT REAL*8 (A-H,O-Y)
      COMMON/MATER/RLA,RMU,RNU,RHO,C1,C2,TIM,fvisc,nmode
c
      tol=1.0d-9
      zc1=c1
      zc2=c2
c     Kelvin solid
      if(dabs(fvisc).lt.tol) then
        zmu=rmu
        zla=rla
        go to 5
      endif
c     Kelvin solid
      if(nmode.eq.1) then 
c        zmu=rmu+rmu*fvisc*z
c        zla=rla-(2.0d0/3.0d0)*rmu*fvisc*z
        zmu=rmu+fvisc*z
        zla=rla-(2.0d0/3.0d0)*fvisc*z
cc       take the real and imaginery parts
c        zmu=rmu+fvisc*z*dcmplx(0.0d0,1.0d0)
c        zla=rla-(2.0d0/3.0d0)*fvisc*z*dcmplx(0.0d0,1.0d0)
      endif
c     Maxwell fluid
      if(nmode.eq.2) then
c        zmu=rmu*z/(z+1.0d0/fvisc)
c        zla=(rla*z+(1.0d0/fvisc)*(rla+(2.0d0/3.0d0)*rmu))/
c     &        (z+1.0d0/fvisc)
        zmu=rmu*z/(z+rmu/fvisc)
        zla=(rla*z+(rmu/fvisc)*(rla+(2.0d0/3.0d0)*rmu))/
     &        (z+rmu/fvisc)
c        zmu=rmu*z/(z+rmu/fvisc)*dcmplx(0.0d0,1.0d0)
c        zla=(rla*z+(rmu/fvisc)*(rla+(2.0d0/3.0d0)*rmu))/
c     &        (z+rmu/fvisc)*dcmplx(0.0d0,1.0d0)
      endif
  5   continue
      zc1=cdsqrt((zla+2.0d0*zmu)/rho)
      zc2=cdsqrt(zmu/rho)
      znu=zla/(2.0d0*(zla+zmu))
      return
      end
