PROJ = LAPLACEW
PROJFILE = LAPLACEW.MAK
DEBUG = 1

CC  = cl
CFLAGS_G  = /W2 /Gw /Zp /BATCH /FR$*.sbr
CFLAGS_D  = /Zi /Gi$(PROJ).mdt /Od /Gs /FPa
CFLAGS_R  = /Oe /Og /Os /Gs /FPa
FOR  = fl
FORFLAGS_G  = /MW1 /G2 /FR$*.sbr
FORFLAGS_R  = /Ot /Ol /Gs
FORFLAGS_D  = /Zi /Od /4Yb
MAPFILE_D  = NUL
MAPFILE_R  = NUL
LFLAGS_G  = /STACK:10000  /BATCH
LFLAGS_D  = /CO /FAR /PACKC
LFLAGS_R  = /NOF
LINKER  = link
ILINK  = ilink
LRF  = echo > NUL
RC  = rc
LLIBS_R  = /NOD:LLIBFOR7 LLIBFEW
LLIBS_D  = /NOD:LLIBFOR7 LLIBFEW
BRFLAGS  =  /o $(PROJ).bsc
PWBRMAKE  = pwbrmake
NMAKEBSC1  = set
NMAKEBSC2  = nmake
BROWSE  = 1

DEF_FILE  = FL.DEF
OBJS  = ..\BEMLAPL\LAPLBLOC.obj ..\BEMLAPL\LAPLB1.obj ..\BEMLAPL\LAPLB2.obj\
        ..\BEMLAPL\LAPLB3.obj
LIBS_EXT  = ..\..\..\FORTRAN\LIB\LLIBFEW.LIB
LIBS  = $(LIBS_EXT)

all: $(PROJ).exe

.SUFFIXES:
.SUFFIXES: .obj .for

..\BEMLAPL\LAPLBLOC.obj : ..\BEMLAPL\LAPLBLOC.FOR

..\BEMLAPL\LAPLB1.obj : ..\BEMLAPL\LAPLB1.FOR

..\BEMLAPL\LAPLB2.obj : ..\BEMLAPL\LAPLB2.FOR

..\BEMLAPL\LAPLB3.obj : ..\BEMLAPL\LAPLB3.FOR


$(PROJ).exe : $(DEF_FILE) $(OBJS) $(LIBS)
!IF $(DEBUG)
        $(LRF) @<<$(PROJ).lrf
$(RT_OBJS: = +^
) $(OBJS: = +^
)
$@
$(MAPFILE_D)
$(LLIBS_G: = +^
) +
$(LLIBS_D: = +^
) +
$(LIBS: = +^
)
$(FOR_DEFS) $(DEF_FILE) $(LFLAGS_G) $(LFLAGS_D);
<<
!ELSE
        $(LRF) @<<$(PROJ).lrf
$(RT_OBJS: = +^
) $(OBJS: = +^
)
$@
$(MAPFILE_R)
$(LLIBS_G: = +^
) +
$(LLIBS_R: = +^
) +
$(LIBS: = +^
)
$(FOR_DEFS) $(DEF_FILE) $(LFLAGS_G) $(LFLAGS_R);
<<
!ENDIF
        $(LINKER) @$(PROJ).lrf
        $(NMAKEBSC1) MAKEFLAGS=
        $(NMAKEBSC2) $(NMFLAGS) -f $(PROJFILE) $(PROJ).bsc

$(PROJ).bsc : 


.for.obj :
!IF $(DEBUG)
        $(FOR) /c $(FORFLAGS_G) $(FORFLAGS_D) /Fo$@ $<
!ELSE
        $(FOR) /c $(FORFLAGS_G) $(FORFLAGS_R) /Fo$@ $<
!ENDIF


run: $(PROJ).exe
        WIN $(PROJ).exe $(RUNFLAGS)

debug: $(PROJ).exe
        WIN CVW $(CVFLAGS) $(PROJ).exe $(RUNFLAGS)
