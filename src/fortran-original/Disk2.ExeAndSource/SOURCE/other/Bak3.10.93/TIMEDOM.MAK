PROJ = TIMEDOM
PROJFILE = TIMEDOM.MAK
DEBUG = 1

CC  = cl
CFLAGS_G  = /W2 /Gw /Zp /BATCH /FR$*.sbr
CFLAGS_D  = /Zi /Gi$(PROJ).mdt /Od /Gs /FPa
CFLAGS_R  = /Oe /Og /Os /Gs /FPa
FOR  = fl
FORFLAGS_G  = /AH /MW1 /G2 /FR$*.sbr
FORFLAGS_R  = /Ot /Ol /Gs
FORFLAGS_D  = /Od /Zi /Ge /4Yb
MAPFILE_D  = NUL
MAPFILE_R  = NUL
LFLAGS_G  =  /NOE /STACK:5000  /BATCH
LFLAGS_D  = /CO /FAR /PACKC
LFLAGS_R  = /NOF
LINKER  = link
ILINK  = ilink
LRF  = echo > NUL
RC  = rc
LLIBS_R  =  /NOD:LLIBFOR7 LLIBFEW
LLIBS_D  =  /NOD:LLIBFOR7 LLIBFEW
BRFLAGS  =  /o $(PROJ).bsc
PWBRMAKE  = pwbrmake
NMAKEBSC1  = set
NMAKEBSC2  = nmake
BROWSE  = 1

DEF_FILE  = FL.DEF
OBJS  = ..\BEMTIM\TIMB1.obj ..\BEMTIM\TIMB2.obj ..\BEMTIM\TIMBLOC.obj\
        ..\BEMTIM\TIMB3.obj
SBRS  = ..\BEMTIM\TIMB1.sbr ..\BEMTIM\TIMB2.sbr ..\BEMTIM\TIMBLOC.sbr\
        ..\BEMTIM\TIMB3.sbr

all: $(PROJ).exe

.SUFFIXES:
.SUFFIXES: .obj .sbr .for

..\BEMTIM\TIMB1.obj : ..\BEMTIM\TIMB1.FOR

..\BEMTIM\TIMB1.sbr : ..\BEMTIM\TIMB1.FOR

..\BEMTIM\TIMB2.obj : ..\BEMTIM\TIMB2.FOR

..\BEMTIM\TIMB2.sbr : ..\BEMTIM\TIMB2.FOR

..\BEMTIM\TIMBLOC.obj : ..\BEMTIM\TIMBLOC.FOR

..\BEMTIM\TIMBLOC.sbr : ..\BEMTIM\TIMBLOC.FOR

..\BEMTIM\TIMB3.obj : ..\BEMTIM\TIMB3.FOR

..\BEMTIM\TIMB3.sbr : ..\BEMTIM\TIMB3.FOR


$(PROJ).exe : $(DEF_FILE) $(OBJS)
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

$(PROJ).bsc : $(SBRS)
        $(PWBRMAKE) @<<
$(BRFLAGS) $(SBRS)
<<


.for.obj :
!IF $(DEBUG)
        $(FOR) /c $(FORFLAGS_G) $(FORFLAGS_D) /Fo$@ $<
!ELSE
        $(FOR) /c $(FORFLAGS_G) $(FORFLAGS_R) /Fo$@ $<
!ENDIF

.for.sbr :
!IF $(DEBUG)
        $(FOR) /Zs $(FORFLAGS_G) $(FORFLAGS_D) /FR$@ $<
!ELSE
        $(FOR) /Zs $(FORFLAGS_G) $(FORFLAGS_R) /FR$@ $<
!ENDIF


run: $(PROJ).exe
        WIN $(PROJ).exe $(RUNFLAGS)

debug: $(PROJ).exe
        WIN CVW $(CVFLAGS) $(PROJ).exe $(RUNFLAGS)
