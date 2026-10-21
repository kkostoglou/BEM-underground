PROJ = BEMZ
PROJFILE = BEMZ.MAK
DEBUG = 1

PWBRMAKE  = pwbrmake
NMAKEBSC1  = set
NMAKEBSC2  = nmake
CC  = cl
CFLAGS_G  = /W2 /Gw /Zp /BATCH
CFLAGS_D  = /Zi /Gi$(PROJ).mdt /Od /Gs /FPa
CFLAGS_R  = /Oe /Og /Os /Gs /FPa
FOR  = fl
FORFLAGS_G  = /AH /MW1 /G2
FORFLAGS_R  = /Ot /Ol /Gs
FORFLAGS_D  = /Od /Zi /4Yb
MAPFILE_D  = NUL
MAPFILE_R  = NUL
LFLAGS_G  =  /STACK:40000  /BATCH
LFLAGS_D  =  /CO /FAR /PACKC
LFLAGS_R  = /NOF
LINKER  = link
ILINK  = ilink
LRF  = echo > NUL
RC  = rc
LLIBS_R  =  /NOD:LLIBFOR7 LLIBFEW
LLIBS_D  =  /NOD:LLIBFOR7 LLIBFEW

DEF_FILE  = FL.DEF
OBJS  = ..\BEMZ\BEMZ1.obj ..\BEMZ\BEMZ2.obj ..\BEMZ\BEMZ3.obj\
        ..\BEMZ\BEMZBLOC.obj ..\BEMZ\BEMZ4.obj

all: $(PROJ).exe

.SUFFIXES:
.SUFFIXES: .obj .for

..\BEMZ\BEMZ1.obj : ..\BEMZ\BEMZ1.FOR

..\BEMZ\BEMZ2.obj : ..\BEMZ\BEMZ2.FOR

..\BEMZ\BEMZ3.obj : ..\BEMZ\BEMZ3.FOR

..\BEMZ\BEMZBLOC.obj : ..\BEMZ\BEMZBLOC.FOR

..\BEMZ\BEMZ4.obj : ..\BEMZ\BEMZ4.FOR


$(PROJ).bsc : 

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
