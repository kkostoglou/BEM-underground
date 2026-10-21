PROJ = DYNAMICS
PROJFILE = DYNAMICS.MAK
DEBUG = 1

PWBRMAKE  = pwbrmake
NMAKEBSC1  = set
NMAKEBSC2  = nmake
CC  = cl
CFLAGS_G  = /W2 /BATCH
CFLAGS_D  = /Gi$(PROJ).mdt /Zi /Od
CFLAGS_R  = /Ot /Oi /Ol /Oe /Og /Gs
FOR  = fl
FORFLAGS_G  = /AH
FORFLAGS_R  = /Ot /Ol /Gs
FORFLAGS_D  = /Od /Zi /4Yb /Gt
MAPFILE_D  = NUL
MAPFILE_R  = NUL
LFLAGS_G  = /BATCH
LFLAGS_D  = /CO /FAR /PACKC
LFLAGS_R  = /EXE /FAR /PACKC
LINKER  = link
ILINK  = ilink
LRF  = echo > NUL
CVFLAGS  = /25 /X /F

OBJS  = MENU.obj
SBRS  = MENU.sbr

all: $(PROJ).exe

.SUFFIXES:
.SUFFIXES: .obj .sbr .for

MENU.obj : MENU.FOR EXEC.FI

MENU.sbr : MENU.FOR EXEC.FI


$(PROJ).bsc : $(SBRS)
        $(PWBRMAKE) @<<
$(BRFLAGS) $(SBRS)
<<

$(PROJ).exe : $(OBJS)
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
$(DEF_FILE) $(LFLAGS_G) $(LFLAGS_D);
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
$(DEF_FILE) $(LFLAGS_G) $(LFLAGS_R);
<<
!ENDIF
        $(LINKER) @$(PROJ).lrf


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
        $(PROJ).exe $(RUNFLAGS)

debug: $(PROJ).exe
        CV $(CVFLAGS) $(PROJ).exe $(RUNFLAGS)
