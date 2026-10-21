PROJ = BEMTIM
PROJFILE = BEMTIM.MAK
DEBUG = 1

ORIGIN  = PWB
ORIGIN_VER  = 2.0
CC  = cl
CFLAGS_G  = /W2 /BATCH
CFLAGS_D  = /f /Zi /Od
CFLAGS_R  = /f- /Ot /Oi /Ol /Oe /Og /Gs
CXX  = cl
CXXFLAGS_G  = /W2 /BATCH
CXXFLAGS_D  = /f /Zi /Od
CXXFLAGS_R  = /f- /Ot /Oi /Ol /Oe /Og /Gs
FOR  = fl
FORFLAGS_G  = /AH
FORFLAGS_R  = /Ot /Ol /Gs
FORFLAGS_D  = /Od /Zi /Ge /4Yb /Gt
ASM  = ml
AFLAGS_G  = /Cx /W2 /WX
AFLAGS_D  = /Zi
AFLAGS_R  = /nologo
MAPFILE_D  = NUL
MAPFILE_R  = NUL
LFLAGS_G  = /STACK:5000  /BATCH
LFLAGS_D  = /CO /FAR /PACKC
LFLAGS_R  = /EXE /FAR /PACKC
LINKER  = link
ILINK  = ilink
LRF  = echo > NUL
ILFLAGS  = /a /e
BC  = bc
BCFLAGS_R  = /Ot
BCFLAGS_D  = /D /Zi
BCFLAGS_G  = /O /G2 /Fpi /Lr
FILES  = ..\BEMTIM\BLOCKTIM.FOR ..\BEMTIM\BEMTIM1.FOR ..\BEMTIM\BEMTIM2.FOR
CVFLAGS  = /25 /F

OBJS  = ..\BEMTIM\TIMB1.obj ..\BEMTIM\TIMB2.obj ..\BEMTIM\TIMBLOC.obj

all: $(PROJ).exe

.SUFFIXES:
.SUFFIXES: .obj .for

..\BEMTIM\TIMB1.obj : ..\BEMTIM\TIMB1.FOR

..\BEMTIM\TIMB2.obj : ..\BEMTIM\TIMB2.FOR

..\BEMTIM\TIMBLOC.obj : ..\BEMTIM\TIMBLOC.FOR


$(PROJ).exe : $(OBJS)
!IF $(DEBUG)
        $(LRF) @<<$(PROJ).lrf
$(RT_OBJS: = +^
) $(OBJS: = +^
)
$@
$(MAPFILE_D)
$(LIBS: = +^
) +
$(LLIBS_G: = +^
) +
$(LLIBS_D: = +^
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
$(LIBS: = +^
) +
$(LLIBS_G: = +^
) +
$(LLIBS_R: = +^
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


run: $(PROJ).exe
        $(PROJ).exe $(RUNFLAGS)

debug: $(PROJ).exe
        CV $(CVFLAGS) $(PROJ).exe $(RUNFLAGS)
