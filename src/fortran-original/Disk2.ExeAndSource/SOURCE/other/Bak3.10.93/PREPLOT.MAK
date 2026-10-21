PROJ = PREPLOT
PROJFILE = PREPLOT.MAK
DEBUG = 1

PWBRMAKE  = pwbrmake
NMAKEBSC1  = set
NMAKEBSC2  = nmake
BRFLAGS  = /o $(PROJ).bsc
BROWSE  = 1
CC  = cl
CFLAGS_G  = /W2 /BATCH /FR$*.sbr
CFLAGS_D  = /Gi$(PROJ).mdt /Zi /Od
CFLAGS_R  = /Ot /Oi /Ol /Oe /Og /Gs
FOR  = fl
FORFLAGS_G  = /AH /FR$*.sbr
FORFLAGS_R  = /Ot /Ol /Gs
FORFLAGS_D  = /Od /Zi /Ge /4Yb
MAPFILE_D  = NUL
MAPFILE_R  = NUL
LFLAGS_G  = /STACK:20000  /BATCH
LFLAGS_D  = /CO /FAR /PACKC
LFLAGS_R  = /EXE /FAR /PACKC
LINKER  = link
ILINK  = ilink
LRF  = echo > NUL
LLIBS_G  = GRAPHICS
CVFLAGS  = /25 /S

OBJS  = ..\PREPLOT\PREPLOT1.obj ..\PREPLOT\PREPLOT2.obj\
        ..\PREPLOT\PREPLOT3.obj

all: $(PROJ).exe

.SUFFIXES:
.SUFFIXES: .obj .for

..\PREPLOT\PREPLOT1.obj : ..\PREPLOT\PREPLOT1.FOR C:\fortran\INCLUDE\FGRAPH.FI\
        C:\fortran\INCLUDE\FGRAPH.FD

..\PREPLOT\PREPLOT2.obj : ..\PREPLOT\PREPLOT2.FOR C:\fortran\INCLUDE\FGRAPH.FD

..\PREPLOT\PREPLOT3.obj : ..\PREPLOT\PREPLOT3.FOR C:\fortran\INCLUDE\FGRAPH.FD


$(PROJ).bsc : 

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
        $(NMAKEBSC1) MAKEFLAGS=
        $(NMAKEBSC2) $(NMFLAGS) -f $(PROJFILE) $(PROJ).bsc


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
