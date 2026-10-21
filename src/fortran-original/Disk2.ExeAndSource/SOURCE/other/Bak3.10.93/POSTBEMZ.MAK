PROJ = POSTBEMZ
PROJFILE = POSTBEMZ.MAK
DEBUG = 1

PWBRMAKE  = pwbrmake
NMAKEBSC1  = set
NMAKEBSC2  = nmake
CC  = cl
CFLAGS_G  = /W2 /BATCH
CFLAGS_D  = /Gi$(PROJ).mdt /Zi /Od
CFLAGS_R  = /Ot /Oi /Ol /Oe /Og /Gs
FOR  = fl
FORFLAGS_R  = /Ot /Ol /Gs
FORFLAGS_D  = /Od /Zi /4Yb
MAPFILE_D  = NUL
MAPFILE_R  = NUL
LFLAGS_G  = /NOI  /BATCH
LFLAGS_D  = /CO /FAR /PACKC
LFLAGS_R  = /EXE /FAR /PACKC
LINKER  = link
ILINK  = ilink
LRF  = echo > NUL
LLIBS_G  = GRAPHICS
CVFLAGS  = /25 /S

OBJS  = ..\POSTBEMZ\POSTBLOZ.obj ..\POSTBEMZ\POSTBZ1.obj\
        ..\POSTBEMZ\POSTBZ2.obj ..\POSTBEMZ\POSTBZ3.obj

all: $(PROJ).exe

.SUFFIXES:
.SUFFIXES: .obj .for

..\POSTBEMZ\POSTBLOZ.obj : ..\POSTBEMZ\POSTBLOZ.FOR

..\POSTBEMZ\POSTBZ1.obj : ..\POSTBEMZ\POSTBZ1.FOR C:\fortran\INCLUDE\FGRAPH.FI\
        C:\fortran\INCLUDE\FGRAPH.FD

..\POSTBEMZ\POSTBZ2.obj : ..\POSTBEMZ\POSTBZ2.FOR C:\fortran\INCLUDE\FGRAPH.FI\
        C:\fortran\INCLUDE\FGRAPH.FD

..\POSTBEMZ\POSTBZ3.obj : ..\POSTBEMZ\POSTBZ3.FOR C:\fortran\INCLUDE\FGRAPH.FI\
        C:\fortran\INCLUDE\FGRAPH.FD


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
