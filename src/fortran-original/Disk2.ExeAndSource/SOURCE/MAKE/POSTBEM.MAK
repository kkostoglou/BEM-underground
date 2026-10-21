PROJ = POSTBEM
PROJFILE = POSTBEM.MAK
DEBUG = 1

PWBRMAKE  = pwbrmake
NMAKEBSC1  = set
NMAKEBSC2  = nmake
CC  = cl
CFLAGS_G  = /BATCH
CFLAGS_D  = /Gi$(PROJ).mdt /Zi /Od
CFLAGS_R  = /Ot
FOR  = fl
FORFLAGS_G  = /AH
FORFLAGS_R  = /Ot
FORFLAGS_D  = /Od /Zi /Ge /4Yb /Gt
MAPFILE_D  = NUL
MAPFILE_R  = NUL
LFLAGS_G  = /STACK:5000  /BATCH
LFLAGS_D  = /CO /FAR /PACKC
LFLAGS_R  = /EXE /FAR /PACKC
LINKER  = link
ILINK  = ilink
LRF  = echo > NUL
LLIBS_G  = GRAPHICS
CVFLAGS  = /25 /X /S

OBJS  = ..\POSTBEM\POSTBLOC.obj ..\POSTBEM\POSTB1.obj ..\POSTBEM\POSTB2.obj\
        ..\POSTBEM\POSTB3.obj
SBRS  = ..\POSTBEM\POSTBLOC.sbr ..\POSTBEM\POSTB1.sbr ..\POSTBEM\POSTB2.sbr\
        ..\POSTBEM\POSTB3.sbr

all: $(PROJ).exe

.SUFFIXES:
.SUFFIXES: .obj .sbr .for

..\POSTBEM\POSTBLOC.obj : ..\POSTBEM\POSTBLOC.FOR

..\POSTBEM\POSTBLOC.sbr : ..\POSTBEM\POSTBLOC.FOR

..\POSTBEM\POSTB1.obj : ..\POSTBEM\POSTB1.FOR C:\fortran\INCLUDE\FGRAPH.FI\
        C:\fortran\INCLUDE\FGRAPH.FD

..\POSTBEM\POSTB1.sbr : ..\POSTBEM\POSTB1.FOR C:\fortran\INCLUDE\FGRAPH.FI\
        C:\fortran\INCLUDE\FGRAPH.FD

..\POSTBEM\POSTB2.obj : ..\POSTBEM\POSTB2.FOR C:\fortran\INCLUDE\FGRAPH.FI\
        C:\fortran\INCLUDE\FGRAPH.FD

..\POSTBEM\POSTB2.sbr : ..\POSTBEM\POSTB2.FOR C:\fortran\INCLUDE\FGRAPH.FI\
        C:\fortran\INCLUDE\FGRAPH.FD

..\POSTBEM\POSTB3.obj : ..\POSTBEM\POSTB3.FOR C:\fortran\INCLUDE\FGRAPH.FI\
        C:\fortran\INCLUDE\FGRAPH.FD

..\POSTBEM\POSTB3.sbr : ..\POSTBEM\POSTB3.FOR C:\fortran\INCLUDE\FGRAPH.FI\
        C:\fortran\INCLUDE\FGRAPH.FD


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
