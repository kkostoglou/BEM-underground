PROJ = BEMLAPL
PROJFILE = BEMLAPL.MAK
DEBUG = 1

PWBRMAKE  = pwbrmake
NMAKEBSC1  = set
NMAKEBSC2  = nmake
CC  = cl
CFLAGS_G  = /BATCH
CFLAGS_D  = /Gi$(PROJ).mdt /Zi /Od
CFLAGS_R  = /Ot
FOR  = fl
FORFLAGS_R  = /Ot
FORFLAGS_D  = /Od /Zi /Ge /4Yb
MAPFILE_D  = NUL
MAPFILE_R  = NUL
LFLAGS_G  =  /STACK:30000  /BATCH
LFLAGS_D  = /CO /FAR /PACKC
LFLAGS_R  = /EXE /FAR /PACKC
LINKER  = link
ILINK  = ilink
LRF  = echo > NUL

OBJS  = ..\BEMLAPL\LAPLBLOC.obj ..\BEMLAPL\LAPLB1.obj ..\BEMLAPL\LAPLB2.obj\
        ..\BEMLAPL\LAPLB3.obj

all: $(PROJ).exe

.SUFFIXES:
.SUFFIXES: .obj .for

..\BEMLAPL\LAPLBLOC.obj : ..\BEMLAPL\LAPLBLOC.FOR

..\BEMLAPL\LAPLB1.obj : ..\BEMLAPL\LAPLB1.FOR

..\BEMLAPL\LAPLB2.obj : ..\BEMLAPL\LAPLB2.FOR

..\BEMLAPL\LAPLB3.obj : ..\BEMLAPL\LAPLB3.FOR


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
