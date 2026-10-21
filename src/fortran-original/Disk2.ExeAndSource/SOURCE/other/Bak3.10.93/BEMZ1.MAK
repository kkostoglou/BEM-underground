PROJ = BEMZ1
PROJFILE = BEMZ1.MAK
DEBUG = 1

CC  = cl
CFLAGS_G  = /BATCH
CFLAGS_D  = /Gi$(PROJ).mdt /Zi /Od
CFLAGS_R  = /Ot
FOR  = fl
FORFLAGS_G  = /AH /4Yb
FORFLAGS_R  = /Ot
FORFLAGS_D  = /Od /Zi /Ge
MAPFILE_D  = NUL
MAPFILE_R  = NUL
LFLAGS_G  =  /STACK:30000  /BATCH
LFLAGS_D  =  /CO /FAR /PACKC
LFLAGS_R  = /EXE /FAR /PACKC
LINKER  = link
ILINK  = ilink
LRF  = echo > NUL
BRFLAGS  =  /o $(PROJ).bsc
PWBRMAKE  = pwbrmake
NMAKEBSC1  = set
NMAKEBSC2  = nmake
BROWSE  = 0

OBJS  = ..\BEMZ\BEMZ1.obj ..\BEMZ\BEMZ2.obj ..\BEMZ\BEMZ3.obj\
        ..\BEMZ\BEMZBLOC.obj ..\BEMZ\BEMZ4.obj
SBRS  = ..\BEMZ\BEMZ1.sbr ..\BEMZ\BEMZ2.sbr ..\BEMZ\BEMZ3.sbr\
        ..\BEMZ\BEMZBLOC.sbr ..\BEMZ\BEMZ4.sbr

all: $(PROJ).exe

.SUFFIXES:
.SUFFIXES: .obj .sbr .for

..\BEMZ\BEMZ1.obj : ..\BEMZ\BEMZ1.FOR

..\BEMZ\BEMZ1.sbr : ..\BEMZ\BEMZ1.FOR

..\BEMZ\BEMZ2.obj : ..\BEMZ\BEMZ2.FOR

..\BEMZ\BEMZ2.sbr : ..\BEMZ\BEMZ2.FOR

..\BEMZ\BEMZ3.obj : ..\BEMZ\BEMZ3.FOR

..\BEMZ\BEMZ3.sbr : ..\BEMZ\BEMZ3.FOR

..\BEMZ\BEMZBLOC.obj : ..\BEMZ\BEMZBLOC.FOR

..\BEMZ\BEMZBLOC.sbr : ..\BEMZ\BEMZBLOC.FOR

..\BEMZ\BEMZ4.obj : ..\BEMZ\BEMZ4.FOR

..\BEMZ\BEMZ4.sbr : ..\BEMZ\BEMZ4.FOR


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
        $(PROJ).exe $(RUNFLAGS)

debug: $(PROJ).exe
        CV $(CVFLAGS) $(PROJ).exe $(RUNFLAGS)
