PROJ = LAPLACE
PROJFILE = LAPLACE.MAK
DEBUG = 1

CC  = cl
CFLAGS_G  = /BATCH /FR$*.sbr
CFLAGS_D  = /Gi$(PROJ).mdt /Zi /Od
CFLAGS_R  = /Ot
FOR  = fl
FORFLAGS_G  = /AH /4Yb /FR$*.sbr
FORFLAGS_R  = /Ot
FORFLAGS_D  = /Od /Zi /Ge
MAPFILE_D  = NUL
MAPFILE_R  = NUL
LFLAGS_G  =  /STACK:40000  /BATCH
LFLAGS_D  =  /CO /FAR /PACKC
LFLAGS_R  = /EXE /FAR /PACKC
LINKER  = link
ILINK  = ilink
LRF  = echo > NUL
BRFLAGS  =  /o $(PROJ).bsc
PWBRMAKE  = pwbrmake
NMAKEBSC1  = set
NMAKEBSC2  = nmake
BROWSE  = 1

OBJS  = ..\BEMLAPL\LAPLBLOC.obj ..\BEMLAPL\LAPLB1.obj ..\BEMLAPL\LAPLB2.obj\
        ..\BEMLAPL\LAPLB3.obj
SBRS  = ..\BEMLAPL\LAPLBLOC.sbr ..\BEMLAPL\LAPLB1.sbr ..\BEMLAPL\LAPLB2.sbr\
        ..\BEMLAPL\LAPLB3.sbr

all: $(PROJ).exe

.SUFFIXES:
.SUFFIXES: .obj .sbr .for

..\BEMLAPL\LAPLBLOC.obj : ..\BEMLAPL\LAPLBLOC.FOR

..\BEMLAPL\LAPLBLOC.sbr : ..\BEMLAPL\LAPLBLOC.FOR

..\BEMLAPL\LAPLB1.obj : ..\BEMLAPL\LAPLB1.FOR

..\BEMLAPL\LAPLB1.sbr : ..\BEMLAPL\LAPLB1.FOR

..\BEMLAPL\LAPLB2.obj : ..\BEMLAPL\LAPLB2.FOR

..\BEMLAPL\LAPLB2.sbr : ..\BEMLAPL\LAPLB2.FOR

..\BEMLAPL\LAPLB3.obj : ..\BEMLAPL\LAPLB3.FOR

..\BEMLAPL\LAPLB3.sbr : ..\BEMLAPL\LAPLB3.FOR


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
