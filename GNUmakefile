ifndef NEUT_ROOT
      NEUT_ROOT = ../../
endif
ifndef GIBUU_ROOT
      GIBUU_ROOT = /home/tfeusels/T2K/GiBUU/release1.6
endif

include $(NEUT_ROOT)/src/neutsmpl/config_neut.gmk


NEUT_SRC = $(NEUT_ROOT)/src
#from compiling GiBUU all .o and .mods are there!!
GIBUU_SRC = $(GIBUU_ROOT)/objects

#FC = g77
FC = gfortran

ROOTINC = -I`root-config --incdir`
NEUTINC = -I../neutcore -I../nuccorspl -I../skmcsvc -I../nuceff
CLASSDIR = $(PWD)/../neutclass
CLASSINC= -I${CLASSDIR}
# NEW
GIBUUINC = -I${GIBUU_SRC}
INCDIRS = ${ROOTINC} ${NEUTINC} ${CLASSINC} ${GIBUUINC}

FCOPTIONS = -g

CXXFLAGS  = `root-config --cflags`
COPTFLAGS = -g -Wno-write-strings -fpermissive


COREDIR =	../neutcore/Linux_pc

LIBDIRS =	-L${COREDIR} -L../nuccorspl/Linux_pc \
		-L../nuceff/Linux_pc -L../partnuck/Linux_pc \
		-L../skmcsvc/Linux_pc -L../specfunc/Linux_pc -L../specfunc

NEUTCOREVER	= 5.3.2
INEUTCOREVER= 532

NUCEFFVER   = 5.3.1
INUCEFFVER  = 531

NUCCORVER   = 5.3.1
INUCCORVER  = 531

PARTNUCKVER = 5.0.5
IPARTNUCKVER= 505

SKMCSVCVER  = 5.0.5
ISKMCSVCVER = 505

SPECFUNCVER = 5.3.2
ISPECFUNCVER = 532

FORTRANDEFINES += -DNECORE_VERSION=$(INEUTCOREVER)
FORTRANDEFINES += -DNENUCE_VERSION=$(INUCEFFVER)
FORTRANDEFINES += -DNENUCC_VERSION=$(INUCCORVER)

MCLIB = ${LIBDIRS}   -lneutcore_${NEUTCOREVER} -lnuceff_${NUCEFFVER} \
		     -lneutcore_${NEUTCOREVER} -lnuccorrspl_${NUCCORVER} \
                     -lpartnuck_${PARTNUCKVER} -lskmcsvc_${SKMCSVCVER} -lspecfunc_${SPECFUNCVER}



CLIBS = ${CERN}/${CERN_LEVEL}/lib/libjetset74.a \
		${CERN}/${CERN_LEVEL}/lib/libpdflib804.a \
		`cernlib jetset74 photos202 mathlib packlib kernlib`

ROOTLIBS  = `root-config --libs` -lCore -lCint -lEG -lPhysics -lRIO -lNet \
        -lTree -lGeom -lGraf -lGraf3d -lHist -lMatrix -lMinuit -lPostscript \
	-lTree -lMathCore -lGpad -lGui -lGX11 -lRint -lThread

LIBS      = $(ROOTLIBS) ${MCLIB} ${CLIBS} -lstdc++

LDOPTFLAGS= -g

FINCDIRS  = ${NEUTINC}

SOBJS = ${CLASSDIR}/neutctrl.so ${CLASSDIR}/neutpart.so \
		${CLASSDIR}/neutvect.so ${CLASSDIR}/neutvtx.so

ROBJS = ${CLASSDIR}/neutfill.o \
		${CLASSDIR}/event_ratefortwrapper.o \
		${CLASSDIR}/NeutRootHandlers.o 

POBJS = ${COREDIR}/structm.o ${COREDIR}/pdfset.o ${COREDIR}/grv94di.o \
		${COREDIR}/grv98_lo.o
MOBJS = ${NEUT_SRC}/neutsmpl/nevecgen.o ${NEUT_SRC}/neutsmpl/grndmq.o ${NEUT_SRC}/neutsmpl/rmarin_dum.o

#OBJS = neutev.o neutxs.o TNuTrajectory.o TNeutOutput.o TNuFlux.o
OBJS = neutev.o neutxs.o gibuuxs.o TNuTrajectory.o TNeutOutput.o TNuFlux.o


# HERE I NEED TO PUT ALL *.o to make sure the linking works fine
GIBUU_OBJS := $(wildcard ${GIBUU_SRC}/*.o )
GIBUU_LIBS = -lbz2 #-lPDF

.SUFFIXES:	.so

GENROOTSO = env COPTFLAGS="${COPTFLAGS}" INCDIRS="${INCDIRS}" \
			${NEUT_SRC}/neutsmpl/bin/gen_root_so.sh

.cc.o:
	$(CXX) -c $(COPTFLAGS) $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cxx.o:
	$(CXX) -c $(COPTFLAGS) $(CXXFLAGS) $(INCDIRS) -o $@ $<

.cc.so:
	$(GENROOTSO) $(basename $<)

.F.o:
	$(FC) $(DEFS) -c $(FCOPTIONS) $(ALLDEFINES) $(FINCDIRS) -o $@ $<

## GIBUU (.f90.o does NOT work!)
%.o: %.f90
	$(FC) $(DEFS) -c $(FCOPTIONS) $(ALLDEFINES) $(GIBUUINC) -o $@ $<

all: event_rate genev

event_rate: event_rate.o $(OBJS) 
	$(FC) -O3 -o $@ event_rate.o -Xlinker -R`pwd`  $(OBJS) ${MOBJS} $(GIBUU_OBJS) $(LIBS) $(GIBUU_LIBS)

genev: genev.o $(OBJS)
	$(FC) -O3 -o $@ genev.o -Xlinker -R`pwd`  $(OBJS) ${MOBJS} $(GIBUU_OBJS) $(LIBS) $(GIBUU_LIBS)

dumptotpau_nd280: dumptotpau_nd280.o $(OBJS)
	$(FC) -O3 -o $@ dumptotpau_nd280.o -Xlinker -R`pwd` ${OBJS} ${MOBJS} $(LIBS)

clean:
	$(RM) -f *.o *~ ${OBJS} event_rate genev dumptotpau_nd280







#cd objects for both
#$(MAKE) $(noPrintDirectory) compileOBJ
#$(MAKE) $(noPrintDirectory) buildLIBS
