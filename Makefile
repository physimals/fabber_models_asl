include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_asl
XFILES   = fabber_asl
SOFILES  = libfsl-fabber_models_asl.so
AFILES   = libfabber_models_asl.a

# The FSL build system changed
# substantially in FSL 6.0.6
# FSL >= 6.0.6
ifeq (${FSL_GE_606}, true)
  LIBS = -lfsl-fabberexec -lfsl-fabbercore -lfsl-newimage \
         -lfsl-miscmaths -lfsl-NewNifti -lfsl-utils \
         -lfsl-cprob -lfsl-znz -ldl
# FSL <= 6.0.5
else
  ifeq ($(shell uname -s), Linux)
	MATLIB := -lopenblas
  endif

  USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_CPROB} \
                -I${INC_BOOST} -I.. -I${FSLDIR}/extras/include/armawrap
  USRLDFLAGS  = -L${LIB_NEWMAT} -L${LIB_PROB} -L../fabber_core  \
                -lnewimage -lmiscmaths -lutils -lprob ${MATLIB} \
                -lNewNifti -lznz -lz -ldl
endif


# Forward models
OBJS = fwdmodel_asl_multiphase.o fwdmodel_asl_grase.o asl_models.o fwdmodel_asl_rest.o \
       fwdmodel_asl_quasar.o fwdmodel_asl_satrecov.o fwdmodel_asl_satrecovdualfa.o fwdmodel_asl_turboquasar.o \
       fwdmodel_asl_2compartment.o fwdmodel_asl_multite.o fwdmodel_asl_velocityselective.o

# For debugging:
#OPTFLAGS = -ggdb

# Pass Git revision details
GIT_SHA1 := $(shell git describe --dirty)
GIT_DATE := $(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

#
# Build
#

# FSL >=606 uses dynamic linking
ifeq (${FSL_GE_606}, true)
all: ${XFILES} ${SOFILES}

# models in a library
libfsl-fabber_models_asl.so : ${OBJS}
	${CXX} ${CXXFLAGS} -shared -o $@ $^ ${LDFLAGS}

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber_asl : fabber_client.o | libfsl-fabber_models_asl.so
	${CXX} ${CXXFLAGS} -o $@ $< -lfsl-fabber_models_asl ${LDFLAGS}

# FSL <=605 uses static linking
else
all: ${XFILES} ${AFILES}

libfabber_models_asl.a : ${OBJS}
	${AR} -r $@ ${OBJS}

fabber_asl : fabber_client.o ${OBJS}
	${CXX} ${CXXFLAGS} -o $@ $^ -lfabbercore -lfabberexec ${LDFLAGS}
endif
