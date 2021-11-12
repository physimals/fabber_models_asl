include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_asl
LIBS = -lfsl-fabberexec -lfsl-fabbercore -lfsl-newimage \
       -lfsl-miscmaths -lfsl-NewNifti -lfsl-utils \
       -lfsl-cprob -lfsl-znz -ldl
XFILES = fabber_asl
SOFILES = libfsl-fabber_models_asl.so

# Forward models
OBJS =  fwdmodel_asl_multiphase.o fwdmodel_asl_grase.o asl_models.o fwdmodel_asl_rest.o \
        fwdmodel_asl_quasar.o fwdmodel_asl_satrecov.o fwdmodel_asl_satrecovdualfa.o fwdmodel_asl_turboquasar.o \
	fwdmodel_asl_2compartment.o fwdmodel_asl_multite.o fwdmodel_asl_velocityselective.o

# For debugging:
#OPTFLAGS = -ggdb

# Pass Git revision details
GIT_SHA1:=$(shell git describe --dirty)
GIT_DATE:=$(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

#
# Build
#

all: ${XFILES} ${SOFILES}

# models in a library
libfsl-fabber_models_asl.so : ${OBJS}
	${CXX} ${CXXFLAGS} -shared -o $@ $^ ${LDFLAGS}

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber_asl : fabber_client.o | libfsl-fabber_models_asl.so
	${CXX} ${CXXFLAGS} -o $@ $< -lfsl-fabber_models_asl ${LDFLAGS}
