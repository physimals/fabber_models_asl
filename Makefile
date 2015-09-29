include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_models_asl

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB}

LIBS = -lutils -lnewimage -lmiscmaths -lprob -lnewmat -lfslio -lniftiio -lznz -lz -lfabbercore 

# Forward models
OBJS =  fwdmodel_asl_multiphase.o fwdmodel_asl_grase.o asl_models.o fwdmodel_asl_rest.o

# For debugging:
OPTFLAGS = -ggdb
#OPTFLAGS =

#
# Build
#

all:	libfabbermodels_asl.a

# models only in a library
libfabbermodels_asl.a : ${OBJS}
	${AR} -r $@ ${OBJS}

# DO NOT DELETE
