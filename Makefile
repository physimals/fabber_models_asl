include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_asl

USRINCFLAGS = -I${INC_NEWMAT} -I${INC_PROB} -I${INC_BOOST}
USRLDFLAGS = -L${LIB_NEWMAT} -L${LIB_PROB}

LIBS = -lutils -lnewimage -lmiscmaths -lprob -lnewmat -lfslio -lniftiio -lznz -lz -lfabbercore 

XFILES = fabber

# Forward models
OBJS =  fwdmodel_asl_multiphase.o fwdmodel_asl_grase.o asl_models.o fwdmodel_asl_rest.o

# For debugging:
OPTFLAGS = -ggdb
#OPTFLAGS =

#
# Build
#

all:	${XFILES} libfabbermodels_asl.a

# models in a library
libfabbermodels_asl.a : ${OBJS}
	${AR} -r $@ ${OBJS}

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber : fabber_client.o ${OBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o $@ $< ${OBJS} -lfabbercore ${LIBS}

# DO NOT DELETE
