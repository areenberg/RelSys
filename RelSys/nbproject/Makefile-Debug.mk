#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux
CND_DLIB_EXT=so
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/Combinatorial.o \
	${OBJECTDIR}/Customer.o \
	${OBJECTDIR}/CustomerData.o \
	${OBJECTDIR}/EntireSystem.o \
	${OBJECTDIR}/Experiments.o \
	${OBJECTDIR}/HeuristicQueue.o \
	${OBJECTDIR}/HyperQueue.o \
	${OBJECTDIR}/LinSolver.o \
	${OBJECTDIR}/PhaseFitter.o \
	${OBJECTDIR}/QueueData.o \
	${OBJECTDIR}/QueuePerformance.o \
	${OBJECTDIR}/RelocEvaluation.o \
	${OBJECTDIR}/RelocSimulation.o \
	${OBJECTDIR}/StatusBar.o \
	${OBJECTDIR}/SystemParameters.o \
	${OBJECTDIR}/eval_time_dep_main.o


# C Compiler Flags
CFLAGS=-O3

# CC Compiler Flags
CCFLAGS=-O3
CXXFLAGS=-O3

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/relsys

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/relsys: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/relsys ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/Combinatorial.o: Combinatorial.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Combinatorial.o Combinatorial.cpp

${OBJECTDIR}/Customer.o: Customer.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Customer.o Customer.cpp

${OBJECTDIR}/CustomerData.o: CustomerData.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/CustomerData.o CustomerData.cpp

${OBJECTDIR}/EntireSystem.o: EntireSystem.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/EntireSystem.o EntireSystem.cpp

${OBJECTDIR}/Experiments.o: Experiments.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/Experiments.o Experiments.cpp

${OBJECTDIR}/HeuristicQueue.o: HeuristicQueue.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/HeuristicQueue.o HeuristicQueue.cpp

${OBJECTDIR}/HyperQueue.o: HyperQueue.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/HyperQueue.o HyperQueue.cpp

${OBJECTDIR}/LinSolver.o: LinSolver.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/LinSolver.o LinSolver.cpp

${OBJECTDIR}/PhaseFitter.o: PhaseFitter.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/PhaseFitter.o PhaseFitter.cpp

${OBJECTDIR}/QueueData.o: QueueData.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/QueueData.o QueueData.cpp

${OBJECTDIR}/QueuePerformance.o: QueuePerformance.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/QueuePerformance.o QueuePerformance.cpp

${OBJECTDIR}/RelocEvaluation.o: RelocEvaluation.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/RelocEvaluation.o RelocEvaluation.cpp

${OBJECTDIR}/RelocSimulation.o: RelocSimulation.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/RelocSimulation.o RelocSimulation.cpp

${OBJECTDIR}/StatusBar.o: StatusBar.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/StatusBar.o StatusBar.cpp

${OBJECTDIR}/SystemParameters.o: SystemParameters.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/SystemParameters.o SystemParameters.cpp

${OBJECTDIR}/eval_time_dep_main.o: eval_time_dep_main.cpp
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/eval_time_dep_main.o eval_time_dep_main.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
