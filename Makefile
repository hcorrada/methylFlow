DEPSDIR=${PWD}/deps

CPPFLAGS = 
CXXFLAGS = -g -O0 -Wall
LDFLAGS = -L${DEPSDIR} -lemon -lglpk
CXX = g++

OBJS = main.o MFGraph.o MethylRead.o MFGraph_solve.o MFRegionPrinter.o

all: methylFlow

test: run1 run2

main.o: main.cpp MFGraph.hpp
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c main.cpp

MFGraph.o: MFGraph.cpp MFGraph.hpp MethylRead.hpp MFRegionPrinter.hpp
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c MFGraph.cpp

MethylRead.o: MethylRead.cpp MethylRead.hpp 
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c MethylRead.cpp

MFGraph_solve.o: MFGraph_solve.cpp MFGraph.hpp MethylRead.hpp MFRegionPrinter.hpp 
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c MFGraph_solve.cpp

MFRegionPrinter.o: MFRegionPrinter.cpp MFRegionPrinter.hpp MFGraph.hpp MethylRead.hpp
	${CXX} ${CPPFLAGS} ${CXXFLAGS} -c MFRegionPrinter.cpp

test: MethylRead.o
	make -C testing all

clean:
	rm lemonTest *.o

methylFlow: Lemon.ts ${OBJS}
	${CXX} ${CXXFLAGS} ${LDFLAGS} -o methylFlow ${OBJS}

run1: methylFlow
	./methylFlow -i testing/sim1.tsv -o testing -l 1.0 -s 10.0

run2: methylFlow
	./methylFlow -i testing/sim2.tsv -o testing -l 1.0 20.0

Glpk.ts:
	(cd glpk && \
	./configure --prefix=${DEPSDIR} --disable-dependency-tracking && \
	make && \
	make install && \
	touch $@)

Lemon.ts: Glpk.ts
	(cd lemon && \
	mkdir build && \
	cd build && \
	cmake -DCMAKE_INSTALL_PREFIX=${DEPSDIR} -DLEMON_ENABLE_GLPK=YES -DLEMON_GLPK_ROOT_DIR=${DEPSDIR} .. && \
	make && \
	make install && \
	touch $@)

