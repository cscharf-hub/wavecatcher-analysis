CXX             =g++
LD              =g++

INCROOT         =$(shell root-config --incdir)

CPPFLAGS        =-I $(INCROOT)/ -I $(INCDIR)/
CXXFLAGS        =-fPIC -g -O2 -Wall -Wextra -Wshadow -pedantic -fopenmp

# Get the C++ version of ROOT
ROOT_CFLAGS := $(shell root-config --cflags)
CXX_STD := $(shell echo $(ROOT_CFLAGS) | grep -oP '(?<=-std=c\+\+)\d+')

# Default
ifeq ($(CXX_STD),)
CXX_STD := 20
endif

CXXFLAGS += -std=c++$(CXX_STD)

DOXYGEN_DIR = html

DICTB           =ReadRunDictUX
DICTH           =${DICTB}.h
DICT            =${DICTB}.cc
DICTO           =${DICTB}.o

HDRS            =src/PMT.h src/CosmicsBox.h src/FFT_WF.h src/Experimental.h src/ReadSampic.h src/Legacy_functions.h src/utils/Helpers.h src/utils/Filters.h src/ReadRun.h 
OBJS            =src/PMT.o src/CosmicsBox.o src/FFT_WF.o src/Experimental.o src/ReadSampic.o src/Legacy_functions.o src/utils/Helpers.o src/utils/Filters.o src/ReadRun.o $(DICTO)

LIBSLIN			=$(shell root-config --glibs)

SLL             =ReadRunLib.sl

#__________________________________________________________

.PHONY: docs upload clean-docs all clean-intermediate clean test_rr

all:		${OBJS} ${DICTO}
		${LD} -shared ${CXXFLAGS} -o ${SLL} ${OBJS} ${LIBSLIN}	
		$(MAKE) clean-intermediate
		@echo "Success!"

${DICT}:    ${HDRS}
		rootcint -rml=ReadRun -f ${DICT} -c ${HDRS} -I. misc/LinkDef.h

test_rr:
		root examples/read_exampledata.cc -b -q > tst_read_exampledata.txt 
		root examples/use_functions_wo_measurement.cc -b -q > tst_use_functions_wo_measurement.txt
		root examples/timing_example.cc -b -q > tst_timing_example.txt
		root examples/timing_example_rebin.cc -b -q > tst_timing_example_rebin.txt
		python examples/read_exampledata.py 0 1 > tst_read_exampledata.py.txt 

docs:
		@doxygen

# update web page: make upload-docs CERNBOX="/path/to/cernbox"
upload-docs:
		rsync -av --delete $(DOXYGEN_DIR)/ "$(CERNBOX)/"

clean-docs:
		rm -rf $(DOXYGEN_DIR)

clean-intermediate:
		@echo "Cleaning up files: ${OBJS} ${DICTB}"
		@rm ${OBJS} ${DICTB}.cc

clean:
		@echo "Cleaning up files: ${SLL}"
		@rm ${SLL}

%.o:        %.cc
		${LD} ${CPPFLAGS} ${CXXFLAGS} -c $< -o $@          