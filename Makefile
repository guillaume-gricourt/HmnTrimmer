TARGET=HmnTrimmer

DIRSRC=./src
DIROBJ=./obj
DIRLIB=./lib
DIRTEST=./test

SRC=$(wildcard $(DIRSRC)/*.cpp)
OBJ=$(patsubst %.cpp,$(DIROBJ)/%.o,$(notdir $(SRC)))

CXX=g++
#/save/ggricourt/opt/gcc-5.3.0/build/bin/g++
LIBS=-I lib/seqan-2.4.0/include -I lib/spdlog-1.5.0/include -I lib/rapidjson-1.1.0/include -I lib/igzip-042/igzip/c_code -I lib/igzip-042/include -L lib/igzip-042/igzip
CXXFLAGS=-std=c++14 -O3 -W -Wall -pedantic -lrt -DNDEBUG -DSEQAN_ENABLE_DEBUG=0 -DSEQAN_ENABLE_TESTING=0 -DSEQAN_HAS_ZLIB=1 -lz -DSEQAN_HAS_OPENMP=1 -lgomp -fopenmp -lpthread -ligzip0c -DSEQAN_HAS_IGZIP=1

FTEST=$(DIRTEST)/run_tests.py

.PHONY: test clean

###########
## Rules ##
###########
all: clean $(TARGET) 

$(TARGET):$(OBJ)
	$(CXX) $(OBJ) $(LIBS) $(CXXFLAGS) -o $@

$(DIROBJ)/%.o:$(DIRSRC)/%.cpp
	@[ -d $(DIROBJ) ] || mkdir $(DIROBJ)
	$(CXX) $(LIBS) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(TARGET)
	$(RM) -r obj

test:
	@[ -x $(FTEST) ] || chmod +x $(FTEST)
	$(FTEST)
