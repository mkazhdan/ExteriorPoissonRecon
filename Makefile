EXTERIOR_POISSON_RECON_TARGET=ExteriorPoissonRecon
EXTERIOR_POISSON_RECON_SOURCE=ExteriorPoissonRecon/ExteriorPoissonRecon.cpp
SAMPLE_TARGET=Sample
SAMPLE_SOURCE=Sample/Sample.cpp
STEREO_TARGET=Stereo
STEREO_SOURCE=Stereo/Stereo.cpp
DILATE_TARGET=Dilate
DILATE_SOURCE=Dilate/Dilate.cpp
EXTRACT_TARGET=Extract
EXTRACT_SOURCE=Extract/Extract.cpp

COMPILER ?= gcc
#COMPILER ?= clang

ifeq ($(COMPILER),gcc)
	CFLAGS += -fopenmp -Wno-deprecated -std=c++17 -Wno-invalid-offsetof
	LFLAGS += -lgomp -lstdc++ -lpthread
	CC=gcc
	CXX=g++
else
	CFLAGS += -Wno-deprecated -std=c++17 -Wno-invalid-offsetof -Wno-dangling-else -Wno-null-dereference
	LFLAGS += -lstdc++
	CC=clang
	CXX=clang++
endif

CFLAGS += -O3 -DRELEASE -funroll-loops -ffast-math -g
LFLAGS += -O3 -g

BIN = Bin/Linux/
BIN_O = Obj/Linux/
INCLUDE = ./Include/ -I .


MD=mkdir

EXTERIOR_POISSON_RECON_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(EXTERIOR_POISSON_RECON_SOURCE))))
EXTERIOR_POISSON_RECON_OBJECT_DIR=$(dir $(EXTERIOR_POISSON_RECON_OBJECTS))
SAMPLE_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(SAMPLE_SOURCE))))
SAMPLE_OBJECT_DIR=$(dir $(SAMPLE_OBJECTS))
STEREO_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(STEREO_SOURCE))))
STEREO_OBJECT_DIR=$(dir $(STEREO_OBJECTS))
DILATE_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(DILATE_SOURCE))))
DILATE_OBJECT_DIR=$(dir $(DILATE_OBJECTS))
EXTRACT_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(EXTRACT_SOURCE))))
EXTRACT_OBJECT_DIR=$(dir $(EXTRACT_OBJECTS))

all: make_dirs
all: $(BIN)$(EXTERIOR_POISSON_RECON_TARGET)
all: $(BIN)$(SAMPLE_TARGET)
all: $(BIN)$(STEREO_TARGET)
all: $(BIN)$(DILATE_TARGET)
all: $(BIN)$(EXTRACT_TARGET)

exteriorpoissonrecon: make_dirs
exteriorpoissonrecon: $(BIN)$(EXTERIOR_POISSON_RECON_TARGET)

sample: make_dirs
sample: $(BIN)$(SAMPLE_TARGET)

stereo: make_dirs
stereo: $(BIN)$(STEREO_TARGET)

dilate: make_dirs
dilate: $(BIN)$(DILATE_TARGET)

extract: make_dirs
extract: $(BIN)$(EXTRACT_TARGET)

clean:
	rm -rf $(BIN)$(EXTERIOR_POISSON_RECON_TARGET)
	rm -rf $(BIN)$(SAMPLE_TARGET)
	rm -rf $(BIN)$(STEREO_TARGET)
	rm -rf $(BIN)$(DILATE_TARGET)
	rm -rf $(BIN)$(EXTRACT_TARGET)
	rm -rf $(BIN_O)

make_dirs: FORCE
	$(MD) -p $(BIN)
	$(MD) -p $(BIN_O)
	$(MD) -p $(EXTERIOR_POISSON_RECON_OBJECT_DIR)
	$(MD) -p $(SAMPLE_OBJECT_DIR)
	$(MD) -p $(STEREO_OBJECT_DIR)
	$(MD) -p $(DILATE_OBJECT_DIR)
	$(MD) -p $(EXTRACT_OBJECT_DIR)

$(BIN)$(EXTERIOR_POISSON_RECON_TARGET): $(EXTERIOR_POISSON_RECON_OBJECTS)
	$(CXX) -o $@ $(EXTERIOR_POISSON_RECON_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN)$(SAMPLE_TARGET): $(SAMPLE_OBJECTS)
	$(CXX) -o $@ $(SAMPLE_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN)$(STEREO_TARGET): $(STEREO_OBJECTS)
	$(CXX) -o $@ $(STEREO_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN)$(DILATE_TARGET): $(DILATE_OBJECTS)
	$(CXX) -o $@ $(DILATE_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN)$(EXTRACT_TARGET): $(EXTRACT_OBJECTS)
	$(CXX) -o $@ $(EXTRACT_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN_O)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

FORCE: