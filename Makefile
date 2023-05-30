EXTERIOR_POISSON_RECON_TARGET=ExteriorPoissonRecon
EXTERIOR_POISSON_RECON_SOURCE=ExteriorPoissonRecon/ExteriorPoissonRecon.cpp
SAMPLE_TARGET=Sample
SAMPLE_SOURCE=Sample/Sample.cpp
VISUALIZE_3D_TARGET=Visualize3D
VISUALIZE_3D_SOURCE=Visualize3D/Visualize3D.cpp
VISUALIZE_4D_TARGET=Visualize4D
VISUALIZE_4D_SOURCE=Visualize4D/Visualize4D.cpp

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
VISUALIZE_3D_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(VISUALIZE_3D_SOURCE))))
VISUALIZE_3D_OBJECT_DIR=$(dir $(VISUALIZE_3D_OBJECTS))
VISUALIZE_4D_OBJECTS=$(addprefix $(BIN_O), $(addsuffix .o, $(basename $(VISUALIZE_4D_SOURCE))))
VISUALIZE_4D_OBJECT_DIR=$(dir $(VISUALIZE_4D_OBJECTS))

all: make_dirs
all: $(BIN)$(EXTERIOR_POISSON_RECON_TARGET)
all: $(BIN)$(SAMPLE_TARGET)
all: $(BIN)$(VISUALIZE_3D_TARGET)
all: $(BIN)$(VISUALIZE_4D_TARGET)

exteriorpoissonrecon: make_dirs
exteriorpoissonrecon: $(BIN)$(EXTERIOR_POISSON_RECON_TARGET)

sample: make_dirs
sample: $(BIN)$(SAMPLE_TARGET)

visualize3d: make_dirs
visualize3d: $(BIN)$(VISUALIZE_3D_TARGET)

visualize4d: make_dirs
visualize4d: $(BIN)$(VISUALIZE_4D_TARGET)

clean:
	rm -rf $(BIN)$(EXTERIOR_POISSON_RECON_TARGET)
	rm -rf $(BIN)$(SAMPLE_TARGET)
	rm -rf $(BIN)$(VISUALIZE_3D_TARGET)
	rm -rf $(BIN)$(VISUALIZE_4D_TARGET)
	rm -rf $(BIN_O)

make_dirs: FORCE
	$(MD) -p $(BIN)
	$(MD) -p $(BIN_O)
	$(MD) -p $(EXTERIOR_POISSON_RECON_OBJECT_DIR)
	$(MD) -p $(SAMPLE_OBJECT_DIR)
	$(MD) -p $(VISUALIZE_3D_OBJECT_DIR)
	$(MD) -p $(VISUALIZE_4D_OBJECT_DIR)

$(BIN)$(EXTERIOR_POISSON_RECON_TARGET): $(EXTERIOR_POISSON_RECON_OBJECTS)
	$(CXX) -o $@ $(EXTERIOR_POISSON_RECON_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN)$(SAMPLE_TARGET): $(SAMPLE_OBJECTS)
	$(CXX) -o $@ $(SAMPLE_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN)$(VISUALIZE_3D_TARGET): $(VISUALIZE_3D_OBJECTS)
	$(CXX) -o $@ $(VISUALIZE_3D_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN)$(VISUALIZE_4D_TARGET): $(VISUALIZE_4D_OBJECTS)
	$(CXX) -o $@ $(VISUALIZE_4D_OBJECTS) -L$(BIN) $(LFLAGS)

$(BIN_O)%.o: $(SRC)%.cpp
	$(CXX) -c -o $@ $(CFLAGS) -I$(INCLUDE) $<

FORCE: