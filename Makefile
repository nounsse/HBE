# C source code
CSRC	= Utilities/io_png.c \
Utilities/mt19937ar.c \
initialization/randmt.c \
Utilities/io_exr.c 
# C++ source code
CXXSRC	= main.cpp \
hbe/hbe.cpp \
hbe/LibMatrix.cpp \
Utilities/LibImages.cpp \
Utilities/Utilities.cpp \
Utilities/io_pgm.cpp \
initialization/DataProvider.cpp \
initialization/PLE_lib.cpp \
initialization/util.cpp 

# all source code
SRC	= $(CSRC) $(CXXSRC)

INCPATH   = -I/usr/include/OpenEXR -I./eigen -I./newmat10 -I./Utilities -I./hbe -I./initialization

# C objects
COBJ	= $(CSRC:.c=.o)
# C++ objects
CXXOBJ	= $(CXXSRC:.cpp=.o)
# all objects
OBJ	= $(COBJ) $(CXXOBJ)
# binary target
BIN	= HBE

# C optimization flags
COPT	= -O3 -ftree-vectorize -funroll-loops
#COPT	= -ftree-vectorize -funroll-loops -g -pg

# C++ optimization flags
CXXOPT	= $(COPT)

# C compilation flags
CFLAGS	= $(COPT) -Wall -Wextra \
	-Wno-write-strings -ansi
# C++ compilation flags
CXXFLAGS	= $(CXXOPT) -Wall \
	-Wno-write-strings -Wno-deprecated -ansi -I/usr/include/OpenEXR -I./eigen -I./newmat10 -I./Utilities -I./hbe -I./initialization
# link flags
LDFLAGS	= -L./newmat10 -lnewmat -lpng -lm -lIlmImf -lHalf
LIBS    = -L./newmat10 -lnewmat -lpng -lm -lIlmImf -lHalf
#LDFLAGS	= -L./newmat10 -lnewmat -lpng -lm -lIlmImf -lHalf -pg

# use openMP with `make OMP=1`
ifdef OMP
CFLAGS	+= -fopenmp
CXXFLAGS	+= -fopenmp
LDFLAGS += -lgomp
else
CFLAGS	+= -Wno-unknown-pragmas
CXXFLAGS  += -Wno-unknown-pragmas
endif

# link all the object code
$(BIN): $(OBJ) $(LIBDEPS) libnewmat.a
	$(CXX) -o $@ $(OBJ) $(LDFLAGS)
	
# partial compilation of C source code
%.o: %.c %.h
	$(CC) -c -o $@  $< $(CFLAGS)
# partial compilation of C++ source code
%.o: %.cpp %.h
	$(CXX) -c -o $@  $< $(CXXFLAGS)

libnewmat.a:
	$(MAKE) -C newmat10 -f nm_gnu.mak

# housekeeping
clean	:
	$(RM) $(OBJ)
