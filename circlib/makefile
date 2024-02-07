ifeq ($(OS),Windows_NT)
  ifeq ($(shell uname -s),) # not in a bash-like shell
	CLEANUP = del /F /Q
	MKDIR = mkdir
  else # in a bash-like shell, like msys
	CLEANUP = rm -f
	MKDIR = mkdir -p
  endif
	EXTENSION=exe
else
	CLEANUP = rm -f
	MKDIR = mkdir -p
	EXTENSION=out
endif

C_COMPILER=gcc
ifeq ($(shell uname -s), Darwin)
C_COMPILER=clang
endif

UNITY_ROOT = ../../unitysrc
QLIB = ..
STATELIB = ../statelib
GATELIB = ../gatelib
CIRCLIB = .

CFLAGS = -std=c11
CFLAGS += -Wall
CFLAGS += -Wextra
CFLAGS += -Wpointer-arith
CFLAGS += -Wcast-align
CFLAGS += -Wwrite-strings
CFLAGS += -Wswitch-default
CFLAGS += -Wunreachable-code
CFLAGS += -Winit-self
CFLAGS += -Wmissing-field-initializers
CFLAGS += -Wno-unknown-pragmas
CFLAGS += -Wstrict-prototypes
CFLAGS += -Wundef
CFLAGS += -Wold-style-definition

INC_PATHS = -I$(QLIB)/include -I$(STATELIB)/include -I$(GATELIB)/include -I$(CIRCLIB)/include
SRC_FILES = $(UNITY_ROOT)/unity.c $(QLIB)/src/cplxutil.c $(QLIB)/src/linalg.c \
$(STATELIB)/src/pstate.c $(GATELIB)/src/qhipster.c $(CIRCLIB)/src/pauli.c $(CIRCLIB)/src/diag.c \
$(CIRCLIB)/src/operator.c

SOURCES_DIR = $(CIRCLIB)/test
SOURCE_FILES = $(wildcard $(SOURCES_DIR)/*.c)
FILENAMES = $(patsubst $(SOURCES_DIR)/%.c,%,$(SOURCE_FILES))

default:
	@$(C_COMPILER) $(CFLAGS) -c -fPIC src/diag.c $(INC_PATHS)
	@$(C_COMPILER) $(CFLAGS) -c -fPIC src/pauli.c $(INC_PATHS)
	@$(C_COMPILER) $(CFLAGS) -c -fPIC src/operator.c $(INC_PATHS)

lib:
	@$(MAKE) -C $(QLIB)
	@$(MAKE) -C $(STATELIB)
	@$(MAKE) -C $(GATELIB)
	@$(C_COMPILER) $(CFLAGS) -c -fPIC $(CIRCLIB)/src/diag.c $(INC_PATHS)
	@$(C_COMPILER) $(CFLAGS) -c -fPIC $(CIRCLIB)/src/pauli.c $(INC_PATHS)
	@$(C_COMPILER) $(CFLAGS) -c -fPIC $(CIRCLIB)/src/operator.c $(INC_PATHS)
	@$(C_COMPILER) -shared -o libcirc.so operator.o pauli.o diag.o $(QLIB)/linalg.o \
	$(STATELIB)/pstate.o $(GATELIB)/qhipster.o -install_name @rpath/libcirc.so
	@$(MAKE) -C $(QLIB) cleanobject
	@$(MAKE) -C $(STATELIB) cleanobject
	@$(MAKE) -C $(GATELIB) cleanobject
	@$(CLEANUP) *.o

unittest: $(FILENAMES)

.SECONDARY:

$(FILENAMES):
	@$(C_COMPILER) $(CFLAGS) $(SRC_FILES) $(CIRCLIB)/test/$@.c -o $@.$(EXTENSION) $(INC_PATHS) \
	-I$(QLIB)/test -I$(UNITY_ROOT)
	- ./$@.$(EXTENSION)

cleanobject:
	@$(CLEANUP) $(CIRCLIB)/*.o

cleantest:
	@$(CLEANUP) $(CIRCLIB)/*.$(EXTENSION)

cleanlib:
	@$(CLEANUP) $(CIRCLIB)/*.so

cleanall:
	@$(CLEANUP) $(CIRCLIB)/*.o
	@$(CLEANUP) $(CIRCLIB)/*.$(EXTENSION)
	@$(CLEANUP) $(CIRCLIB)/*.so

ci: CFLAGS += -Werror
ci: default
