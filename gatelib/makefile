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
GATELIB = .

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

INC_PATHS = -I$(QLIB)/include -I$(STATELIB)/include -I$(GATELIB)/include
SRC_FILES = $(UNITY_ROOT)/unity.c $(QLIB)/src/cplxutil.c $(QLIB)/src/linalg.c \
$(STATELIB)/src/pstate.c $(GATELIB)/src/qhipster.c

SOURCES_DIR = $(GATELIB)/test
SOURCE_FILES = $(wildcard $(SOURCES_DIR)/*.c)
FILENAMES = $(patsubst $(SOURCES_DIR)/%.c,%,$(SOURCE_FILES))

default:
	@$(C_COMPILER) $(CFLAGS) -c -fPIC $(GATELIB)/src/qhipster.c $(INC_PATHS)

lib:
	@$(MAKE) -C $(STATELIB)
	@$(C_COMPILER) $(CFLAGS) -c -fPIC $(GATELIB)/src/qhipster.c $(INC_PATHS)
	@$(C_COMPILER) -shared -o libgate.so qhipster.o $(STATELIB)/pstate.o \
	-install_name @rpath/libgate.so
	@$(MAKE) -C $(STATELIB) cleanobject
	@$(CLEANUP) *.o

unittest: $(FILENAMES)

.SECONDARY:

$(FILENAMES):
	@$(C_COMPILER) $(CFLAGS) $(SRC_FILES) $(GATELIB)/test/$@.c -o $@.$(EXTENSION) $(INC_PATHS) \
	-I$(QLIB)/test -I$(UNITY_ROOT)
	- ./$@.$(EXTENSION)

cleanobject:
	@$(CLEANUP) $(GATELIB)/*.o

cleantest:
	@$(CLEANUP) $(GATELIB)/*.$(EXTENSION)

cleanlib:
	@$(CLEANUP) $(GATELIB)/*.so

cleanall:
	@$(CLEANUP) $(GATELIB)/*.o
	@$(CLEANUP) $(GATELIB)/*.$(EXTENSION)
	@$(CLEANUP) $(GATELIB)/*.so

ci: CFLAGS += -Werror
ci: default
