ifeq ($(OS),Windows_NT)
  ifeq ($(shell uname -s),) # not in a bash-like shell
	CLEANUP = del /F /Q
	MKDIR = mkdir
  else # in a bash-like shell, like msys
	CLEANUP = rm -f
	MKDIR = mkdir -p
  endif
	TARGET_EXTENSION=exe
else
	CLEANUP = rm -f
	MKDIR = mkdir -p
	TARGET_EXTENSION=out
endif

C_COMPILER=gcc
ifeq ($(shell uname -s), Darwin)
C_COMPILER=clang
endif

UNITY_ROOT = ../../unitysrc
QLIB = ..
STATELIB = .

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

INC_PATHS = -I$(QLIB)/include -I$(STATELIB)/include

default:
	@$(C_COMPILER) $(CFLAGS) -c -fPIC $(STATELIB)/src/pstate.c $(INC_PATHS)

lib:
	@$(C_COMPILER) $(CFLAGS) -c -fPIC $(STATELIB)/src/pstate.c $(INC_PATHS)
	@$(C_COMPILER) -shared -o libstate.so pstate.o -install_name @rpath/libstate.so
	@$(CLEANUP) pstate.o

unittest:

cleanobject:
	@$(CLEANUP) $(STATELIB)/*.o

cleantest:
	@$(CLEANUP) $(STATELIB)/*.$(EXTENSION)

cleanlib:
	@$(CLEANUP) $(STATELIB)/*.so

cleanall:
	@$(CLEANUP) $(STATELIB)/*.o
	@$(CLEANUP) $(STATELIB)/*.$(EXTENSION)
	@$(CLEANUP) $(STATELIB)/*.so

ci: CFLAGS += -Werror
ci: default
