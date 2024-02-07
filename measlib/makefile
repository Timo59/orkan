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
CIRCLIB = ../circlib
MEASLIB = .
LAPACK = /usr/local/opt/lapack
CBLAS = /usr/local/opt/openblas

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

TARGET_BASE = testExact
TARGET = $(TARGET_BASE)$(TARGET_EXTENSION)
INC_PATHS = -I$(QLIB)/include -I$(STATELIB)/include -I$(GATELIB)/include -I$(CIRCLIB)/include \
-I$(MEASLIB)/include -I$(CBLAS)/include -I$(LAPACK)/include
LIB_PATHS = -L$(QLIB) -llinalg -L$(STATELIB) -lpstate -L$(CIRCLIB) -ldiag -L$(CIRCLIB) -lpauli \
-L$(CIRCLIB) -loperator -L$(CBLAS)/lib -L$(LAPACK)/lib
LDFLAGS = -L$(LAPACK)/lib -L$(CBLAS)/lib
SRC_FILES = $(UNITY_ROOT)/unity.c $(QLIB)/src/cplxutil.c $(QLIB)/src/linalg.c \
$(STATELIB)/src/pstate.c $(GATELIB)/src/qhipster.c $(CIRCLIB)/src/operator.c $(CIRCLIB)/src/diag.c \
$(CIRCLIB)/src/pauli.c $(MEASLIB)/src/exact.c

SOURCES_DIR = $(MEASLIB)/test
SOURCE_FILES = $(wildcard $(SOURCES_DIR)/*.c)
FILENAMES = $(patsubst $(SOURCES_DIR)/%.c,%,$(SOURCE_FILES))

default:
	@$(C_COMPILER) $(CFLAGS) -c -fPIC $(MEASLIB)/src/exact.c $(INC_PATHS)

lib:
	@$(MAKE) -C $(QLIB)
	@$(MAKE) -C $(STATELIB)
	@$(MAKE) -C $(GATELIB)
	@$(MAKE) -C $(CIRCLIB)
	@$(C_COMPILER) $(CFLAGS) -c -fPIC $(MEASLIB)/src/exact.c $(INC_PATHS)
	@$(C_COMPILER) -shared -o libmeas.so exact.o $(QLIB)/linalg.o $(STATELIB)/pstate.o \
	$(GATELIB)/qhipster.o $(CIRCLIB)/diag.o $(CIRCLIB)/pauli.o $(CIRCLIB)/operator.o \
	-install_name @rpath/libmeas.so
	@$(MAKE) -C $(QLIB) cleanobject
	@$(MAKE) -C $(STATELIB) cleanobject
	@$(MAKE) -C $(GATELIB) cleanobject
	@$(MAKE) -C $(CIRCLIB) cleanobject
	@$(CLEANUP) *.o

unittest: $(FILENAMES)

.SECONDARY:

$(FILENAMES):
	@$(C_COMPILER) $(CFLAGS) $(SRC_FILES) $(MEASLIB)/test/$@.c -o $@.$(EXTENSION) $(INC_PATHS) \
	-I$(QLIB)/test -I$(UNITY_ROOT)
	- ./$@.$(EXTENSION)

cleanobject:
	@$(CLEANUP) $(MEASLIB)/*.o

cleantest:
	@$(CLEANUP) $(MEASLIB)/*.$(EXTENSION)

cleanlib:
	@$(CLEANUP) $(MEASLIB)/*.so

cleanall:
	@$(CLEANUP) $(MEASLIB)/*.o
	@$(CLEANUP) $(MEASLIB)/*.$(EXTENSION)
	@$(CLEANUP) $(MEASLIB)/*.so

ci: CFLAGS += -Werror
ci: default
