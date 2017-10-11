# Makefile for project

OS := $(shell uname)
HOST := $(shell hostname)
INC_PATH = includes
LIB_PATH = lib
OBJ_PATH = obj
TEST_PATH = tests

TEST_FILES = $(wildcard $(TEST_PATH)/*.hpp)
SOURCE_FILES = $(wildcard $(LIB_PATH)/*.cpp)
OBJS = $(addprefix $(OBJ_PATH)/, $(notdir $(SOURCE_FILES:.cpp=.o)))
CPPFLAGS = -Ofast -I${MKLROOT}/include -ipo

# Check if mac
ifeq ($(OS),Darwin)
LDFLAGS =  -L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -ipo
else
# Check if ARCHER
ifneq (,$(findstring eslogin, $(HOST)))
LDFLAGS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -ipo
else
LDFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl -ipo
endif
endif

default: $(OBJS)
	$(CC) $(OBJS) -o run $(LDFLAGS)

$(OBJ_PATH)/%.o: $(LIB_PATH)/%.cpp
	$(CC) $(CPPFLAGS) -c -o $@ $<

clean:
	-rm -f $(OBJS)
	-rm -f run
