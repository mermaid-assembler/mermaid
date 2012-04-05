PROG_NAME = mermaid

CXX = mpic++
CXXFLAGS = -c -g -Wall -O3 -fshort-enums
LDFLAGS = 
ifeq ($(shell uname), Darwin)
LIBS = -lboost_mpi -lboost_serialization-mt -lboost_filesystem-mt -lboost_system-mt
else
LIBS = -lboost_mpi -lboost_serialization -lboost_filesystem -lboost_system
endif
DEPFLAGS = -MM

SRCS := $(wildcard *.cpp)
OBJS := $(SRCS:.cpp=.o)
DEPS := $(OBJS:.o=.d)

all : $(PROG_NAME)

$(PROG_NAME) : $(OBJS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)

%.d : %.cpp
	@$(CXX) $(DEPFLAGS) -o $@ $^
#	@echo -e '\t$$(CXX) $$(CXXFLAGS) -o $$@ $$<' >> $@

-include $(DEPS)

.PHONY : clean

clean :
	rm -rf *.o *.d $(PROG_NAME)
