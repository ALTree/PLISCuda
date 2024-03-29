
NVCC = nvcc
OBJDIR = build
SRCDIR = src
NVCCFLAGS = -std=c++11 -rdc=true -arch=sm_30

objs = $(OBJDIR)/leap.o $(OBJDIR)/rates.o $(OBJDIR)/log.o \
	$(OBJDIR)/ssa.o $(OBJDIR)/tau.o $(OBJDIR)/init.o \
	$(OBJDIR)/driver.o $(OBJDIR)/main.o 

##  build configurations  ##

# release simulates the system (no print-logging)
release: NVCCFLAGS += -O3
release: all

# log simulates while print-logging ALL events
log: NVCCFLAGS += -O2 -DLOG
log: all

# debugs logs everything and prints all the internal data
debug: NVCCFLAGS += -O0 -g -G -lineinfo --compiler-options -Wall -DLOG -DDEBUG
debug: all

# profile is the same as release but with a few profile flags
profile: NVCCFLAGS += -O3 -DPROFILE -lineinfo
profile: all

# test defines the TEST preprocessor variable
test: NVCCFLAGS += -O2 --compiler-options -Wall -DTEST
test: all


##  rules  ##

all: $(objs)
	$(NVCC) $^ -o pliscuda

# build object files for cuda code
$(OBJDIR)/%.o: $(SRCDIR)/%.cu
	mkdir -p build
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# build object file for main.cpp 
$(OBJDIR)/main.o: $(SRCDIR)/main.cpp
	$(NVCC) --compile -x c++ $(NVCCFLAGS) -c $< -o $@

clean:
	@echo "Cleaning up..."
	rm -f $(OBJDIR)/*
	rm -f pliscuda
