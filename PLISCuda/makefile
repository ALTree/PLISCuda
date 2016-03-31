
NVCC = nvcc
OBJDIR = build
SRCDIR = src
NVCCFLAGS = -std=c++11 -rdc=true -gencode arch=compute_20,code=compute_20 -gencode arch=compute_20,code=sm_20
GCCFLAGS = -std=c++11 

objs = $(OBJDIR)/leap.o $(OBJDIR)/rates.o $(OBJDIR)/log.o \
	$(OBJDIR)/ssa.o $(OBJDIR)/driver.o $(OBJDIR)/main.o

##  build configurations  ##

# relase simulates the system (no print-logging)
release: NVCCFLAGS += -O3
release: GCCFLAGS += -O3
release: all

# log simulates while print-logging ALL events
log: NVCCFLAGS += -O2 -DLOG
log: GCCFLAGS += -O2 -DLOG
log: all

# debugs logs everything and prints all the internal data
debug: NVCCFLAGS += -O0 --compiler-options -Wall -DLOG -DDEBUG
debug: GCCFLAGS += -O0 -Wall -DLOG -DDEBUG
debug: all




##  rules  ##

all: $(objs)
	$(NVCC) $^ -o PLISCuda

# build object files for cuda code
$(OBJDIR)/%.o: $(SRCDIR)/%.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# build object file for main.cpp 
$(OBJDIR)/main.o: $(SRCDIR)/main.cpp
	$(NVCC) --compile -x c++ $(GCCFLAGS) -c $< -o $@

clean:
	@echo "Cleaning up..."
	rm -f $(OBJDIR)/*
	rm -f PLISCuda
