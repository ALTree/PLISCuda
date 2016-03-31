
NVCC = nvcc
OBJDIR = build
SRCDIR = src
NVCCFLAGS = -std=c++11 -gencode="arch=compute_30, code=sm_30"

objs = $(OBJDIR)/leap.o $(OBJDIR)/rates.o $(OBJDIR)/log.o \
	$(OBJDIR)/ssa.o $(OBJDIR)/driver.o $(OBJDIR)/main.o

##  build configurations  ##

debug: NVCCFLAGS += -O0 --compiler-options -Wall
debug: all

release: NVCCFLAGS += -O3
release: all


##  rules  ##

all: $(objs)
	$(NVCC) $^ -o PLISCuda

# build object files for cuda code
$(OBJDIR)/%.o: $(SRCDIR)/%.cu
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

# build object file for main.cpp 
$(OBJDIR)/main.o: $(SRCDIR)/main.cpp
	c++ -std=c++11 -c $< -o $@

clean:
	@echo "Cleaning up..."
	rm -f $(OBJDIR)/*
	rm -f PLISCuda
