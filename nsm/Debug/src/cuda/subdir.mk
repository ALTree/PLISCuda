################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CU_SRCS += \
../src/cuda/driver.cu \
../src/cuda/leap.cu \
../src/cuda/log.cu \
../src/cuda/rates.cu \
../src/cuda/ssa.cu 

CU_DEPS += \
./src/cuda/driver.d \
./src/cuda/leap.d \
./src/cuda/log.d \
./src/cuda/rates.d \
./src/cuda/ssa.d 

OBJS += \
./src/cuda/driver.o \
./src/cuda/leap.o \
./src/cuda/log.o \
./src/cuda/rates.o \
./src/cuda/ssa.o 


# Each subdirectory must supply rules for building sources it contributes
src/cuda/%.o: ../src/cuda/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/bin/nvcc -G -g -lineinfo -O0 -ccbin /usr/bin/gcc-4.8 -gencode arch=compute_30,code=sm_30  -odir "src/cuda" -M -o "$(@:%.o=%.d)" "$<"
	../nvcc.hack -G -g -lineinfo -O0 -ccbin /usr/bin/gcc-4.8 --compile --relocatable-device-code=true -gencode arch=compute_30,code=compute_30 -gencode arch=compute_30,code=sm_30  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


