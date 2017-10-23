#SYSTEM = FX10
#SYSTEM = INTEL
#SYSTEM = XC30
SYSTEM = PGI

OPT = 8

#FX10
ifeq ($(SYSTEM),FX10)
OPTFLAGS = -fs
CC=mpifccpx
F90=mpifrtpx -Kfast,openmp
#F90=mpifrtpx -Kopenmp
CCFLAGS = $(OPTFLAGS)
F90FLAGS = $(OPTFLAGS) -Cfpp
LDFLAGS = -SSL2
endif

#intel
ifeq ($(SYSTEM),INTEL)
#OPTFLAGS = -O3 -traceback -ip -heap-arrays -qopenmp
OPTFLAGS = -qopenmp -O3 -ip
CC=mpiicc
F90=mpiifort
CCFLAGS = $(OPTFLAGS)
#F90FLAGS = $(OPTFLAGS) -fpp -assume nounderscore -names uppercase
F90FLAGS = $(OPTFLAGS) -fpp
#F90FLAGS = $(OPTFLAGS) -fpp -check all
#F90FLAGS = -fpe0 -traceback -g -CB -assume nounderscore -names lowercase -fpp -check all
#LDFLAGS = -mkl -trace
LDFLAGS = -mkl
endif

#XC30
ifeq ($(SYSTEM),XC30)
OPTFLAGS = -O2 -homp
CC=cc
F90=ftn
CCFLAGS = $(OPTFLAGS)
F90FLAGS = $(OPTFLAGS)
endif

#PGI
ifeq ($(SYSTEM),PGI)

#MODE = ACC_UNIFIED
#MODE = ACC
#MODE = ACC_CUDA_UNIFIED
MODE = ACC_CUDA
GPUOPTIMIZED=1
ALIGNED=0
BLOCKING=0

ifeq ($(MODE),ACC_UNIFIED)
OFLAGS = -acc -ta=tesla:cc60,managed -O3 -fastsse -Minfo=accel -mcmodel=medium
LFLAGS = -acc -ta=tesla:cc60,managed
endif
ifeq ($(MODE),ACC)
OFLAGS = -acc -ta=tesla:cc60 -O3 -fastsse -Minfo=accel -mcmodel=medium
LFLAGS = -acc -ta=tesla:cc60
endif
ifeq ($(MODE),ACC_CUDA_UNIFIED)
OFLAGS = -acc -ta=tesla:cc60,managed -O3 -fastsse -Minfo=accel -Mcuda=loadcache:L2,ptxinfo,fastmath,fma -DUSECUDA -mcmodel=medium 
LFLAGS = -acc -ta=tesla:cc60,managed -Mcuda 
OBJ2 = m_ppohBEM_user_func_cuda.o m_HACApK_base_cuda.o 
endif
ifeq ($(MODE),ACC_CUDA)
#OFLAGS = -acc -ta=tesla:cc60 -O3 -fastsse -Minfo=accel -Mcuda=loadcache:L2,ptxinfo,fastmath,fma,maxregcount:168 -mcmodel=medium -DUSECUDA #-Minline=name:face_integral_cuda,cross_product
#OFLAGS = -acc -ta=tesla:cc60 -O0 -Minfo=accel -Mcuda=loadcache:L2,ptxinfo,fastmath,fma,maxregcount:168 -mcmodel=medium -DUSECUDA
#OFLAGS = -acc -ta=tesla:cc60 -O3 -Minfo=accel -Mcuda=O0,ptxinfo,cc60,maxregcount:168 -mcmodel=medium -DUSECUDA
OFLAGS = -acc -ta=tesla:cc60 -O3 -Minfo=accel -Mcuda=ptxinfo,cc60 -mcmodel=medium -DUSECUDA -DOPT=${OPT}
LFLAGS = -acc -ta=tesla:cc60 -Mcuda
#OBJ2 = m_ppohBEM_user_func_cuda.o m_HACApK_base_cuda.o
OBJ2 = m_HACApK_base_cuda.o
endif

OPTFLAGS = $(OFLAGS) -DGPUOPTIMIZED=$(GPUOPTIMIZED) -DALIGNED=$(ALIGNED) -DBLOCKING=$(BLOCKING)

CC=mpicc
F90=mpif90
CCFLAGS = $(OPTFLAGS)
F90FLAGS = $(OPTFLAGS) -Mpreprocess
MKLPATH=/lustre/pz0108/z30108/ppohBEM_0.5.0/HACApK_1.0.0/src/HACApK_with_BEM-BB-framework_1.0.0/MKL_PGI
LDFLAGS = $(LFLAGS) -I${MKLPATH}/include -L${MKLROOT}/lib/intel64 -lmkl_blas95_lp64
endif


LINK=$(F90)

OBJS= HACApK_lib.o m_ppohBEM_user_func.o m_ppohBEM_matrix_element_ij.o m_HACApK_calc_entry_ij.o $(OBJ2)\
	 m_HACApK_base.o m_HACApK_solve_cuda.o m_HACApK_solve.o m_HACApK_use.o m_ppohBEM_bembb2hacapk.o bem-bb-fw-HACApK-0.5.0.o \

TARGET=bem-bb-SCM${OPT}.out

.SUFFIXES: .o .c .f90 .cuf

$(TARGET): $(OBJS)
			$(LINK) -o $@ $(OBJS) $(LDFLAGS)

.c.o: *.c
			$(CC) -c $(CCFLAGS) $<
.f90.o: *.f90
#			echo 'f90 complile'
			$(F90) -c $< $(F90FLAGS)

.cuf.o: *.cuf
#			echo 'f90 complile'
			$(F90) -c $< $(F90FLAGS)
clean:
#	rm -f *.o *.mod $(TARGET) *~
	rm -f *.o *.mod bem-bb-SCM.out *~

rmod:
	rm -f m_*.o *.mod

