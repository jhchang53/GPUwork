NVCC=/usr/local/cuda-12.4/bin/nvcc
CUDAH = /usr/local/cuda-12.4/targets/x86_64-linux/include
CUDAL = /usr/local/cuda-12.4/targets/x86_64-linux/lib

CXX = nvcc
CXXFLAGS = -g -rdc=true -G -I$(CUDAH)

matgpu: matgpu.o Mater.o Mathost.o Matdev.o XSLIB.o  XSiso.o XSiso_collect.o
	nvcc -rdc=true -o matgpu -g -G $^
matlib: matlib.o Mater.o XSLIB.o  XSiso.o XSiso_collect.o
	nvcc -rdc=true -o matlib -g -G $^
mater: mater.o Mater.o XSLIB.o  XSiso.o XSiso_collect.o
	nvcc -rdc=true -o mater -g -G $^
xslib: xslib.o XSLIB.o XSiso.o XSiso_collect.o \
  XSLhost.o XSLhost2.o XSLhost3.o XSLhost4.o XSLdev.o
	nvcc -rdc=true -o xslib -g -G $^
#       g++ -o xslib -g $^ -lm
# xslib: xslib.o XSLIB.o XSiso.o XSiso_collect.o  XSLIB_host.cpp
#	g++ -o xslib -g $^ -lm
clear:
	rm -f *.o
%.o: %.cu
	nvcc -g -rdc=true -G -c $< -I $(CUDAH)

