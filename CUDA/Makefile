NVCC=/usr/local/cuda-12.4/bin/nvcc
CUDAH = /usr/local/cuda-12.4/targets/x86_64-linux/include
CUDAL = /usr/local/cuda-12.4/targets/x86_64-linux/lib

CXX = nvcc
CXXFLAGS = -G -I$(CUDAH)

testrk: testrk.o track.o CGhost.o CGdev.o SRChost.o SRCdev.o GEOMdev.o BOX.o
	nvcc -rdc=true -o testrk -G $^
rantest: rantest.o
	nvcc -o rantest -g $^
b10kerma: blkdata.o B10kerma.o
	nvcc -rdc=true -o b10kerma -G $^
clear:
	rm -f *.o
%.o: %.cu
	nvcc -rdc=true -G -c $< -I $(CUDAH)

