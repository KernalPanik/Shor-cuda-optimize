all: program

program: cudacode.o
	nvcc -o program -I/usr/local/cuda-10.2/include -L/usr/local/cuda-10.2/lib64 -lcuda -lcudart -lcufft main.cpp unitary.cpp Shor.cpp quantumRegister.cpp util.cpp rand.cpp cudacode.o

cudacode.o:
	nvcc qft.cu -lcufft -I/usr/local/cuda-10.2/include -L/usr/local/cuda-10.2/lib64 -lcuda -lcudart


