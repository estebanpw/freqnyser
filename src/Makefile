CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall #-DVERBOSE
BIN=../bin

all: freqgen



freqgen: freqgen.c
	$(CXX) $(CFLAGS) commonFunctions.c -lm freqgen.c -lm -o $(BIN)/freqgen

clean:
	rm -rf $(BIN)/freqgen
