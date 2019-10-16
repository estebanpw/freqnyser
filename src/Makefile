CC=gcc
CXX=g++
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE -Wall #-DVERBOSE
BIN=../bin

all: freqgen massive-gen compare-gen



freqgen: freqgen.c
	$(CXX) $(CFLAGS) commonFunctions.c -lm freqgen.c -lm -o $(BIN)/freqgen

massive-gen: massive-gen.c
	$(CXX) $(CFLAGS) commonFunctions.c -lm massive-gen.c -lm -o $(BIN)/massive-gen

compare-gen: compare-gen.c
	$(CXX) $(CFLAGS) commonFunctions.c -lm compare-gen.c -lm -o $(BIN)/compare-gen


clean:
	rm -rf $(BIN)/freqgen $(BIN)/massive-gen $(BIN)/compare-gen
