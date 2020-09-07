

all: mclcomp

mclcomp: mclcomp.cpp mcl.h
	$(CXX) -Iblaze -march=native -std=c++17 $< -o $@
