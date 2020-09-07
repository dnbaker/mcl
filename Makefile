

mclcomp: mclcomp.cpp
	$(CXX) -Iblaze -march=native -std=c++17 $< -o $@
