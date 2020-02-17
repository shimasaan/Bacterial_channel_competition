FILE = ./source/voter_phase.cpp
OUT = cell_growth_chemnet_MPI

all:$(FILE).cpp
	mpiicpc $(FILE).cpp -o $(FILE) -O3 -lm -fomit-frame-pointer -funroll-loops -ipo -std=c++11