makesolution:
	g++ Solution.cpp -o Solution.o -c -Wall -DPRINT -std=c++14 -O2 
	g++ molecule.cpp -o molecule.o -c -Wall -DPRINT -std=c++14 -O2 
	g++ -o P1Solution Solution.o molecule.o -O2 -lm -llapack
