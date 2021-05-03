SpinodalDecomposition.o : src/SpinodalDecomposition.cpp src/SpinodalDecompositio.hpp
	g++ -g -O2 -Wall -c src/SpinodalDecomposition.cpp
CaseStudyI : main.cpp
	g++ -g -O2 -Wall -o CaseStudyI main.cpp SpinodalDecomposition.o
	
run : CaseStudyI
	./CaseStudyI 
	
clean : CaseStudyI Output/* Results/*
	rm CaseStudyI Output/* Results/*