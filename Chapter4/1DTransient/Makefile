1DTransient: main.cpp transient.o
	g++ -o 1DTransient main.cpp transient.o
	
transient.o : transient.cpp transient.hpp structures.hpp
	g++ -c src/transient.cpp
	
run : 1DTransient
	./1DTransient 
	
clean : Initial_Profile.* Output-*.* 1DTransient transient.o
	rm Output/Initial_Profile.* Output/Output-*.* 1DTransient transient.o
