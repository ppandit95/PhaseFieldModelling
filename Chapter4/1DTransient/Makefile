1DTransient: main.cpp transient.o
	g++ -o 1DTransient main.cpp transient.o
	
transient.o : src/transient.cpp src/transient.hpp src/structures.hpp
	g++ -c src/transient.cpp
	
run : 1DTransient
	./1DTransient 
	
clean : Output/Initial_Profile.* Output/Output-*.* 1DTransient transient.o Tests/Test1
	rm Output/Initial_Profile.* Output/Output-*.* 1DTransient transient.o Tests/Test1
	
tests : Tests/Tests.cpp
	g++ -o Tests/Test1 Tests/Tests.cpp -L"/usr/lib" -lboost_unit_test_framework
