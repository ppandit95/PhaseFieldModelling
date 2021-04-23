	
run : 1DTransient
	./1DTransient 
	
test/run : Tests/Test1
	cd Tests/ && ./Test1 --log_level=test_suite && cd ..
	
clean : Output/Initial_Profile.* Output/Output-*.* 1DTransient transient.o Tests/Test1
	rm Output/Initial_Profile.* Output/Output-*.* 1DTransient transient.o Tests/Test1
	
tests : Tests/Tests.cpp
	g++ -o Tests/Test1 Tests/Tests.cpp transient.o -L"/usr/lib" -lboost_unit_test_framework
	

	
1DTransient: main.cpp transient.o
	g++ -o 1DTransient main.cpp transient.o
	
transient.o : src/transient.cpp src/transient.hpp src/structures.hpp
	g++ -c src/transient.cpp
