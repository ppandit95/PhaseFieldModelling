all : main.cpp
	g++ -o 1DTransient main.cpp
	
run : 1DTransient
	./1DTransient
	
clean : Initial_Profile.* Output-*.* 1DTransient
	rm Initial_Profile.* Output-*.* 1DTransient
