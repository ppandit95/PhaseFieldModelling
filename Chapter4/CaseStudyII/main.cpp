/*
 * CaseStudyII.cpp
 *
 *  Created on: 21-May-2021
 *      Author: pushkar
 */
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////   CASE STUDY II :-> Phase Field Modelling of Grain Growth with Finite Difference Algorithm     ////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<filesystem>
#include <unistd.h>
#include<fstream>
#include<cassert>
#include<string>
#include<vector>
#include<cassert>
#include<stdio.h>
int main(int argc,char** argv){

	//Get initial wall time
	struct timespec begin,end;
	clock_gettime(CLOCK_REALTIME,&begin);
	//Simulation Cell Parameters
	unsigned int Nx = 512;
	unsigned int Ny = 512;
	unsigned int Nz = 1;
	float dx = 1.0;
	float dy = 1.0;

	//Time Integration parameters
	unsigned int nstep = 100;
	float dtime = 0.005;
	float ttime = 0.0;

	//Material Parameters
	float mobil = 5.0;
	float grcoef = 0.1;

	//Task1 : Take grain region from input file and modify the corresponding array of order parameters
	//Task 1.1:Setup the grid points so that its convenient to compare grid point with the equation of edges

	//Setting up the grid
	double xgrid[Nx],ygrid[Ny];
	for(unsigned int i=0;i<Nx;i++)
		xgrid[i] = i*dx;
	for(unsigned int j=0;j<Ny;j++)
		ygrid[j] = j*dy;

	//Read the input data from input folder

	/*std::ifstream input_file("voronoi_input.inp");
	double voronoi_edge_x,voronoi_edge_y;
	unsigned int lines=0,n = 25;
	std::string buffer;
	std::vector<double> x_edge[n],y_edge[n];
	assert(input_file.is_open());
	input_file.seekg(8,input_file.cur);
	while(!input_file.eof()){
		input_file >> voronoi_edge_x;
		input_file >> voronoi_edge_y;
		x_edge[lines].push_back(voronoi_edge_x);
		y_edge[lines].push_back(voronoi_edge_y);
		//input_file.seekg(1,input_file.cur);
		std::getline(input_file,buffer);
		if(buffer.length()==0)
			lines++;
	}
	std::cout << lines << std::endl;
	input_file.seekg(8,input_file.cur);
	input_file >> voronoi_edge_x;
	input_file >> voronoi_edge_y;
	x_edge[0].push_back(voronoi_edge_x);
	y_edge[0].push_back(voronoi_edge_y);
	input_file.seekg(1,input_file.cur);*/

	double voronoi_edge_x,voronoi_edge_y;
	std::vector<double> x_edge[25],y_edge[25];
	unsigned int grainNo = 1,count;
	std:: string line;
	while(grainNo < 26){
		std::string filename = "/media/pushkar/Data/Personal Files/PhaseFieldModelling/Chapter4/CaseStudyII/input_files/Grain-"+std::to_string(grainNo)+".inp";
		std::ifstream input(filename);
		assert(input.is_open());
		count = 0;
		while(std::getline(input,line)){
			count++;
		}
		input.close();
		input.open(filename);
		for(unsigned int i=0;i<(count-1);i++){
			input >> voronoi_edge_x;
			input >> voronoi_edge_y;
			x_edge[grainNo-1].push_back(voronoi_edge_x);
			y_edge[grainNo-1].push_back(voronoi_edge_y);
			input.seekg(1,input.cur);
		}
		input.close();
		grainNo ++;

	}

	return 0;
}



