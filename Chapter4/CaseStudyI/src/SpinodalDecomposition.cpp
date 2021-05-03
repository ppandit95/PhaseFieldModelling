/*
 * SpinodalDecomposition.cpp
 *
 *  Created on: 03-May-2021
 *      Author: pushkar
 */
#include"SpinodalDecomposition.hpp"
#include"structures.hpp"
#include<iostream>
#include<fstream>
#include<cmath>

SpinodalDecomposition::SpinodalDecomposition(){
	p.Nx = 64;
	p.Ny = 64;
	p.dx = 1.0;
	p.dy = 1.0;

	t.nstep = 10000;
	t.dtime = 0.01;
	t.ttime = 0.0;

	double c0 = 0.4;
	double mobility = 1.0;
	double grad_coef = 0.5;
	double A = 1.0;
	double energy;

	p.con = new double*[p.Nx];
	p.mu = new double* [p.Nx];
	p.dfdc = new double* [p.Nx];
	p.laplace_con = new double* [p.Nx];
	p.laplace_dfdc = new double* [p.Nx];
	p.con_temp = new double* [p.Nx];
	for(int i=0;i<p.Nx;i++){
		p.con[i] = new double[p.Ny];
		p.mu[i] = new double[p.Ny];
		p.dfdc[i] = new double[p.Ny];
		p.laplace_con[i] = new double[p.Ny];
		p.laplace_dfdc[i] = new double[p.Ny];
		p.con_temp[i] = new double[p.Ny];
	}
}

SpinodalDecomposition::~SpinodalDecomposition(){
	for(int i=0;i<p.Nx;i++){
		delete[] p.con[i];
		delete[] p.mu[i];
		delete[] p.dfdc[i];
		delete[] p.laplace_con[i];
		delete[] p.laplace_dfdc[i];
		delete[] p.con_temp[i];

	}

	delete[] p.con_temp;
	delete[] p.laplace_dfdc;
	delete[] p.laplace_con;
	delete[] p.dfdc;
	delete[] p.mu;
	delete[] p.con;
}
SimParams SpinodalDecomposition::GetSimParams(){
	return p;
}
TimeParams SpinodalDecomposition::GetTimeParams(){
	return t;
}
MatParams SpinodalDecomposition::GetMatParams(){
	return m;
}
void SpinodalDecomposition::Initial_Profile(){
	srand( (unsigned)time( NULL ) );
	double noise = 0.02;
	for(int i=0;i<p.Nx;i++){
		for(int j=0;j<p.Ny;j++){
			p.con[i][j] = m.c0 + noise*(0.5 - (float) rand()/RAND_MAX);
		}
	}
}
void SpinodalDecomposition::write_output(std::string filename){
	std::ofstream output(filename);
	for(int i=0;i<p.Nx;i++){
		for(int j=0;j<p.Ny;j++)
			output << i <<"\t"<<j<<"\t"<<p.con[i][j] << "\n";
		}
	output.close();
}
void SpinodalDecomposition::CalcBulkEnergy(unsigned int i,unsigned int j){
	if(i<(p.Nx-1) && j<(p.Ny-1))
		m.energy += p.con[i][j]*p.con[i][j]*(1-p.con[i][j])*(1-p.con[i][j])
					+ (m.grad_coef/2.0)*((p.con[i+1][j]-p.con[i][j])*(p.con[i+1][j]-p.con[i][j]) + (p.con[i][j+1]-p.con[i][j])*(p.con[i][j+1]-p.con[i][j]));
	else
		m.energy += p.con[i][j]*p.con[i][j]*(1-p.con[i][j])*(1-p.con[i][j]);

}
void SpinodalDecomposition::CalcChemicalPotential(unsigned int i,unsigned int j){
	p.mu[i][j] = m.A*(2.0*p.con[i][j]*std::pow((1-p.con[i][j]),2) - 2.0*std::pow(p.con[i][j],2)*(1-p.con[i][j]));
}
void SpinodalDecomposition::Set_Periodic_BCs(unsigned int i,unsigned int j,unsigned int count){
	if(i==0&&j==0)
		laplace_con[i][j] = std::pow(1/p.dx,2)*(con[i+1][j]-2*con[i][j]+con[i+p.Nx-1][j])
							+ std::pow(1/p.dy,2)*(con[i][j+1]-2*con[i][j]+con[i][j+p.Ny-1]);
	else if(i==(p.Nx-1)&&j==0)
		laplace_con[i][j] = std::pow(1/p.dx,2)*(con[i+1-p.Nx][j]-2*con[i][j]+con[i-1][j])
							+ std::pow(1/p.dy,2)*(con[i][j+1]-2*con[i][j]+con[i][j+p.Ny-1]);
	else if(i==(p.Nx-1)&&j==(p.Ny-1))
		laplace_con[i][j] = std::pow(1/p.dx,2)*(con[i+1-p1.Nx][j]-2*con[i][j]+p1.con[i-1][j])
											+ std::pow(1/p1.dy,2)*(p1.con[i][j+1-p1.Ny]-2*p1.con[i][j]+p1.con[i][j-1]);
					else if(i==0&&j==(p1.Ny-1))
						p1.laplace_con[i][j] = std::pow(1/p1.dx,2)*(p1.con[i+1][j]-2*p1.con[i][j]+p1.con[i-1+p1.Nx][j])
																+ std::pow(1/p1.dy,2)*(p1.con[i][j+1-p1.Ny]-2*p1.con[i][j]+p1.con[i][j-1]);
					else if(i==0)
						p1.laplace_con[i][j] = std::pow(1/p1.dx,2)*(p1.con[i+1][j]-2*p1.con[i][j]+p1.con[i+p1.Nx-1][j]) +
											std::pow(1/p1.dy,2)*(p1.con[i][j+1]-2*p1.con[i][j]+p1.con[i][j-1]);
					else if(i==(p1.Nx-1))
						p1.laplace_con[i][j] = std::pow(1/p1.dx,2)*(p1.con[i+1-p1.Nx][j]-2*p1.con[i][j]+p1.con[i-1][j]) +
											std::pow(1/p1.dy,2)*(p1.con[i][j+1]-2*p1.con[i][j]+p1.con[i][j-1]);
					else if(j==0)
						p1.laplace_con[i][j] = std::pow(1/p1.dx,2)*(p1.con[i+1][j]-2*p1.con[i][j]+p1.con[i-1][j]) +
											std::pow(1/p1.dy,2)*(p1.con[i][j+1]-2*p1.con[i][j]+p1.con[i][j-1+p1.Ny]);
					else if(j==(p1.Ny-1))
						p1.laplace_con[i][j] = std::pow(1/p1.dx,2)*(p1.con[i+1][j]-2*p1.con[i][j]+p1.con[i-1][j]) +
											std::pow(1/p1.dy,2)*(p1.con[i][j+1-p1.Ny]-2*p1.con[i][j]+p1.con[i][j-1]);
					else
						p1.laplace_con[i][j] = std::pow(1/p1.dx,2)*(p1.con[i+1][j]-2*p1.con[i][j]+p1.con[i-1][j]) +
											std::pow(1/p1.dy,2)*(p1.con[i][j+1]-2*p1.con[i][j]+p1.con[i][j-1]);

}


