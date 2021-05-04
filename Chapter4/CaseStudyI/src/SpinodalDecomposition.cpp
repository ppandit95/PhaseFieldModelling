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
#include<string>
#include<filesystem>
#include <unistd.h>

SpinodalDecomposition::SpinodalDecomposition(){
	p.Nx = 64;
	p.Ny = 64;
	p.dx = 1.0;
	p.dy = 1.0;

	t.nstep = 10000;
	t.dtime = 0.01;
	t.ttime = 0.0;

	m.c0 = 0.4;
	m.mobility = 1.0;
	m.grad_coef = 0.5;
	m.A = 1.0;


	p.con = new double*[p.Nx];
	p.mu = new double* [p.Nx];
	p.dfdc = new double* [p.Nx];
	p.laplace_con = new double* [p.Nx];
	p.laplace_dfdc = new double* [p.Nx];
	p.con_temp = new double* [p.Nx];
	for(unsigned int i=0;i<p.Nx;i++){
		p.con[i] = new double[p.Ny];
		p.mu[i] = new double[p.Ny];
		p.dfdc[i] = new double[p.Ny];
		p.laplace_con[i] = new double[p.Ny];
		p.laplace_dfdc[i] = new double[p.Ny];
		p.con_temp[i] = new double[p.Ny];
	}
}

SpinodalDecomposition::~SpinodalDecomposition(){
	for(unsigned int i=0;i<p.Nx;i++){
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
	for(unsigned int i=0;i<p.Nx;i++){
		for(unsigned int j=0;j<p.Ny;j++){
			p.con[i][j] = m.c0 + noise*(0.5 - (float) rand()/RAND_MAX);
		}
	}
}
void SpinodalDecomposition::write_output(std::string filename,double** con){
	std::ofstream output(filename);
	for(unsigned int i=0;i<p.Nx;i++){
		for(unsigned int j=0;j<p.Ny;j++)
			output << i <<"\t"<<j<<"\t"<<con[i][j] << "\n";
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

	if(count==0){
		if(i==0&&j==0)
			p.laplace_con[i][j] = std::pow(1/p.dx,2)*(p.con[i+1][j]-2*p.con[i][j]+p.con[i+p.Nx-1][j])
								+ std::pow(1/p.dy,2)*(p.con[i][j+1]-2*p.con[i][j]+p.con[i][j+p.Ny-1]);
		else if(i==(p.Nx-1)&&j==0)
			p.laplace_con[i][j] = std::pow(1/p.dx,2)*(p.con[i+1-p.Nx][j]-2*p.con[i][j]+p.con[i-1][j])
								+ std::pow(1/p.dy,2)*(p.con[i][j+1]-2*p.con[i][j]+p.con[i][j+p.Ny-1]);
		else if(i==(p.Nx-1)&&j==(p.Ny-1))
			p.laplace_con[i][j] = std::pow(1/p.dx,2)*(p.con[i+1-p.Nx][j]-2*p.con[i][j]+p.con[i-1][j])
								+ std::pow(1/p.dy,2)*(p.con[i][j+1-p.Ny]-2*p.con[i][j]+p.con[i][j-1]);
		else if(i==0&&j==(p.Ny-1))
			p.laplace_con[i][j] = std::pow(1/p.dx,2)*(p.con[i+1][j]-2*p.con[i][j]+p.con[i-1+p.Nx][j])
								+ std::pow(1/p.dy,2)*(p.con[i][j+1-p.Ny]-2*p.con[i][j]+p.con[i][j-1]);
		else if(i==0)
			p.laplace_con[i][j] = std::pow(1/p.dx,2)*(p.con[i+1][j]-2*p.con[i][j]+p.con[i+p.Nx-1][j]) +
								  std::pow(1/p.dy,2)*(p.con[i][j+1]-2*p.con[i][j]+p.con[i][j-1]);
		else if(i==(p.Nx-1))
			p.laplace_con[i][j] = std::pow(1/p.dx,2)*(p.con[i+1-p.Nx][j]-2*p.con[i][j]+p.con[i-1][j]) +
								  std::pow(1/p.dy,2)*(p.con[i][j+1]-2*p.con[i][j]+p.con[i][j-1]);
		else if(j==0)
			p.laplace_con[i][j] = std::pow(1/p.dx,2)*(p.con[i+1][j]-2*p.con[i][j]+p.con[i-1][j]) +
								  std::pow(1/p.dy,2)*(p.con[i][j+1]-2*p.con[i][j]+p.con[i][j-1+p.Ny]);
		else if(j==(p.Ny-1))
			p.laplace_con[i][j] = std::pow(1/p.dx,2)*(p.con[i+1][j]-2*p.con[i][j]+p.con[i-1][j]) +
								  std::pow(1/p.dy,2)*(p.con[i][j+1-p.Ny]-2*p.con[i][j]+p.con[i][j-1]);
		else
			p.laplace_con[i][j] = std::pow(1/p.dx,2)*(p.con[i+1][j]-2*p.con[i][j]+p.con[i-1][j]) +
												std::pow(1/p.dy,2)*(p.con[i][j+1]-2*p.con[i][j]+p.con[i][j-1]);
	}
	else if(count==1){
		if(i==0&&j==0)
			p.laplace_dfdc[i][j] = std::pow(1/p.dx,2)*(p.dfdc[i+1][j]-2*p.dfdc[i][j]+p.dfdc[i+p.Nx-1][j])
									+ std::pow(1/p.dy,2)*(p.dfdc[i][j+1]-2*p.dfdc[i][j]+p.dfdc[i][j+p.Ny-1]);
		else if(i==(p.Nx-1)&&j==0)
			p.laplace_dfdc[i][j] = std::pow(1/p.dx,2)*(p.dfdc[i+1-p.Nx][j]-2*p.dfdc[i][j]+p.dfdc[i-1][j])
									+ std::pow(1/p.dy,2)*(p.dfdc[i][j+1]-2*p.dfdc[i][j]+p.dfdc[i][j+p.Ny-1]);
		else if(i==(p.Nx-1)&&j==(p.Ny-1))
			p.laplace_dfdc[i][j] = std::pow(1/p.dx,2)*(p.dfdc[i+1-p.Nx][j]-2*p.dfdc[i][j]+p.dfdc[i-1][j])
								+ std::pow(1/p.dy,2)*(p.dfdc[i][j+1-p.Ny]-2*p.dfdc[i][j]+p.dfdc[i][j-1]);
		else if(i==0&&j==(p.Ny-1))
			p.laplace_dfdc[i][j] = std::pow(1/p.dx,2)*(p.dfdc[i+1][j]-2*p.dfdc[i][j]+p.dfdc[i-1+p.Nx][j])
								+ std::pow(1/p.dy,2)*(p.dfdc[i][j+1-p.Ny]-2*p.dfdc[i][j]+p.dfdc[i][j-1]);
		else if(i==0)
			p.laplace_dfdc[i][j] = std::pow(1/p.dx,2)*(p.dfdc[i+1][j]-2*p.dfdc[i][j]+p.dfdc[i+p.Nx-1][j]) +
									std::pow(1/p.dy,2)*(p.dfdc[i][j+1]-2*p.dfdc[i][j]+p.dfdc[i][j-1]);
		else if(i==(p.Nx-1))
			p.laplace_dfdc[i][j] = std::pow(1/p.dx,2)*(p.dfdc[i+1-p.Nx][j]-2*p.dfdc[i][j]+p.dfdc[i-1][j]) +
									std::pow(1/p.dy,2)*(p.dfdc[i][j+1]-2*p.dfdc[i][j]+p.dfdc[i][j-1]);
		else if(j==0)
			p.laplace_dfdc[i][j] = std::pow(1/p.dx,2)*(p.dfdc[i+1][j]-2*p.dfdc[i][j]+p.dfdc[i-1][j]) +
									std::pow(1/p.dy,2)*(p.dfdc[i][j+1]-2*p.dfdc[i][j]+p.dfdc[i][j-1+p.Ny]);
		else if(j==(p.Ny-1))
			p.laplace_dfdc[i][j] = std::pow(1/p.dx,2)*(p.dfdc[i+1][j]-2*p.dfdc[i][j]+p.dfdc[i-1][j]) +
								   std::pow(1/p.dy,2)*(p.dfdc[i][j+1-p.Ny]-2*p.dfdc[i][j]+p.dfdc[i][j-1]);
		else
			p.laplace_dfdc[i][j] = std::pow(1/p.dx,2)*(p.dfdc[i+1][j]-2*p.dfdc[i][j]+p.dfdc[i-1][j]) +
								   std::pow(1/p.dy,2)*(p.dfdc[i][j+1]-2*p.dfdc[i][j]+p.dfdc[i][j-1]);
	}
}
void SpinodalDecomposition::CalcEnergyDer(unsigned int i,unsigned int j){
	p.dfdc[i][j] = p.mu[i][j] - m.grad_coef*p.laplace_con[i][j];
}
void SpinodalDecomposition::EvolveConc(unsigned int i,unsigned int j){
	p.con_temp[i][j] = p.con[i][j] + t.dtime*m.mobility*p.laplace_dfdc[i][j];
}
void SpinodalDecomposition::BackTransfer(){
	for(unsigned int i=0;i<p.Nx;i++){
		for(unsigned int j=0;j<p.Ny;j++)
			p.con[i][j] = p.con_temp[i][j];
	}
}
void SpinodalDecomposition::solve(){
	std::string filename0 = getcwd(p.path,256) + (std::string)"/Output/Gibbs_Energy.dat";
	std::ofstream energ(filename0);
	//Evolve with Cahn Hilliard Formulation
	for(unsigned int istep=0;istep<t.nstep;istep++){
			m.energy = 0.0;
			t.ttime += t.dtime;
			for(unsigned int i=0;i<p.Nx;i++){
				for(unsigned int j=0;j<p.Ny;j++){
					//Calculate bulk energy
					CalcBulkEnergy(i, j);

					//Compute the chemical potential at each time step and each point
					CalcChemicalPotential(i, j);
					//Compute laplacian of concentration in case of periodic BCs
					Set_Periodic_BCs(i, j, 0);

					//Compute functional derivative of free energy with concentration
					CalcEnergyDer(i, j);

					//Compute laplacian of functional derivative
					Set_Periodic_BCs(i, j, 1);

					//Evolve the concentration profile
					EvolveConc(i, j);
				}
			}
			//Backtransfering from temporary array to con array
			BackTransfer();
			energ<<m.energy<<std::endl;

			//Write the Chemical Potential in a dat file at specific intervels
			if(istep%1000 == 0){
				std::string filename = getcwd(p.path,256) + (std::string)"/Output/Chemical_Potential_at_t_"+std::to_string(istep)+(std::string)".dat";
				write_output(filename,p.mu);
			}
			//Write the Derivative in a dat file at specific intervels
			if(istep%1000 == 0){
			std::string filename = getcwd(p.path,256) +(std::string)"/Output/Derivative_at_t_"+std::to_string(istep)+(std::string)".dat";
			write_output(filename,p.dfdc);
			}
			//Write the concentration profile at specific intervel
			if(istep%1000==0){
				std::string filename = getcwd(p.path,256) + (std::string)"/Output/Concentration_Profile_t_" + std::to_string(istep) + (std::string)".dat";
				write_output(filename, p.con);
			}
		}
		energ.close();
}
void SpinodalDecomposition::write_initial_profile(){
	std::string filename = getcwd(p.path,256) + (std::string)"/Output/Initial_Profile.dat";
	write_output(filename,p.con);
}


