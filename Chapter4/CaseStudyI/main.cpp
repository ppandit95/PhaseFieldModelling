/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////             CASE STUDY I : Simulation of the Spinodal Decomposition of a binary alloy                   ////////
////////                with explicit Euler Finite Difference Algorithm                                          ////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<time.h>
#include<fstream>
#include<cmath>
#include<string>
#include<filesystem>
#include <unistd.h>

#include"src/structures.hpp"
#include"src/SpinodalDecomposition.hpp"


int main(int argc, char **argv) {
	//Get initial wall time
	struct timespec begin,end;
	clock_gettime(CLOCK_REALTIME,&begin);
	SpinodalDecomposition sd;


	//Simulation cell parameters
	SimParams p1 = sd.GetSimParams();


	//Write energy functional to dat file
	std::string filename0 = getcwd(p1.path,256) + (std::string)"/Output/Gibbs_Energy.dat";
	std::ofstream energ(filename0);


	//Time Integration Parameters
	TimeParams t1 = sd.GetTimeParams();


	//Material Specific Parameters
	MatParams m1 = sd.GetMatParams();

	//Initial Concentration Profile
	sd.Initial_Profile();

	std::string filename = getcwd(p1.path,256) + (std::string)"/Output/Initial_Profile.dat";
	//Write the initial Profile in a dat file
	sd.write_output(filename);

	//Evolve with Cahn Hilliard Formulation
	for(unsigned int istep=0;istep<t1.nstep;istep++){
		m1.energy = 0.0;
		t1.ttime += t1.dtime;
		for(int i=0;i<p1.Nx;i++){
			for(int j=0;j<p1.Ny;j++){
				//Calculate bulk energy
				if(i<(p1.Nx-1) && j<(p1.Ny-1))
					m1.energy += p1.con[i][j]*p1.con[i][j]*(1-p1.con[i][j])*(1-p1.con[i][j])
							+ (m1.grad_coef/2.0)*((p1.con[i+1][j]-p1.con[i][j])*(p1.con[i+1][j]-p1.con[i][j]) + (p1.con[i][j+1]-p1.con[i][j])*(p1.con[i][j+1]-p1.con[i][j]));
				else
					m1.energy += p1.con[i][j]*p1.con[i][j]*(1-p1.con[i][j])*(1-p1.con[i][j]);

				//Compute the chemical potential at each time step and each point
				p1.mu[i][j] = m1.A*(2.0*p1.con[i][j]*std::pow((1-p1.con[i][j]),2) - 2.0*std::pow(p1.con[i][j],2)*(1-p1.con[i][j]));

				//Compute laplacian of concentration in case of periodic BCs
				if(i==0&&j==0)
					p1.laplace_con[i][j] = std::pow(1/p1.dx,2)*(p1.con[i+1][j]-2*p1.con[i][j]+p1.con[i+p1.Nx-1][j])
										+ std::pow(1/p1.dy,2)*(p1.con[i][j+1]-2*p1.con[i][j]+p1.con[i][j+p1.Ny-1]);
				else if(i==(p1.Nx-1)&&j==0)
					p1.laplace_con[i][j] = std::pow(1/p1.dx,2)*(p1.con[i+1-p1.Nx][j]-2*p1.con[i][j]+p1.con[i-1][j])
										+ std::pow(1/p1.dy,2)*(p1.con[i][j+1]-2*p1.con[i][j]+p1.con[i][j+p1.Ny-1]);
				else if(i==(p1.Nx-1)&&j==(p1.Ny-1))
					p1.laplace_con[i][j] = std::pow(1/p1.dx,2)*(p1.con[i+1-p1.Nx][j]-2*p1.con[i][j]+p1.con[i-1][j])
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

				//Compute functional derivative of free energy with concentration
				p1.dfdc[i][j] = p1.mu[i][j] - m1.grad_coef*p1.laplace_con[i][j];

				//Compute laplacian of functional derivative
				if(i==0&&j==0)
					p1.laplace_dfdc[i][j] = std::pow(1/p1.dx,2)*(p1.dfdc[i+1][j]-2*p1.dfdc[i][j]+p1.dfdc[i+p1.Nx-1][j])
										+ std::pow(1/p1.dy,2)*(p1.dfdc[i][j+1]-2*p1.dfdc[i][j]+p1.dfdc[i][j+p1.Ny-1]);
				else if(i==(p1.Nx-1)&&j==0)
					p1.laplace_dfdc[i][j] = std::pow(1/p1.dx,2)*(p1.dfdc[i+1-p1.Nx][j]-2*p1.dfdc[i][j]+p1.dfdc[i-1][j])
														+ std::pow(1/p1.dy,2)*(p1.dfdc[i][j+1]-2*p1.dfdc[i][j]+p1.dfdc[i][j+p1.Ny-1]);
				else if(i==(p1.Nx-1)&&j==(p1.Ny-1))
					p1.laplace_dfdc[i][j] = std::pow(1/p1.dx,2)*(p1.dfdc[i+1-p1.Nx][j]-2*p1.dfdc[i][j]+p1.dfdc[i-1][j])
														+ std::pow(1/p1.dy,2)*(p1.dfdc[i][j+1-p1.Ny]-2*p1.dfdc[i][j]+p1.dfdc[i][j-1]);
				else if(i==0&&j==(p1.Ny-1))
					p1.laplace_dfdc[i][j] = std::pow(1/p1.dx,2)*(p1.dfdc[i+1][j]-2*p1.dfdc[i][j]+p1.dfdc[i-1+p1.Nx][j])
														+ std::pow(1/p1.dy,2)*(p1.dfdc[i][j+1-p1.Ny]-2*p1.dfdc[i][j]+p1.dfdc[i][j-1]);
				else if(i==0)
					p1.laplace_dfdc[i][j] = std::pow(1/p1.dx,2)*(p1.dfdc[i+1][j]-2*p1.dfdc[i][j]+p1.dfdc[i+p1.Nx-1][j]) +
														std::pow(1/p1.dy,2)*(p1.dfdc[i][j+1]-2*p1.dfdc[i][j]+p1.dfdc[i][j-1]);
				else if(i==(p1.Nx-1))
					p1.laplace_dfdc[i][j] = std::pow(1/p1.dx,2)*(p1.dfdc[i+1-p1.Nx][j]-2*p1.dfdc[i][j]+p1.dfdc[i-1][j]) +
														std::pow(1/p1.dy,2)*(p1.dfdc[i][j+1]-2*p1.dfdc[i][j]+p1.dfdc[i][j-1]);
				else if(j==0)
					p1.laplace_dfdc[i][j] = std::pow(1/p1.dx,2)*(p1.dfdc[i+1][j]-2*p1.dfdc[i][j]+p1.dfdc[i-1][j]) +
														std::pow(1/p1.dy,2)*(p1.dfdc[i][j+1]-2*p1.dfdc[i][j]+p1.dfdc[i][j-1+p1.Ny]);
				else if(j==(p1.Ny-1))
					p1.laplace_dfdc[i][j] = std::pow(1/p1.dx,2)*(p1.dfdc[i+1][j]-2*p1.dfdc[i][j]+p1.dfdc[i-1][j]) +
														std::pow(1/p1.dy,2)*(p1.dfdc[i][j+1-p1.Ny]-2*p1.dfdc[i][j]+p1.dfdc[i][j-1]);
				else
					p1.laplace_dfdc[i][j] = std::pow(1/p1.dx,2)*(p1.dfdc[i+1][j]-2*p1.dfdc[i][j]+p1.dfdc[i-1][j]) +
														std::pow(1/p1.dy,2)*(p1.dfdc[i][j+1]-2*p1.dfdc[i][j]+p1.dfdc[i][j-1]);

				//Evolve the concentration profile

				p1.con_temp[i][j] = p1.con[i][j] + t1.dtime*m1.mobility*p1.laplace_dfdc[i][j];
			}
		}
		//Backtransfering from temporary array to con array
		for(int i=0;i<p1.Nx;i++){
			for(int j=0;j<p1.Ny;j++)
				p1.con[i][j] = p1.con_temp[i][j];
		}
		energ<<m1.energy<<std::endl;

		//Write the Chemical Potential in a dat file at specific intervels
		if(istep%1000 == 0){
			std::string filename = getcwd(p1.path,256) + (std::string)"/Output/Chemical_Potential_at_t_"+std::to_string(istep)+(std::string)".dat";
			std::ofstream MuProf(filename);
			for(int i=0;i<p1.Nx;i++){
				for(int j=0;j<p1.Ny;j++)
					MuProf << i <<"\t"<<j<<"\t"<<p1.mu[i][j] << "\n";
			}
			MuProf.close();
		}
		//Write the Derivative in a dat file at specific intervels
		if(istep%1000 == 0){
		std::string filename = getcwd(p1.path,256) +(std::string)"/Output/Derivative_at_t_"+std::to_string(istep)+(std::string)".dat";
		std::ofstream DerProf(filename);
		for(int i=0;i<p1.Nx;i++){
				for(int j=0;j<p1.Ny;j++)
					DerProf << i <<"\t"<<j<<"\t"<<p1.dfdc[i][j] << "\n";
			}
			DerProf.close();
		}
		//Write the concentration profile at specific intervel
		if(istep%1000==0){
			std::string filename = getcwd(p1.path,256) + (std::string)"/Output/Concentration_Profile_t_" + std::to_string(istep) + (std::string)".dat";
			std::ofstream ConcProf(filename);
			for(int i=0;i<p1.Nx;i++){
				for(int j=0;j<p1.Ny;j++)
					ConcProf << i<<"\t"<<j<<"\t"<<p1.con[i][j]<<std::endl;
			}
			ConcProf.close();
		}
	}
	energ.close();
	for(int i=0;i<p1.Nx;i++){
		delete[] p1.con[i];
		delete[] p1.mu[i];
		delete[] p1.dfdc[i];
		delete[] p1.laplace_con[i];
		delete[] p1.laplace_dfdc[i];
		delete[] p1.con_temp[i];

	}

	delete[] p1.con_temp;
	delete[] p1.laplace_dfdc;
	delete[] p1.laplace_con;
	delete[] p1.dfdc;
	delete[] p1.mu;
	delete[] p1.con;

	clock_gettime(CLOCK_REALTIME,&end);
	std::cout<<"Elapsed Time : "<<(end.tv_sec-begin.tv_sec)<<" seconds "<<(end.tv_nsec - begin.tv_nsec)<<" nanoseconds"<<std::endl;
	return 0;
}
