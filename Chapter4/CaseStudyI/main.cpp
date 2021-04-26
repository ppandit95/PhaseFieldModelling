/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////             CASE STUDY I : Simulation of the Spinodal Decomposition of a binary alloy                   ////////
////////                with explicit Euler Finite Difference Algorithm                                          ////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include<iostream>
#include<time.h>
#include<fstream>
#include<cmath>
int main(int argc, char **argv) {
	//Get initial wall time
	struct timespec begin,end;
	clock_gettime(CLOCK_REALTIME,&begin);

	//Open output file to write energy functional values at each timestep
	//std::ofstream out("time_energ.dat");

	//Simulation cell parameters
	int Nx = 64;
	int Ny = 64;
	double dx = 1.0;
	double dy = 1.0;
	double** con,**mu,**dfdc,**laplace_con;

	//Time Integration Parameters
	int nstep = 10000;
	int nprint = 50;
	double dtime = 0.01;
	double ttime = 0.0;


	//Material Specific Parameters
	double c0 = 0.4;
	double mobility = 1.0;
	double grad_coef = 0.5;
	double A = 1.0;

	//Initial concentration profile
	srand( (unsigned)time( NULL ) );
	con = new double*[Nx];
	mu = new double* [Nx];
	dfdc = new double* [Nx];
	laplace_con = new double* [Nx];
	for(int i=0;i<Nx;i++){
		con[i] = new double[Ny];
		mu[i] = new double[Ny];
		dfdc[i] = new double[Ny];
		laplace_con[i] = new double[Ny];
	}
	double noise = 0.02;
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			con[i][j] = c0 + noise*(0.5 - (float) rand()/RAND_MAX);
		}
	}

	//Write the initial Profile in a dat file
	std::ofstream initProf("Initial_Profile.dat");
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++)
			initProf << i <<"\t"<<j<<"\t"<<con[i][j] << "\n";
		//initProf << std::endl;
	}

	//Evolve with Cahn Hilliard Formulation
	for(int istep=0;istep<nstep;istep++){
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				//Compute the chemical potential at each time step and each point
				mu[i][j] = A*(2.0*con[i][j]*con[i][j]*std::pow((1-con[i][j]),2) + 2.0*std::pow(con[i][j],2)*(1-con[i][j]));

				//Compute laplacian of concentration in case of periodic BCs
				if(i==0&&j==0)
					laplace_con[i][j] = std::pow(1/dx,2)*(con[i+1][j]-2*con[i][j]+con[i+Nx-1][j])
										+ std::pow(1/dy,2)*(con[i][j+1]-2*con[i][j]+con[i][j+Ny-1]);
				else if(i==(Nx-1)&&j==0)
					laplace_con[i][j] = std::pow(1/dx,2)*(con[i+1-Nx][j]-2*con[i][j]+con[i-1][j])
										+ std::pow(1/dy,2)*(con[i][j+1]-2*con[i][j]+con[i][j+Ny-1]);
				else if(i==(Nx-1)&&j==(Ny-1))
					laplace_con[i][j] = std::pow(1/dx,2)*(con[i+1-Nx][j]-2*con[i][j]+con[i-1][j])
										+ std::pow(1/dy,2)*(con[i][j+1-Ny]-2*con[i][j]+con[i][j-1]);
				else if(i==0&&j==(Ny-1))
					laplace_con[i][j] = std::pow(1/dx,2)*(con[i+1][j]-2*con[i][j]+con[i-1+Nx][j])
															+ std::pow(1/dy,2)*(con[i][j+1-Ny]-2*con[i][j]+con[i][j-1]);
				else if(i==0)
					laplace_con[i][j] = std::pow(1/dx,2)*(con[i+1][j]-2*con[i][j]+con[i+Nx-1][j]) +
										std::pow(1/dy,2)*(con[i][j+1]-2*con[i][j]+con[i][j-1]);
				else if(i==(Nx-1))
					laplace_con[i][j] = std::pow(1/dx,2)*(con[i+1-Nx][j]-2*con[i][j]+con[i-1][j]) +
										std::pow(1/dy,2)*(con[i][j+1]-2*con[i][j]+con[i][j-1]);
				else if(j==0)
					laplace_con[i][j] = std::pow(1/dx,2)*(con[i+1][j]-2*con[i][j]+con[i-1][j]) +
										std::pow(1/dy,2)*(con[i][j+1]-2*con[i][j]+con[i][j-1+Ny]);
				else if(j==(Ny-1))
					laplace_con[i][j] = std::pow(1/dx,2)*(con[i+1][j]-2*con[i][j]+con[i-1][j]) +
										std::pow(1/dy,2)*(con[i][j+1-Ny]-2*con[i][j]+con[i][j-1]);
				else
					laplace_con[i][j] = std::pow(1/dx,2)*(con[i+1][j]-2*con[i][j]+con[i-1][j]) +
										std::pow(1/dy,2)*(con[i][j+1]-2*con[i][j]+con[i][j-1]);

				//Compute functional derivative of free energy with concentration
				dfdc[i][j] = mu[i][j] - grad_coef*laplace_con[i][j];
			}
		}
		//Write the Chemical Potential in a dat file at specific intervels
		if(istep%1000 == 0){
			std::string filename = "Chemical_Potential_at_t_"+std::to_string(istep)+".dat";
			std::ofstream MuProf(filename);
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++)
					MuProf << i <<"\t"<<j<<"\t"<<con[i][j] << "\n";
			}
			MuProf.close();
		}
		//Write the Derivative in a dat file at specific intervels
		if(istep%1000 == 0){
		std::string filename = "Derivative_at_t_"+std::to_string(istep)+".dat";
		std::ofstream DerProf(filename);
		for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++)
					DerProf << i <<"\t"<<j<<"\t"<<con[i][j] << "\n";
			}
			DerProf.close();
		}

	}
	for(int i=0;i<Nx;i++){
		delete[] con[i];
		delete[] mu[i];
		delete[] dfdc[i];
		delete[] laplace_con[i];
	}
	delete[] laplace_con;
	delete[] dfdc;
	delete[] mu;
	delete[] con;
	return 0;
}
