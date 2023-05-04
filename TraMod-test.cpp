/**
	Test for the position trajectory modification for position, velocity and acceleration constraints.
	Note that the algorithm needs Eigen (https://eigen.tuxfamily.org/index.php?title=Main_Page).
	@author: Hisayoshi Muramatsu
	@date: 2023.05.02
*/

#include <iostream>
#include <fstream>
#include <math.h>

#define DtoInt6(x) ((int)round(1000000*x))

#include "TraMod.cpp"

int main(){
	std::cout << "=====================================" << std::endl;
	std::cout << "       Trajectory Modification       " << std::endl;
	std::cout << "=====================================" << std::endl;

	const double tend=10; // Test time [s]
	const double T=0.0005; // Test time [s]

	// No position constraint mode if x_max=x_min=0
	const double x_max=1; // Maximum position constraint
	const double x_min=-2; // Minimum position constraint
	const double dx_max=1.5; // Maximum velocity constraint
	const double dx_min=-1.5; // Minimum velocity constraint
	const double ddx_max=10; // Maximum acceleration constraint
	const double ddx_min=-10; // Minimum acceleration constraint

	const double w1(700), w2(700); // Weights for QP
	const double epth(1e-8), kth(10); // Thresholds for QP

	// Variables
	double t=0;
	double x_k=-0.5, x_k1=0, dx_k=0, dx_k1=0, ddx_k=0;
	double MOD_x_k=0, MOD_x_k1=0, MOD_dx_k=0, MOD_dx_k1=0, MOD_ddx_k=0;
	int iflag=0, pflag=0;

	TraMod* TM;
	TM = new TraMod(x_k, w1, w2, x_max, x_min, dx_max, dx_min, ddx_max, ddx_min, epth, kth, T);

	std::ofstream ofs("DATA.dat");

	while((int)round(t/T)<=(int)round(tend/T)){
		/**
			Input trajectory
		*/
		if(t<0.5){
			x_k=-0.5;
		}else if(t<2.5){
			x_k=1.5;
		}else if(t<5){
			x_k=0.5;
		}else{
			x_k=0.5*sin(2*t);
		}

		/**
			Position trajectory modification
		*/
		MOD_x_k = TM->TrajectoryModification(x_k);

		/**
			Computation of velocity and acceleration
		*/
		if(iflag==0){ // Initial operation
			x_k1=x_k;
			MOD_x_k1=MOD_x_k;
			if(DtoInt6(x_min)==DtoInt6(x_min)&&DtoInt6(x_min)==0) pflag=1;
			iflag++;
		}
		dx_k =(x_k-x_k1)/T;
			x_k1=x_k;
		ddx_k=(dx_k-dx_k1)/T;
			dx_k1=dx_k;
		MOD_dx_k =(MOD_x_k-MOD_x_k1)/T;
			MOD_x_k1=MOD_x_k;
		MOD_ddx_k=(MOD_dx_k-MOD_dx_k1)/T;
			MOD_dx_k1=MOD_dx_k;

		/**
			Check of modified trajectory
		*/
		if(DtoInt6(MOD_x_k)<DtoInt6(x_min) || DtoInt6(x_max)<DtoInt6(MOD_x_k)){
			if(pflag!=1) std::cout << std::endl << t << "s error" << ", reference position failed to satisfy the constraints as" << x_min << " m <=" << MOD_x_k << " m <= " << x_max << " m" << std::endl;
		}
		if(DtoInt6(MOD_dx_k)<DtoInt6(dx_min) || DtoInt6(dx_max)<DtoInt6(MOD_dx_k)){
			std::cout << std::endl << t << "s error" << ", reference velocity failed to satisfy the constraints as" << dx_min << " m/s <=" << MOD_dx_k << " m/s <= " << dx_max << " m/s" << std::endl;
		}
		if(DtoInt6(MOD_ddx_k)<DtoInt6(ddx_min) || DtoInt6(ddx_max)<DtoInt6(MOD_ddx_k)){
			std::cout << std::endl << t << "s error" << ", reference acceleration failed to satisfy the constraints as" << ddx_min << " m/s^2 <=" << MOD_ddx_k << " m/s^2 <= " << ddx_max << " m/s^2" << std::endl;
		}

		/**
			Output DATA
		*/
		ofs << t << " " << x_k << " " << MOD_x_k << " " << dx_k << " " << MOD_dx_k << " " << ddx_k << " " << MOD_ddx_k << std::endl;
		if((int)round(t/T)%(int)round(0.5/T)==0){
			std::cout << std::setprecision(3) << "Time: " << t << " s, Pos.: " << MOD_x_k << " m, Vel.: " << MOD_dx_k << " m/s, Acc.: " << MOD_x_k << " m/s^2" << std::endl;
		}

		t+=T;

	}

	ofs.close();

	std::cout << "=============== Finish ==============" << std::endl;

	return 0;

}
