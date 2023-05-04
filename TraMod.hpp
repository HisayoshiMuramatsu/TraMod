/**
	Position trajectory modification for position, velocity and acceleration constraints.
	Note that the algorithm needs Eigen (https://eigen.tuxfamily.org/index.php?title=Main_Page).
	@author: Hisayoshi Muramatsu
	@date: 2023.05.02
*/

#ifndef DEF_TRAJECTORY_MODIFICATION
#define DEF_TRAJECTORY_MODIFICATION

#include <math.h>
#include <Eigen/Core>
#include <Eigen/LU>

using namespace Eigen;

class TraMod{
public:
	TraMod(
		double x_k,
		double w1,
		double w2,
		double x_max,
		double x_min,
		double dx_max,
		double dx_min,
		double ddx_max,
		double ddx_min,
		double epth,
		double kth,
		double T
	);
	// Functions
	double TrajectoryModification(double x_k);
	void SetParameters(double _w1, double _w2, double _x_max, double _x_min, double _dx_max, double _dx_min, double _ddx_max, double _ddx_min, double _T);

private:
	// Variables, vectors, and matrices
	double x_max, x_min, dx_max, dx_min, ddx_max, ddx_min;
	double epth, kth, T;
	double sig, gam, alpha_max, ep;
	const int n, l;
	double x_kp2, x_kp1, x_k, x_km1, x_km2, e_km1, e_km2;
	MatrixXd V, D, M;
	VectorXd x, c, r, lam, z, s, ones;
	Matrix2d Ass;
	Vector2d Bss, Css;
	// Functions
	double Update();
	double PositionConst(double phi_k, double phi_km1, double phi_km2);
	double sat(double in, double min, double max);

};

#endif
