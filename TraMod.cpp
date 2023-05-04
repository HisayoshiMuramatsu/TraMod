/**
	Position trajectory modification for position, velocity and acceleration constraints.
	Note that the algorithm needs Eigen (https://eigen.tuxfamily.org/index.php?title=Main_Page).
	@author: Hisayoshi Muramatsu
	@date: 2023.05.02
*/

#include "TraMod.hpp"

TraMod::TraMod(
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
	):
	x_max(x_max), x_min(x_min), dx_max(dx_max), dx_min(dx_min), ddx_max(ddx_max), ddx_min(ddx_min), epth(epth), kth(kth), T(T),
	sig(0.1), gam(0.99), alpha_max(10), ep(0), n(3), l(12),
	x_kp2(x_k), x_kp1(x_k), x_k(x_k), x_km1(x_k), x_km2(x_k), e_km1(0), e_km2(0),
	V(n,n), D(l,3), M(n+l,n+l),
	x(n), c(n), r(l), lam(l), z(n+l), s(l), ones(l),
	Ass(2,2),
	Bss(2), Css(2)
{

	x=VectorXd::Zero(n);
	c=VectorXd::Zero(n);

	MatrixXd a0=MatrixXd::Zero(2,3), a1=MatrixXd::Zero(2,3), a2=MatrixXd::Zero(2,3);
	a0 << 0,0,1/T,
	      0,0,1/T/T;
	a1 << 0,1/T,-1/T,
	      0,1/T/T,-2/T/T;
	a2 << 1/T,-1/T,0,
	      1/T/T,-2/T/T,1/T/T;
	MatrixXd A=MatrixXd::Zero(6,3);
	for(int i(0); i<2; i++){
		for(int j(0); j<3; j++){
			A(  i,j)=a2(i,j);
			A(2+i,j)=a1(i,j);
			A(4+i,j)=a0(i,j);
		}
	}

	V=MatrixXd::Identity(n,n);
	V<<w2,-w2,0,
	   -w2,w1+w2,-w1,
	   0,-w1,1+w1;
	D=MatrixXd::Zero(l,3);
	for(int i(0); i<6; i++){
		for(int j(0); j<3; j++){
			D(  i,j)= A(i,j);
			D(6+i,j)=-A(i,j);
		}
	}
	MatrixXd Dt=D.transpose();
	M=MatrixXd::Zero(n+l,n+l);
	for(int i(0); i<n; i++){
		for(int j(0); j<n; j++) M(i,j)   = V(i,j);
		for(int j(0); j<l; j++) M(i,n+j) = Dt(i,j);
	}

	x=VectorXd::Zero(n);
	c=VectorXd::Zero(n);
	r=VectorXd::Zero(l);
	lam=VectorXd::Ones(l);
	z=VectorXd::Zero(n+l);
		for(int i(0); i<n; i++) z(i)  =x(i);
		for(int i(0); i<l; i++) z(n+i)=lam(i);
	s=VectorXd::Zero(l);
		s=r-D*x;
	ones=VectorXd::Ones(l);

	Ass=Matrix2d::Zero();
		Ass << 1, T,
			   0, 1;
	Bss=Vector2d::Zero();
		Bss << 0, T;
	Css=Vector2d::Zero();
		Css << 1, 0;
}

double TraMod::TrajectoryModification(double x_k_in){

	// -------- Buffer --------
	x_km2 = x_km1;
	x_km1 = x_k;
	x_k   = x_kp1;
	x_kp1 = x_kp2;
	x_kp2 = x_k_in;
	// -------- buffer --------

	// -------- Definitions of vectors --------
	VectorXd phi_max=VectorXd(2), phi_min=VectorXd(2);
		phi_max << dx_max, ddx_max;
		phi_min << dx_min, ddx_min;
	VectorXd PHI_max=VectorXd::Zero(6), PHI_min=VectorXd::Zero(6);
	for(int i(0); i<2; i++){
		PHI_max(  i)=phi_max(i);
		PHI_max(2+i)=phi_max(i);
		PHI_max(4+i)=phi_max(i);
		PHI_min(  i)=phi_min(i);
		PHI_min(2+i)=phi_min(i);
		PHI_min(4+i)=phi_min(i);
	}
	VectorXd b0=VectorXd::Zero(2), b1=VectorXd::Zero(2), b2=VectorXd::Zero(2);
	b0 << -(1/T)*e_km1 + (1/T)*(x_k-x_km1),
	      -(1/T/T)*(2*e_km1-e_km2) + (1/T/T)*(x_k-2*x_km1+x_km2);
	b1 << (1/T)*(x_kp1-x_k),
	      (1/T/T)*e_km1 + (1/T/T)*(x_kp1-2*x_k+x_km1);
	b2 << (1/T)*(x_kp2-x_kp1),
	      (1/T/T)*(x_kp2-2*x_kp1+x_k);
	VectorXd B=VectorXd::Zero(6);
	for(int i(0); i<2; i++){
		B(  i)=b2(i);
		B(2+i)=b1(i);
		B(4+i)=b0(i);
	}
	for(int i(0); i<6; i++){
		r(  i)= PHI_max(i)-B(i);
		r(6+i)=-PHI_min(i)+B(i);
	}
	// -------- Definitions of vectors --------

	// -------- Initial values --------
	// x
	x(2)= - T*T*b0(1); 			 // <-> ddx_{k}+dde_{k}=0 (b0(1)=ddx_{k}+(-2e_{k-1}+e_{k-2})/T^2)
	x(1)= 2*x(2)-T*T*b1(1); 	 // <-> ddx_{k+1}+dde_{k+1}=0 (b1(1)=ddx_{k+1}+e_{k-1}/T^2)
	x(0)= 2*x(1)-x(2)-T*T*b2(1); // <-> ddx_{k+2}+dde_{k+2}=0 (b2(1)=ddx_{k+2})
	// lambda
	Matrix3d DtD = (D.transpose())*D;
	Matrix3d DtDinv = DtD.inverse();
	lam=-D*DtDinv*V*x;
	for(int i(0);i<12;i++){
		if(lam(i)<0) lam(i)=1e-10; // lam needs to be positive; it is modified to be a positive small value if it is negative.
	}
	// -------- Initial values --------

	// -------- Check of initial constraints --------
	VectorXd Dxmr=VectorXd::Zero(l);
	double Dxmr_max(0);
	Dxmr=D*x-r;
	Dxmr_max=Dxmr(0);
	for(int i(0); i<l; i++){
		if(Dxmr_max<Dxmr(i)) Dxmr_max=Dxmr(i); // The constraints are satisfied if Dxmr_max is negative or equal to zero: Dx-r<=0 <-> Dx<=r.
	}
	// -------- Check of initial constraints --------

	// -------- QP with Velocity & Acceleration Constraints --------
	int k(0);
	double e_k(0), ep(0);
	if(Dxmr_max>0){ // If Dxmr_max is not negative, the previous solution: x_k+e_k=x_{k-1} is used.
		std::cout << "Tra. Mod. has an initial state error." << std::endl;
		e_k=x_km1-x_k;
	}else{
		for(int i(0); i<n; i++) z(i)  =x(i);
		for(int i(0); i<l; i++) z(n+i)=lam(i);
		s=r-D*x;
		ep=s.dot(lam)/((double)l);
		do{
			ep=Update();
			k++;
		}while(ep>epth&&k<kth);
		e_k=x(2);
	}
	// -------- QP with Velocity & Acceleration Constraints --------

	// -------- Position Constraints --------
	if((int)round(100000*x_min)==0&&(int)round(100000*x_max)==0){ // If x_min=x_max=0, this trajectory modification only imposes the velocity and acceleration constraints without position constraints.
	}else{
		e_k=PositionConst(x_k+e_k, x_km1+e_km1, x_km2+e_km2);
	}
	// -------- Position Constraints --------

	// -------- buffer --------
	e_km2=e_km1;
	e_km1=e_k;
	// -------- buffer --------

	return x_k + e_k;

}

void TraMod::SetParameters(double _w1, double _w2, double _x_max, double _x_min, double _dx_max, double _dx_min, double _ddx_max, double _ddx_min, double _T){
	double w1 = _w1;
	double w2 = _w2;
	x_max	= _x_max;
	x_min	= _x_min;
	dx_max	= _dx_max;
	dx_min	= _dx_min;
	ddx_max = _ddx_max;
	ddx_min = _ddx_min;
	T		= _T;
	V<<w2,-w2,0,
	   -w2,w1+w2,-w1,
	   0,-w1,1+w1;
}

double TraMod::Update(){

	MatrixXd S=s.asDiagonal();
	MatrixXd Lam=lam.asDiagonal();

	MatrixXd LamD=Lam*D;
	for(int i(0); i<l; i++){
		for(int j(0); j<n; j++) M(n+i,j)   = -LamD(i,j);
		for(int j(0); j<l; j++) M(n+i,n+j) = S(i,j);
	}

	VectorXd v1=-V*x-c-(D.transpose())*lam;
	VectorXd v2=sig*ep*ones-S*lam;

	VectorXd v=VectorXd::Zero(n+l);
	for(int i(0); i<n; i++) v(i)   = v1(i);
	for(int i(0); i<l; i++) v(n+i) = v2(i);

	FullPivLU<MatrixXd> lu(M); // FullPivLU
	VectorXd Delta = lu.solve(v); // Solve the simultaneous equations

	VectorXd Delta_x=VectorXd::Zero(n), Delta_lam=VectorXd::Zero(l);
	for(int i(0); i<n; i++) Delta_x(i)=Delta(i);
	for(int i(0); i<l; i++) Delta_lam(i)=Delta(n+i);
	VectorXd Delta_s=-D*Delta_x;

	double alpha(alpha_max), alpha_d(0);
	for(int i(0); i<l; i++){
		if(Delta_s(i)<0){
			alpha_d = - s(i)/Delta_s(i);
			if(alpha_d < alpha) alpha = alpha_d;
		}
		if(Delta_lam(i)<0){
			alpha_d = - lam(i)/Delta_lam(i);
			if(alpha_d < alpha) alpha = alpha_d;
		}
	}
	alpha*=gam;

	z+=alpha*Delta;
	for(int i(0); i<n; i++) x(i)	=z(i);
	for(int i(0); i<l; i++) lam(i)	=z(n+i);
	s=r-D*x;
	ep=s.dot(lam)/((double)l);

	return ep;

}

double TraMod::PositionConst(double phi_k, double phi_km1, double phi_km2){
	double dphi_k=(phi_k-phi_km1)/T;

	int L(0);
	double bddx_k(0), uddx_k(0);
		if(dphi_k>0){
			L=(int)ceil(-(dphi_k+T*ddx_max-T*ddx_min)/(T*ddx_min));
			bddx_k=ddx_max; uddx_k=ddx_min;
		}else if(dphi_k<0){
			L=(int)ceil(-(dphi_k+T*ddx_min-T*ddx_max)/(T*ddx_max));
			bddx_k=ddx_min; uddx_k=ddx_max;
		}

	double hphi_kpL(0);
	Vector2d Phi(phi_k,dphi_k);
		if(L==0){
			hphi_kpL=phi_k;
		}else if(L==1){
			hphi_kpL=Css.dot(Ass*Phi)+Css.dot(Bss)*bddx_k;
		}else{
			Matrix2d AL=Ass, AL1=Ass, Ai=Matrix2d::Identity();
				for(int i(1);i<L;i++) AL*=Ass;
				for(int i(1);i<L-1;i++) AL1*=Ass;
			hphi_kpL=Css.dot(AL*Phi)+Css.dot(AL1*Bss)*bddx_k;
			for(int i(0);i<L-1;i++){
				hphi_kpL+=Css.dot(Ai*Bss)*uddx_k;
				Ai*=Ass;
			}
		}

	if(x_min<=hphi_kpL&&hphi_kpL<=x_max){
		return phi_k - x_k;
	}else{
		return 2*phi_km1 - phi_km2 + T*T*sat(-(phi_km1-phi_km2)/T/T,ddx_min,ddx_max) - x_k;
	}

}

double TraMod::sat(double in, double min, double max){
	if(in>=max){
		return max;
	}else if(in<=min){
		return min;
	}else{
		return in;
	}
}
