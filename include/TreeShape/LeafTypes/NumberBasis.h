#pragma once
#include "Core/Matrix.h"
#include "Core/Vector.h"
#include "LeafInterface.h"

class NumberBasis :
public LeafInterface
{
 public:
	NumberBasis(int dim, bool fermion);
	~NumberBasis();

	virtual void Initialize(double occ, double oSQR, double dummy2, double dummy3);

	void InitSPF(Tensorcd & phi)const;

	Tensorcd ApplyKin(const Tensorcd & phi)const;   //number operator
	Tensorcd ApplyP(const Tensorcd & phi)const;     //kreator
	Tensorcd applyX(const Tensorcd & phi)const;     //annihilator
	Tensorcd ApplyX2(const Tensorcd & phi)const=0;  //swaped number operator
	const Vectord& GetX()const { return x_; } //does nothing
	Vectord& GetX() { return x_; } //does nothing

	Tensorcd ToGrid(const Tensorcd& Acoeffs)const; //does nothing
	Tensorcd FromGrid(const Tensorcd& Acoeffs)const; //does nothing
  
	int oSQR()const{
		return osqrtree_;}
	bool isfermion()const{return fermion_;}
	bool HasDVR()const {return false; } // This tells that a number basis is not treated by the CDVR

 protected:
	// parameters that define the basis
	Vectord x_;
	int startOcc_;
	int dim_;
	int osqrtree_;
	int minOcc_;
	bool fermion_; //true=fermion_, false=boson
};
