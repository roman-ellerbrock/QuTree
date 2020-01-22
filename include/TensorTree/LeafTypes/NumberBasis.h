#pragma once
#include "Matrix.h"
#include "Vector.h"
#include "PrimitiveBasis.h"

class NumberBasis :
public PrimitiveBasis
{
 public:
	NumberBasis(int dim, bool fermion_);
	~NumberBasis();

	virtual void Initialize(double occ, double oSQR, double dummy2, double dummy3);

	void InitSPF(Tensorcd & phi)const;

	Tensorcd ApplyKin(const Tensorcd & phi)const;   //number operator
	Tensorcd ApplyP(const Tensorcd & phi)const;     //kreator
	Tensorcd applyX(const Tensorcd & phi)const;     //annihilator
	Tensorcd ApplyX2(const Tensorcd & phi)const=0;  //swaped number operator
	const Vectord& GetX()const { return x; } //does nothing
	Vectord& GetX() { return x; } //does nothing

	Tensorcd ToGrid(const Tensorcd& Acoeffs)const; //does nothing
	Tensorcd FromGrid(const Tensorcd& Acoeffs)const; //does nothing
  
	int oSQR()const{
		return osqrtree;}
	bool isfermion()const{return fermion;}
	bool HasDVR()const {return false; } // This tells that a number basis is not treated by the CDVR

 protected:
	// parameters that define the basis
	Vectord x;
	int startocc;
	int dim;
	int osqrtree;
	int minOcc;
	bool fermion; //true=fermion, false=boson
};
