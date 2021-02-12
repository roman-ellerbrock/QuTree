#pragma once
#include "Core/stdafx.h"
#include "AbstractNode.h"
#include "NodePosition.h"

#include "LeafTypes/LeafInterface.h"

// @TODO: Rename PhysPar
// @TODO: Use dictionary to store general leaf-memory?
class PhysPar
	/**
	 * \class Parameter container of a leaf.
	 */
{
public:
	PhysPar()
		: par0(0), par1(0), par2(0), par3(0) {}

	explicit PhysPar(istream& file)
		: par0(0), par1(0), par2(0), par3(0) {
		file >> par0;
		//assert(par0 != 0);
		file >> par1;
		file >> par2;
		file >> par3;
		//assert(par3 >= 0);
	}

	~PhysPar() = default;

	void setOmega(double omega) {
		par0 = omega;
	}

	void setR0(double r0) {
		par1 = r0;
	}

	void setWFR0(double wfr0) {
		par2 = wfr0;
	}

	void setWFOmega(double omega) {
		par3 = omega;
	}

	double Omega() const { return par0; }

	double R0() const { return par1; }

	double WFR0() const { return par2; }

	double WFOmega() const { return par3; }

	void info(ostream& os = cout) const {
		os << "par0=" << par0 << "\t";
		os << "par1=" << par1 << "\t";
		os << "par2=" << par2 << "\t";
		os << "par3=" << par3 << "\t";
		os << endl;
	}

protected:
	double par0, par1, par2, par3;
};

class Leaf
	: public AbstractNode
	/**
	 * \class Leaf
	 * \ingroup TTBasis
	 * \brief This class represents the leaf in the tree of a TTBasis.
	 *
	 * The leafs represent the lowest layer (even below the bottomlayer)
	 * in a tree of a TTBasis. Leaves contain abstract class pointers
	 * to PrimitiveBasis which provides the interface to the problem
	 * under consideration.
	 */
{
public:
	Leaf(istream& file, AbstractNode *up, NodePosition position);
	Leaf(size_t dim, size_t mode, size_t type, size_t subtype,
		PhysPar par);
	Leaf();
	Leaf(const Leaf&);
	Leaf(Leaf&&) = default;
	Leaf& operator=(const Leaf&);
	Leaf& operator=(Leaf&&) = default;
	~Leaf() override = default;

	void CreatePrimitiveBasis(size_t type, size_t subtype, size_t dim);

	void info(ostream& os = cout) const override;
	void Write(ostream& os = cout) const override;

	size_t nTotalNodes() const override { return 1; }

	size_t nNodes() const override { return 0; }

	size_t nLeaves() const override { return 1; }

	AbstractNode *nextNode() override { return this; }

	AbstractNode *nextNodeManthe() override { return this; }

	// Getter & Setter
	size_t Mode() const { return mode_; }

	int& Mode() { return mode_; }

	size_t Type() const { return type_; }

	size_t Dim() const { return dim_; }

	int type() const override { return nodeType_; }

	LeafInterface& PrimitiveGrid() { return *primitiveBasis_; }

	const LeafInterface& PrimitiveGrid() const { return *primitiveBasis_; }

	// This is not a GetNode& to avoid circular dependencies
	AbstractNode& Up() const { return *up_; };

	void SetPar(PhysPar par) { par_ = par; }

	PhysPar Par() const { return par_; }

	double Omega() const { return par_.Omega(); }

	double R0() const { return par_.R0(); }

	double WFR0() const { return par_.WFR0(); }

	double WFOmega() const { return par_.WFOmega(); }

	void update(const NodePosition& p) override;

	void UpdatePosition(const NodePosition& p);

	// Danger zone
	void setParent(AbstractNode *node) override { up_ = node; }

protected:
	int dim_, type_, mode_;
	int subType_;
	int nodeType_;
	// This is not a GetNode& to avoid circular dependencies
	AbstractNode *up_;
	PhysPar par_;
	NodePosition position_;
	unique_ptr<LeafInterface> primitiveBasis_;
};

