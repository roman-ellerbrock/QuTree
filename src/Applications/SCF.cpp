//
// Created by Roman Ellerbrock on 11/14/22.
//

#include "Applications/SCF.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeClasses/HamiltonianRepresentation.h"
#include "TreeClasses/TreeIO.h"
#include "Core/TensorBLAS.h"
#include <iomanip> // for setprecision
#include <chrono>
using namespace std::chrono;

void addNodes(vector<const Node *>& sweep, const Node *p) {
	sweep.push_back(p);
	if (p->isBottomlayer()) {
		return;
	} else {
		for (int k = 0; k < p->nChildren(); ++k) {
			const Node *child = &(p->child(k));
			addNodes(sweep, child);
			sweep.push_back(p);
		}
	}
}

vector<const Node *> scf_sweep(const Tree& tree) {
	const Node *p = &tree.topNode();

	vector<const Node *> sweep;
	addNodes(sweep, p);
	sweep.push_back(nullptr);
	return sweep;
}

void apply(Tensorcd & HA,

const Tensorcd& A,

const SparseMatrixTreescd& hMats,

const SparseMatrixTreescd& hCons,

const SparseMatrixTreecd& hCorr,

const SparseMatrixTreecd& hConCorr,

const MatrixTreecd *rho,

const Hamiltonian& H,

const Node& node,

const vector<size_t>& n_actives, Tensorcd

* work) {

HA.
zero();

static Tensorcd hA;

if (hA.
shape()
!= A.
shape()
) {
hA = Tensorcd(A.shape());

}
for (

size_t l = 0;

l<hMats.
size();
++l) {
hA.
zero();

if (node.
isBottomlayer()
) {
/// only apply part if active hole
if (hMats[l].

isActive(node)

&& hCons[l].

isActive(node)

) {
const MLOcd& M = H[l];

const Leaf& leaf = node.getLeaf();

Tensorcd MA = M.apply(A, leaf);

matrixTensorBLAS(hA, hCons[l][node], MA, node

.

parentIdx(),

false);
HA += H.

coeff(l)

*

hA;

}
} else {
//			if (nActives(hMats[l], H[l], hCons[l], node) > 1) {
if (n_actives[l] > 1) {
TreeFunctions::apply(hA, hMats[l], & hCons[l], rho, A,
	hCons[l]

.

sparseTree(), node,

-1, work);
HA += H.

coeff(l)

*

hA;

}
}
}

if (!node.
isBottomlayer()
) {
/// add correction
for (

size_t k = 0;

k<node.
nChildren();
++k) {
const Node& child = node.child(k);

matrixTensorBLAS(HA, hCorr[child], A, child

.

childIdx(),

false);
}
}
if (!node.
isToplayer()
) {
matrixTensorBLAS(HA, hConCorr[node], A, node

.

parentIdx(),

false);
}
}

void apply(Tensorcd & HA,

const Tensorcd& A,

const HamiltonianRepresentation& hrep,

const Hamiltonian& H,

const Node& node

) {
const MatrixTreecd *rho = nullptr;

const SparseMatrixTreecd& hCorr = hrep.hCorr_;

const SparseMatrixTreecd& hConCorr = hrep.hConCorr_;

const vector<size_t>& n_actives = hrep.nActives_[node];

apply(HA, A, hrep

.hMats_, hrep.hContractions_,
hCorr, hConCorr, rho, H, node, n_actives, &(hrep.mem_.work2_));
}

Tensorcd apply(const Tensorcd& A, const HamiltonianRepresentation& hrep,
	const Hamiltonian& H, const Node& node) {
	Tensorcd HA(A.shape());
	apply(HA, A, hrep, H, node);
	return HA;
}

complex<double> fullContraction(const Tensorcd& A, const Tensorcd& B) {
	auto rho = A.dotProduct(B);
	return rho.trace();
}

double normalize(Tensorcd & A) {
	double norm = real(fullContraction(A, A));
	norm = sqrt(norm);
	A /= norm;
	return norm;
}

void testKrylovSpace(const KrylovSpace& space,
	const HamiltonianRepresentation& hrep,
	const Hamiltonian& H, const Node& node) {

	size_t dim = space.space_.size();
	Matrixcd S(dim, dim);

	for (size_t l = 0; l < dim; ++l) {
		for (size_t m = 0; m < dim; ++m) {
			S(l, m) = fullContraction(space.space_[l], space.space_[m]);
		}
	}
	cout << "|1 - S| = " << residual(identityMatrixcd(dim), S) << endl;

	Matrixcd hm(dim, dim);
	for (size_t l = 0; l < dim; ++l) {
		auto HPsi = apply(space.space_[l], hrep, H, node);
		for (size_t m = 0; m < dim; ++m) {
			hm(l, m) = fullContraction(HPsi, space.space_[m]);
		}
	}
	cout << "<<H>> = \n";
	(hm * QM::cm).print();
}

void solveKrylovSpace(KrylovSpace& krylov, Tensorcd Psi, const HamiltonianRepresentation& hrep,
	const Hamiltonian& H, const Node& node, size_t krylov_size, double conversion) {

	vector<Tensorcd>& krylovSpace = krylov.space_;

	/// normalize wavepacket
	auto norm = normalize(Psi);
	if (abs(norm - 1.) > 1e-5) {
		cerr << "Norm != 1.\n";
		cerr << "norm: " << norm << ", renormalizing." << endl;
//		exit(1);
	}
	krylovSpace[0] = Psi;

	vector<double> beta;
	vector<double> alpha;

	/// Do the first step separately, since no beta exists
	Tensorcd HPsi = apply(Psi, hrep, H, node);

	/// Calculate first alpha
	/// Note: in general <Psi_i|H|Psi_j> is a matrix if multiple states are stored in Psi
	complex<double> alphamat = fullContraction(Psi, HPsi);

	/// We use only a single state, so we read the only element of the 1x1 matrix
	alpha.push_back(real(alphamat));

	/// First step is simpler, because beta = 0
	Psi = HPsi - alpha.back() * Psi;

	/// remaining steps
	for (size_t l = 1; l < krylov_size; ++l) {
		beta.push_back(normalize(Psi));
		krylovSpace[l] = Psi;

		apply(HPsi, Psi, hrep, H, node);
		alphamat = fullContraction(Psi, HPsi);
		alpha.push_back(real(alphamat));

		const Tensorcd& v = krylovSpace[l];
		const Tensorcd& vlast = krylovSpace[l - 1];
		Psi = HPsi - alpha.back() * v - beta.back() * vlast;
	}

	/// Build tri-diagonal hamiltonian matrix
	Matrixcd Hrep(krylov_size, krylov_size);
	for (size_t l = 0; l < krylov_size; ++l) {
		Hrep(l, l) = alpha[l];
		if (l == krylov_size - 1) { break; }
		Hrep(l + 1, l) = beta[l];
		Hrep(l, l + 1) = beta[l];
	}
//	cout << "<H>:\n";
//	(Hrep * QM::cm).print();

	/// Diagonalize Hamiltonian matrix and save transformation matrix U and eigenvalues E.
	krylov.spectrum_ = diagonalize(Hrep);

	/// output
	auto x = krylov.spectrum_.second * conversion;
	cout << x(0) << endl;
}

void imaginaryTimePropagation(Tensorcd & Psi,

const KrylovSpace& kry,

double beta

) {
const Matrixcd& U = kry.spectrum_.first;

const Vectord& E = kry.spectrum_.second;

size_t krylov_size = E.dim();

Vectorcd c(krylov_size);

for (

size_t l = 0;

l<krylov_size;

++l) {
c(l)

+=

U(l,

0) *

conj(U(0, 0)

);
for (

size_t k = 0;

k<krylov_size;

++k) {
complex<double> exp_k = exp(-beta * (E(k) - E(0)));

c(l)

+=

U(l, k

) *

exp_k *conj(U(0, k));

}
}

Psi.
zero();
for (

size_t l = 0;

l<krylov_size;

++l) {
const Tensorcd& psi_l = kry.space_[l];

Psi +=

c(l)

*

psi_l;

}

normalize(Psi);

}

void diagonalizeKrylov(Tensorcd & Psi,

const KrylovSpace& kry

) {
const Matrixcd& U = kry.spectrum_.first;

const Vectord& E = kry.spectrum_.second;

size_t krylov_size = kry.space_.size();

Psi.
zero();
for (

size_t l = 0;

l<krylov_size;

++l) {
const Tensorcd& psi_l = kry.space_[l];

Psi +=

U(l,

0) *

psi_l;

}
}

size_t adjacentIndex(const Node& from, const Node *to) {
	/// returns index pointing from node *from* to node *to*
	/// caution when touching the logic/order of conditions here!
	if (to == nullptr) {
		return from.parentIdx();
	}
	for (size_t k = 0; k < from.nChildren(); ++k) {
		if (from.child(k).address() == to->address()) {
			return k;
		}
	}
	if (from.parent().address() == to->address()) {
		return from.parentIdx();
	}
	cerr << "Node *to* is not adjacent to node *from*.\n";
	exit(1);
	return 0;
}

void outputSCF(Tensorcd Psi, const HamiltonianRepresentation& hrep,
const Hamiltonian& H, const Node& node) {

auto HPsi = apply(Psi, hrep, H, node);

cout << setprecision(12);
cout << "<H> = " << real(fullContraction(Psi, HPsi) ) * QM::cm << endl;

}

void scf(SCF_parameters& par) {
	cout << "# ===============    SCF    ===============\n";
	cout << "# nIter = " << par.nIter << endl;
	cout << "# nKrylov = " << par.nKrylov << endl;
	cout << "# nITP = " << par.nITP << endl;
	cout << "# beta = " << par.beta << endl;
	double time = 0;
	size_t krylov_size = par.nKrylov;
	double beta = par.beta;
	double conversion = par.conversion;
	const Tree& tree = *par.tree;
	const Hamiltonian& H = *par.h;
	TensorTreecd& Psi = *par.psi;

	auto sweeper = scf_sweep(tree);
	double eps = 1. / ((double) (sweeper.size() - 1)); /// just for printing output

	vector<KrylovSpace> spaces;
	for (const Node& node: tree) {
		spaces.emplace_back(KrylovSpace(node.shape(), par.nKrylov));
	}

	HamiltonianRepresentation hrep(H, tree, tree, false);
	cout << setprecision(12);
	for (const Node& node: tree) {
		if (!node.isToplayer()) {
			const Node& next = node.parent();
			hrep.buildSCF(H, Psi, node, next, time);
		}
	}
	mt19937 gen(192931923);

	auto start = chrono::high_resolution_clock::now();

	for (size_t it = 0; it < par.nIter; ++it) {
		for (size_t l = 0; l < sweeper.size() - 1; ++l) {
			if (sweeper[l] == nullptr) { break; }

			/// Build Krylov Space
			const Node& node = *sweeper[l];
			size_t tot = node.shape().totalDimension();
			size_t ksize = krylov_size > tot ? tot : krylov_size;
			cout << it + l * eps << " ";

			KrylovSpace& krylov = spaces[node.address()];
			solveKrylovSpace(krylov, Psi[node], hrep, H, node, ksize, conversion);
			if (it < par.nITP) {
				imaginaryTimePropagation(Psi[node], krylov, beta);
			} else {
				diagonalizeKrylov(Psi[node], krylov);
			}

			/// isometrize Psi towards next node
			const Node *next_ptr = sweeper[l + 1];
			size_t outIdx = adjacentIndex(node, next_ptr);
			Tensorcd PsiW = Psi[node];
			Psi[node] = qr(Psi[node], outIdx);

			/// overlap with <I | A> and multiply into adjacent node
			if (next_ptr == nullptr) { break; }
			const Node& next = *next_ptr;
			auto r = contraction(Psi[node], PsiW, outIdx);
			size_t inIdx = adjacentIndex(next, &node);
			Psi[next] = matrixTensor(r, Psi[next], inIdx);

			/// rebuild matrix elements pointing towards next node
			hrep.buildSCF(H, Psi, node, next, time);
		}
		if (par.output) { TreeIO::output(Psi, tree); }
//		if (krylov_size < par.nKrylov) { krylov_size++; }
	}
	auto stop = chrono::high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	auto duration_sec = duration_cast<milliseconds>(stop - start);
	cout << "Time taken by function: \n"
		 << duration.count() / 1000e0 << " ms\n"
		 << duration_sec.count() / 1000e0 << " s\n"
		 << duration_sec.count() / (1000e0 * par.nIter) << " s/iteraton." << endl;
}
