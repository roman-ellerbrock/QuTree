//
// Created by Roman Ellerbrock on 4/10/20.
//
#include <TreeClasses/TreeIO.h>
#include <TreeOperators/LeafOperator.h>
#include <TreeOperators/MultiLeafOperator.h>
#include "Measurements.h"
#include "TreeClasses/TensorTree.h"
#include "Circuits/GateOperators.h"


namespace Measurements {

	void Update(SparseMatrixTreecd& S, SparseMatrixTreecd& rho,
		const TensorTreecd& Psi,
		size_t next, size_t last, const Tree& tree) {

		MLOcd M;
		const auto& leaf = tree.getLeaf(next);
		M.push_back(identityMatrixcd(leaf.dim()), next);
		const auto& leaf_last = tree.getLeaf(last);
		M.push_back(identityMatrixcd(leaf_last.dim()), last);

		using namespace TreeFunctions;

		TreeFunctions::represent(S, M, Psi, tree);
		TreeFunctions::contraction(rho, Psi, S, tree);
	}

	Measurement measurement(TensorTreecd& Psi, mt19937& gen,
		const vector<size_t>& targets, const Tree& tree) {

		SparseMatrixTreecd overlap(targets, tree);
		SparseMatrixTreecd holeOverlap(targets, tree);
		uniform_real_distribution<double> dist(0., 1.);
		Measurement measures;

		assert(!targets.empty());
		size_t idx = 0;
		for (size_t coord : targets) {
			const Leaf& leaf = tree.getLeaf(coord);
			const Node& node = (Node&) leaf.parent();

			/// calculate occupancy of single-qubit state (|0>, |1>)
			size_t last = targets[idx];
			size_t next = targets[idx];
			if (idx + 1 < targets.size()) { next = targets[idx + 1]; }
			Update(overlap, holeOverlap, Psi, next, last, tree);
//			TreeFunctions::DotProduct(overlap, Psi, Psi, tree);
//			TreeFunctions::Contraction(holeOverlap, Psi, Psi, overlap, tree);
			auto occupancy = TreeIO::leafDensity(Psi, holeOverlap, leaf, tree);
			double norm = 0.;
			for (size_t i = 0; i < occupancy.dim1(); ++i) {
				norm += real(occupancy(i, i));
			}
			occupancy /= norm;

			double random_double = dist(gen);
			double acc = 0.;
			size_t i = 0;
			for (; i < occupancy.dim1(); ++i) {
				acc += real(occupancy(i, i));
				if (acc >= random_double) {
					break;
				}
			}
			assert(i < occupancy.dim1());
			measures.emplace_back(i);

			/// Project on single-qubit state
			shared_ptr<LeafOperatorcd> P = make_shared<Circuits::PrimitiveProjector>(i);
			MLOcd Proj;
			Proj.push_back(P, coord);
			Psi = Proj.apply(Psi, tree);

			idx++;
		}
		return measures;
	}

	Sample sample(const TensorTreecd& Psi, mt19937& gen,
		size_t n_samples, const vector<size_t>& targets, const Tree& tree) {

		Sample Map;
		for (size_t i = 0; i < n_samples; ++i) {
//            cout << " Sample no. :  " << i << endl;
			TensorTreecd Psi2(Psi);
			auto measure = measurement(Psi2, gen, targets, tree);
			Map[measure]++;
		}

		return Map;
	}

	void print(const Measurement& M, ostream& os) {
		for (auto xi : M) {
			cout << xi << "\t";
		}
	}

	ostream& operator<<(ostream& os, const Measurement& M) {
		print(M, os);
		return os;
	}

	void print(const Sample& S, ostream& os) {
		SampleVector sortedmap;
		for (auto it : S) {
			sortedmap.emplace_back(it);
		}
		sort(sortedmap.begin(), sortedmap.end(),
			[](pair<Measurement, size_t> a, pair<Measurement, size_t> b) { return a.second > b.second; });
		for (const auto& x : sortedmap) {
			print(x.first, os);
			os << ":\t" << x.second << endl;
		}
		os << endl;
	}

	Sample sample(const FullRank::Wavefunction& Psi, mt19937& gen,
		size_t number_samples) {

		Sample S;
		uniform_real_distribution<double> dist(0., 1.);

		for (size_t n = 0; n < number_samples; ++n) {
			double r = dist(gen);
			double sum = 0.;
			const TensorShape& shape = Psi.shape();
			size_t idx;
			for (idx = 0; idx < shape.totalDimension(); ++idx) {
				double p = pow(abs(Psi(idx)), 2);
				sum += p;
				if (sum >= r) {
					break;
				}
			}
			Measurement measure = indexMapping(idx, shape);
			S[measure]++;
		}
		return S;
	}
}

/// Cross entropy benchmarking
namespace Measurements {
	/// References:
	/// [1] (https://doi.org/10.1038/s41567-018-0124-x)

	/// gamma in [1]
	double eulers_constant = 0.5772156649;

	double H_U_exp(const Sample& S) {
		/// Eq. (18) in SI of Ref. [1]
		double H = 0;
		for (const auto& s: S) {
			size_t freq = s.second;
			double p = 1. / (1. * freq);
			H -= p * log(p);
		}
		return H;
	}

	size_t nSamples(const Sample& S) {
		/// Number of samples in S
		size_t num = 0;
		for (const auto& s : S) {
			num += s.second;
		}
		return num;
	}

	double sPorterThomas(size_t n_qubits) {
		size_t N = pow(2, n_qubits);
		double S = log((double) N) + eulers_constant - 1.;
		return S;
	}

	complex<double> probablity(const vector<size_t>& configuration, const FullRank::Wavefunction& Psi) {
		auto idx = indexMapping(configuration, Psi.shape());
		return Psi(idx);
	}

	Vectord probabilityDistibution(const FullRank::Wavefunction& Psi) {
		const TensorShape& shape = Psi.shape();
		Vectord p(shape.totalDimension());
		for (size_t i = 0; i < p.dim();++i) {
			p(i) = pow(abs(Psi(i)), 2);
		}
		return p;
	}

	/// Evaluate expected value of x over probability distribution x, i.e. E_p(x)
	double expectX(const Vectord& p, const Vectord& x) {
		double e = 0.;
		for (size_t i = 0; i < p.dim(); ++i) {
			e += p(i) * x(i);
		}
		return e;
	}

	double expectMomentX(const Vectord& p, Vectord x, size_t moment) {
		for (size_t i = 0; i < x.dim(); ++i) {
			x(i) = pow(x(i), moment);
		}
		return expectX(p, x);
	}

	double entropyPorterThomas(const FullRank::Wavefunction& Psi) {
		const TensorShape& shape = Psi.shape();
		size_t N = shape.totalDimension();
		return log((double) N) + eulers_constant - 1.;
	}

	double klDivergence(const Sample& X, const FullRank::Wavefunction& Psi,
		const function<double(const FullRank::Wavefunction& Psi)>& s, bool verbose) {
		double S = s(Psi);
		double H = crossEntropy(X, Psi, verbose);
		return S - H;
	}

	double entropy(const FullRank::Wavefunction& Psi) {
		const TensorShape& shape = Psi.shape();
		double S = 0.;
		for (size_t i = 0; i < shape.totalDimension(); ++i) {
			double p = pow(abs(Psi(i)), 2);
			p = max(p, 1e-16);
			S -= p * log(p);
		}
		return S;
	}

	double crossEntropy(const FullRank::Wavefunction& p, const FullRank::Wavefunction& q) {
		const TensorShape& shape = p.shape();
		double H = 0.;
		for (size_t i = 0; i < shape.totalDimension(); ++i) {
			H -= pow(abs(p(i)), 2) * log( pow(abs(q(i)), 2) );
		}
		return H;
	}

	double crossEntropy(const Sample& X, const FullRank::Wavefunction& Psi, bool verbose) {
		size_t N = 0;
		for (const auto& x : X) {
			N += x.second;
		}

		double H = 0.;
		for (const auto& x : X) {
			vector<size_t> measurement = x.first;
			size_t freq = x.second;

			double p_approx = ((double) freq) / ((double) N);
			double p_accurate = pow(abs(probablity(measurement, Psi)), 2);

			if (verbose) { cout << "p_approx: " << p_approx << ", p_acc: " << p_accurate << endl; }

			p_accurate = max(p_accurate, 1e-16);
			H -= p_approx * log(p_accurate);
		}

		return H;
	}

	double kullbackLeiblerDivergence(const FullRank::Wavefunction& p,
		const FullRank::Wavefunction& q) {
		double H = crossEntropy(p, q);
		double S = entropy(p);
		return H - S;
	}

	double crossEntropyDifference(const FullRank::Wavefunction& p,
		const FullRank::Wavefunction& q) {
		double H = crossEntropy(q, p);
		double S = entropy(p);
		return H - S;
	}

}
