//
// Created by Roman Ellerbrock on 7/27/21.
//

#include "ElectronicStructure.h"
#include "Util/RandomProjector_Implementation.h"

using namespace JordanWigner;

TwoIndex readTwoIndexIntegral(ifstream& file, size_t nOrbital) {
	TwoIndex hpq;
	for (size_t i = 0; i < nOrbital * nOrbital; ++i) {
		size_t p, q = 0;
		double h = 0.;
		file >> p >> q >> h;
		cout << p << " " << q << " | " << h << endl;
		hpq.push_back({p, q, h});
		// check whether hpq ends here by looking for an equal '=' sign
		auto pos = file.tellg();
		string peek;
		file >> peek;
		if (peek == "=") { break; }
		file.seekg(pos);
	}
	string t;
	file >> t;
	return hpq;
}

FourIndex readFourIndexIntegral(ifstream& file, size_t nOrbital) {
	FourIndex hpqrs;
	for (size_t i = 0; i < pow(nOrbital, 4); ++i) {
		size_t p, q, r, s = 0;
		double h = 0.;
		file >> p >> q >> r >> s >> h;
		hpqrs.push_back({p, q, r, s, h});

//		cout << p << " " << q << " " << r << " " << s << " | " << h << endl;

		auto pos = file.tellg();
		string peek;
		file >> peek;
		if (peek == "CI") { break; }
		file.seekg(pos);
	}
	return hpqrs;
}

template<typename T>
int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

size_t countelem(const Matrixd& A) {
	size_t count = 0;
	for (size_t i = 0; i < A.dim1(); ++i) {
		for (size_t j = 0; j < A.dim2(); ++j) {
			if (abs(A(i, j)) > 1e-10) {
				count++;
			}
		}
	}
	return count;
}

Matrixd sparsify(Matrixd A, double sqrtn, double eps) {
	if (A.dim1() != A.dim2()) {
		cerr << "nope.\n";
		exit(1);
	}
	cout << "number of elements before sparsification: " << countelem(A) << endl;
	mt19937 gen(12381823);
	uniform_real_distribution<double> dist(0., 1.);
	size_t num = 0;
	cout << sqrtn << " | n / 2 = " << sqrtn / 2 << endl;
	for (size_t i = 0; i < A.dim1(); ++i) {
		for (size_t j = i; j < A.dim2(); ++j) {
			if (abs(A(i, j)) > (eps / sqrtn) || (abs(A(i, j)) < 1e-15)) {
				continue;
			}
			/// let active space intact
			{
				size_t n = (int) sqrtn;
				size_t p = i % n;
				size_t q = i / n;
				size_t r = j % n;
				size_t s = j / n;
				size_t count = 0; // n_virtual
				n /= 2;
				if (p < n) { count++; }
				if (q < n) { count++; }
				if (r >= n) { count++; }
				if (s >= n) { count++; }
				if (count <= 2) { continue; } // only sparsify if more than 'count' e- are touched in virtual orbitals
				cout << p << " " << " " << q << " " << r << " " << s << " | " << count << " " << A(i, j) << endl;
			}
			double pr = abs(A(i, j)) * sqrtn / eps;
			if (dist(gen) < pr) {
				A(i, j) = sgn(A(i, j)) * eps / sqrtn;
			} else {
				if (abs(A(i, j)) > 1e-8) { num++; }
				A(i, j) = 0.;
			}
			A(j, i) = A(i, j);
		}
	}

	cout << "number of elements after sparsification: " << countelem(A) << endl;
	cout << "number saved: " << num << endl;

	getchar();
	return A;
}

Matrixd buildMat(const FourIndex& f, size_t n) {
	Matrixd h(n * n, n * n);
	for (const auto& x : f) {
		size_t p = get<0>(x);
		size_t q = get<1>(x);
		size_t r = get<2>(x);
		size_t s = get<3>(x);
		double he = get<4>(x);
		h(p + q * n, r + s * n) = he;
	}
	return h;
}

FourIndex fourIndex(const Matrixd& f, size_t n) {
	FourIndex h;
	double eps = 1e-12;
	for (size_t j = 0; j < f.dim2(); ++j) {
		for (size_t i = 0; i < f.dim1(); ++i) {
			if (abs(f(i, j)) > eps) {
				size_t p = i % n;
				size_t q = i / n;
				size_t r = j % n;
				size_t s = j / n;
				h.push_back(tuple<int, int, int, int, double>(p, q, r, s, f(i, j)));
			}
		}
	}
	return h;
}

void eigenvals(const Matrixd& h) {
	Matrixcd h2(h.dim1(), h.dim2());
	for (size_t i = 0; i < h.dim1() * h.dim2(); ++i) { h2[i] = h[i]; }
	mt19937 gen(1239129319);
	auto x = Random::diagonalizeRandom<complex<double>, Matrixcd, void>(h2, 100, 3, gen, nullptr);
	cout << "eigenvalues: " << endl;
	x.second.print();
}

SOPcd electronicStructure(const string& filename) {
	ifstream file(filename);
	if (file.bad()) {
		cerr << "Error: Cannot open integral file for electronic structure calculation.\n";
		exit(1);
	}
	string dump;
	size_t nOrbitals = 0;
	/// Read header, nElectrons, nOrbitals
	file >> dump >> dump >> dump >> nOrbitals;
	cout << "#Electrons = " << nOrbitals << endl;
	file >> dump >> dump >> dump >> nOrbitals;
	cout << "#Orbitals = " << nOrbitals << endl;

	/// Read header, E_HF, E_core
	double E = 0., Ecore = 0.;
	file >> dump >> dump >> dump >> E;
	file >> dump >> dump >> dump >> Ecore;
	cout << "E = " << E << " | E_core = " << Ecore << " | E-Ecore" << E - Ecore << endl;

	/// Read Hpq
	TwoIndex hpq = readTwoIndexIntegral(file, nOrbitals);
	cout << "#hpq = " << hpq.size() << endl;

	/// Read Hpqrs
	FourIndex hpqrs = readFourIndexIntegral(file, nOrbitals);
	cout << "#hpqrs = " << hpqrs.size() << endl;

	/// Sparsify attempt
	/*
	auto h = buildMat(hpqrs, nOrbitals);
	eigenvals(h);
	h = sparsify(h, nOrbitals, 1e-0);
	eigenvals(h);
	hpqrs = fourIndex(h, nOrbitals);
	*/

	return JordanWigner::electronicHamiltonian(hpq, hpqrs);
}

Matrixd convertTwoIndex(const TwoIndex& h) {
	size_t dim = 1;
	for (auto pq : h) {
		size_t p = get<0>(pq);
		if (p >= dim) { dim = p + 1; }
	}
	Matrixd ht(dim, dim);
	for (auto pq : h) {
		size_t p = get<0>(pq);
		size_t q = get<1>(pq);
		double val = get<2>(pq);
		ht(p, q) = val;
	}
	return ht;
}

Tensord convertFourIndex(const FourIndex& h) {
	size_t dim = 1;
	for (auto pqrs : h) {
		size_t p = get<0>(pqrs);
		if (p >= dim) { dim = p + 1; }
	}
	Tensord ht({dim, dim, dim, dim});
	for (auto pqrs : h) {
		size_t p = get<0>(pqrs);
		size_t q = get<1>(pqrs);
		size_t r = get<2>(pqrs);
		size_t s = get<3>(pqrs);
		double val = get<4>(pqrs);
		ht({p, q, r, s}) = val;
	}
	return ht;
}

