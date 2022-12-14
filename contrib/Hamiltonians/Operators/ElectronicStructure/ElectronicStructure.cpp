//
// Created by Roman Ellerbrock on 7/27/21.
//

#include "ElectronicStructure.h"

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

		cout << p << " " << q << " " << r << " " << s << " | " << h << endl;

		auto pos = file.tellg();
		string peek;
		file >> peek;
		if (peek == "CI" || peek == "=") { break; }
		file.seekg(pos);
	}
	return hpqrs;
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
	file >> dump >> dump >> dump >> dump >> nOrbitals;
	cout << "#Electrons = " << nOrbitals << endl;
	file >> dump >> dump >> dump >> dump >> dump;
	file >> dump >> dump >> dump >> nOrbitals;
	cout << "#Orbitals = " << nOrbitals << endl;

	/// Read header, E_HF, E_core
	double E = 0., Ecore = 0.;
	file >> dump >> dump >> dump >> E;
	file >> dump >> dump >> dump >> Ecore;
	cout << "E = " << E << " | E_core = " << Ecore << " | E-Ecore" << E - Ecore << endl;
	double Esinglet = 0.;
	file >> dump >> dump >> dump >> dump >> dump >> Esinglet;
	cout << "E_singlet = " << Esinglet << endl;
	cout << "E_singlet - E_core = " << Esinglet - Ecore << endl;

	/// Read Hpq
	TwoIndex hpq = readTwoIndexIntegral(file, nOrbitals);
	cout << "#hpq = " << hpq.size() << endl;

	/// Read Hpqrs
	FourIndex hpqrs = readFourIndexIntegral(file, nOrbitals);
	cout << "#hpqrs = " << hpqrs.size() << endl;

	cerr << "WARNING:\n";
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

void tucker(Tensord h);

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
//	tucker(ht);
	return ht;
}

void tucker(Tensord h) {
	for (size_t k = 0; k < h.shape().order(); ++k) {
		auto rho = contraction(h, h, k);
		auto xev = diagonalize(rho);
		cout << "=== k : " << k << " ===\n";
		xev.second.print();
	}

	Tensord g = h;
	for (size_t I = 0; I < h.shape().totalDimension(); ++I) {
		auto idx = indexMapping(I, h.shape());
		swap(idx[1], idx[2]);
		h(idx) = g(I);
	}
	size_t n = h.shape()[0];
	TensorShape shape({n*n, n*n});
	h.reshape(shape);
	auto rho = contraction(h, h, 0);
	auto xev = diagonalize(rho);
	xev.second.print();
	getchar();
}