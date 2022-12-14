//
// Created by Roman Ellerbrock on 12/12/21.
//

#include "NDIHamiltonian.h"


namespace Operator {

	// !!! IMPORTANT !!!
	double cutoff_acc = 1e-3;

	std::string fix(const std::string& in) {
		return std::regex_replace(
			in,
			std::regex("(?:([a-zA-Z])([0-9]))|(?:([0-9])([a-zA-Z]))"),
			"\\1\\3 \\2\\4",
			std::regex_constants::format_sed
		);
	}

	void toPauli(MLOcd& M, string op, size_t t) {
		if (op == "X") {
			M.push_back(sigma_x(), t);
		} else if (op == "Z") {
			M.push_back(sigma_z(), t);
		} else {
			cerr << "unknown operator for pauli conversion.\n";
			exit(1);
		}
	}

	void readLineMPO(SOPcd& S, string line, size_t nprod, size_t nFragments) {
		std::replace(line.begin(), line.end(), '*', ' ');
		stringstream ss(line);

		double coeff = 0;
		ss >> coeff;
		bool out = abs(coeff) > cutoff_acc;
		if (out) { cout << coeff; }

		line = line.erase(0, line.find(' '));
		line = fix(line);
		ss = stringstream(line);
		MLOcd M;
		for (size_t i = 0; i < nprod; ++i) {
			string op;
			ss >> op;
			size_t t;
			ss >> t;
			if (out) { cout << " * " << op << t; }
			toPauli(M, op, t);
		}
		for (size_t i = 0; i < nprod; ++i) {
//			if (M.mode(i) >= 13) {
//			if (M.mode(i) >= 52) {
			if (M.mode(i) >= nFragments) {
				return;
			}
			if (out) { cout << endl; }
		}
		if (abs(coeff) > cutoff_acc) {
			S.push_back(M, coeff);
		}
	}

	SOPcd ndi() {
		string filename = "pauli_156q.txt";
//		string filename = "pauli_13q.txt";
		ifstream file(filename);
		size_t nFragments = 156; // number of active fragments
		size_t n_single = 312;
		size_t n_double = 59442;
		size_t nline = 1 + n_single + n_double;
		SOPcd S;
		string line;

		file >> cutoff_acc;
		file >> nFragments;
		file >> n_single;
		file >> n_double;
		cout << "cutoff: " << cutoff_acc << endl;
		cout << "nFragments: " << nFragments << endl;
		cout << "n_single: " << n_single << endl;
		cout << "n_double: " << n_double << endl;

		/// Energy shift
		double shift = 0.;
		file >> shift;
		getline(file, line);
		MLOcd I;
		shift = shift / 156.;
		cout << "E-shift: " << shift << endl;
		for (size_t n = 0; n < nFragments; ++n) {
			I.push_back(identityMatrixcd(2), n);
//			double shift = 0.5914818726393506 / (double) nFragments;
//			double shift = 8.56655 / 52.;
//			double shift = 6.562775130472779 / (double) nFragments;
//			double shift = 7.000638078882831 / (double) nFragments;
			S.push_back(I, shift);
		}

		for (size_t i = 0; i < n_single; ++i) {
			getline(file, line);
			readLineMPO(S, line, 1, nFragments);
		}

		for (size_t i = 0; i < n_double; ++i) {
			getline(file, line);
			readLineMPO(S, line, 2, nFragments);
		}

		cout << "S.size() = " << S.size() << endl;
		return S;
	}
}
