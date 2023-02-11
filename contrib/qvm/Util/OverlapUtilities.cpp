//
// Created by Roman Ellerbrock on 1/19/21.
//
#include "TreeClasses/TensorTreeFunctions.h"
#include "Measurements.h"
#include "Circuits/GateOperators.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "Util/filenames.h"
#include "Util/OverlapUtilities.h"

namespace Utility {

	MLOcd projectMeasurement(const vector<size_t>& measurement) {
		MLOcd M;
		for (size_t k = 0; k < measurement.size(); ++k) {
			auto l = new Circuits::PrimitiveProjector(measurement[k]);
			shared_ptr<Circuits::PrimitiveProjector> Pk(l);
			M.push_back(Pk, k);
		}
		return M;
	}

	double calcXEB(const vector<double>& ps, size_t n_qubits) {
		double xeb;
		size_t M = ps.size();
		for (size_t i = 0; i < M; ++i) {
			xeb += ps[i];
		}
		xeb /= (double) M;
		xeb *= pow(2, n_qubits);
		xeb -= 1.;
		return xeb;
	}

	void wavefunctionOverlap(const Tree& tree, const vector<TensorTreecd>& Psi,
		const vector<TensorTreecd>& Chi) {
		const Node& top = tree.topNode();
		for (size_t i = 0; i < Psi.size(); ++i) {
			for (size_t j = 0; j < Chi.size(); ++j) {
				MatrixTreecd S(tree);
				for (const Node& node : tree) {
					size_t n = Psi[i][node].shape().lastDimension();
					size_t m = Chi[j][node].shape().lastDimension();
					S[node] = Matrixcd(n, m);
				}
				TreeFunctions::dotProduct(S, Psi[i], Chi[j], tree);
				auto& s = S[top];
				cout << i << "\t" << j << "\t" << pow(abs(s(0,0)),2) << "\t";
				s.print();
			}
		}
	}

	void statisticalWavefunctionOverlap(const Tree& tree) {
		/// Read a reference and a sample of wavefunctions and check fidelity
		size_t sample = 4;
		complex<double> s_avg;
		TensorTreecd Psi("Psi.fr.dat");
		for (size_t s = 0; s < sample; ++s) {
			string filename = "Psi." + to_string(s) + ".dat";
			TensorTreecd Psis(filename);
			auto S = TreeFunctions::dotProduct(Psi, Psis, tree);
			const Node& top = tree.topNode();
			complex<double> lap = S[top](0, 0);
			s_avg += lap;
			cout << s << "\t" << lap << endl;
		}
		cout << "S = " << abs(s_avg) / ((double) sample) << endl;
	}

	void xeb(const Tree& fr_tree, const Tree& xeb_tree, const string& file_fr,
		const string& file_xeb) {
		/// Read wavefunctions
		const TensorTreecd Psi(file_fr);
		const TensorTreecd Chi(file_xeb);

		/// some required variables
		mt19937 gen(time(NULL));
		vector<size_t> targets;
		for (size_t i = 0; i < xeb_tree.nLeaves(); ++i) { targets.emplace_back(i); }
		const Node& top = fr_tree.topNode();
		ofstream of("xeb.dat");

		/// Perform sampling
		size_t nsample = 10000;
		vector<double> probs;
		for (size_t s = 0; s < nsample; ++s) {
			auto Chi2 = Chi;
			auto meas = Measurements::measurement(Chi2, gen, targets, xeb_tree);
			MLOcd M = projectMeasurement(meas);
			auto MPsi = M.apply(Psi, fr_tree);
			auto MS = TreeFunctions::dotProduct(MPsi, MPsi, fr_tree);
			double Mprob = abs(MS[top](0, 0));
			probs.push_back(Mprob);
			double xeb = calcXEB(probs, xeb_tree.nLeaves());

			cout << s << ", |<m_i|Psi>|^2 = " << Mprob << ", xeb = " << xeb << endl;
			of << s << "\t" << xeb << "\t" << Mprob << endl;
		}
	}

	double calcProbability(const Measurements::Measurement& meas, const TensorTreecd& Psi,
		const Tree& tree) {
		MLOcd M = projectMeasurement(meas);
		auto MPsi = M.apply(Psi, tree);
		auto MS = TreeFunctions::dotProduct(MPsi, MPsi, tree);
		return abs(MS[tree.topNode()](0, 0));
	}

	void xeb_stat(ostream& os, mt19937& gen, const TensorTreecd& Psi,
		const statistical::Wavefunctions& Chis, const Tree& tree,
		size_t nsample) {
		vector<size_t> targets;
		for (size_t i = 0; i < tree.nLeaves(); ++i) { targets.emplace_back(i); }
		vector<double> ps;
		for (size_t s = 0; s < nsample; ++s) {
			auto Chis2 = Chis;
			auto meas = statistical::measurement(Chis2, gen, targets, tree);
			auto p = calcProbability(meas, Psi, tree);
			ps.push_back(p);
			double xeb = calcXEB(ps, tree.nLeaves());
			cout << s << ", |<m_i|Psi>|^2 = " << p << ", xeb = " << xeb << endl;
			os << s << "\t" << xeb << "\t" << p << endl;
		}
	}


}
