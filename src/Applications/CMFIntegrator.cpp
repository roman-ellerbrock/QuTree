//
// Created by Roman Ellerbrock on 3/8/20.
//

#include "Applications/CMFIntegrator.h"
#include <chrono>
#include "TreeClasses/TreeIO.h"

constexpr bool eom_spf = true;

CMFIntegrator::CMFIntegrator(const Hamiltonian& H,
	const Tree& tree, const Tree& cdvrtree, complex<double> phase)
	: matrices_(H, tree, cdvrtree), tconst_(0.25), max_increase_(2.25),
	  cdvr_nr_(0) {
	for (const Node& node : tree) {
		interfaces_.emplace_back(LayerInterface(H, matrices_, node, phase));

		// create an empty object for the integrator
		const TensorShape& shape = node.shape();
		Tensorcd empty(shape);
		bs_integrators_.emplace_back(bs_integrator(shape.totalDimension(), empty));
	}
}

void CMFIntegrator::Integrate(IntegratorVariables& job, ostream& os) {
	using namespace chrono;
	high_resolution_clock::time_point t1;
	high_resolution_clock::time_point t2;

	microseconds mattime(0);
	microseconds steptime(0);

	// read work parameters
	double time = job.time_now;
	double dt = job.dt;
	double timeend = job.time_end;
	double out = job.out;

	double accuracy_CMF = job.accuracy_root;
	double accuracy_BS = job.accuracy_leaf;

	const Tree& tree = *job.tree;
	Wavefunction& Psi = *job.psi;
	const Hamiltonian& H = *job.h;

	Psi.write(job.ofname, job.append_);

	// Reinitiate initia; steps
	double timepart = 0.25;
	dt_bs_.clear();
	for (const Node& node : tree) {
		dt_bs_.push_back((dt / 2.));
	}

	// Set time for next output
	double t_next = time + out;
	bool calc_mat = true;

	// Save starting Wavefunction
	Wavefunction Psistart = Psi;

	// Initial output
	matrices_.build(H, Psi, tree, time);
	calc_mat = false;
	Output(time, Psi, Psistart, H, tree, os);

	double err = 0.;

	while (time + 1E-6 < timeend) {
		// restrict timelength
		double dtmax = min(timeend - time, t_next - time);
		if (dt > dtmax) {
			dt = dtmax;
		} else {
			if (2 * dt > dtmax) {
				dt = min(dt, dtmax * 0.6);
			}
		}

		// Stepsize output
		cout << "time= " << time << " dt=" << dt << ", err_cmf = " << err << endl;

		// t=0 based calculation
		Wavefunction Psi0 = Psi;

		t1 = high_resolution_clock::now();
		if (calc_mat) {
			matrices_.build(H, Psi, tree, time);
			calc_mat = false;
		}

		t2 = high_resolution_clock::now();
		mattime += duration_cast<microseconds>(t2 - t1);

		double timestep = dt * tconst_;
		t1 = high_resolution_clock::now();
		CMFstep(Psi, time, time + timestep, accuracy_BS, tree);
		t2 = high_resolution_clock::now();
		steptime += duration_cast<microseconds>(t2 - t1);

		Wavefunction Psi2 = Psi;

		t1 = high_resolution_clock::now();
		CMFstep(Psi2, time + timestep, time + 2 * timestep, accuracy_BS, tree);
		t2 = high_resolution_clock::now();
		steptime += duration_cast<microseconds>(t2 - t1);

		// t=dt/2 based calculation
		t1 = high_resolution_clock::now();
		matrices_.build(H, Psi2, tree, time);
		t2 = high_resolution_clock::now();
		mattime += duration_cast<microseconds>(t2 - t1);

		t1 = high_resolution_clock::now();
		CMFstep(Psi, time + dt * tconst_, time + dt * (1 - tconst_), accuracy_BS, tree);
		t2 = high_resolution_clock::now();
		steptime += duration_cast<microseconds>(t2 - t1);

		// Propagate from t -> t + dt in one step
		Psi2 = Psi0;
		CMFstep(Psi2, time, time + dt, accuracy_BS, tree);

		// t=dt based propagation
		t1 = high_resolution_clock::now();
		matrices_.build(H, Psi2, tree, time);
		t2 = high_resolution_clock::now();
		mattime += duration_cast<microseconds>(t2 - t1);

		t1 = high_resolution_clock::now();
		CMFstep(Psi, time + dt * (1 - tconst_), time + dt, accuracy_BS, tree);
		t2 = high_resolution_clock::now();
		steptime += duration_cast<microseconds>(t2 - t1);

		// Check if step is accepted
		err = Error(Psi, Psi2, matrices_.rho_, tree);
		err *= 0.25;

		// Step refused
		if (err > accuracy_CMF) {
			// Decrease the timestep and restart
			Psi = Psi0;
			calc_mat = true;
			dt *= pow(accuracy_CMF / err, 1. / 3.) * 0.8;
		} else {
			// step accepted
			time += dt;

			// call orthogonal
			orthogonal(Psi, tree);

			// Datout
			if (time + 1E-10 >= t_next) {
				Psi.write(job.ofname, job.append_);
				Output(time, Psi, Psistart, H, tree, os);
				double time_tot = mattime.count() + steptime.count();
				cout << "Time in CMF-Integrator: " << time_tot / 1000000. << " s\n";
				cout << "Matrix calculation:     " << mattime.count() / time_tot * 100. << " %\n";
				cout << "Integration:            " << steptime.count() / time_tot * 100. << " %\n" << endl;
				t_next = time + out;
			}

			// Adjust the timestep
			dt = dt * min(pow(accuracy_CMF / err, 1. / 3.), max_increase_);
		}
	}

	//
	double time_tot = mattime.count() + steptime.count();
	cout << "Time in CMF-Integrator: " << time_tot / 1000000. << " s\n";
	cout << "Matrix calculation:     " << mattime.count() / time_tot * 100. << " %\n";
	cout << "Integration:            " << steptime.count() / time_tot * 100. << " %\n" << endl;

	job.dt = dt;
}

void CMFIntegrator::CMFstep(Wavefunction& Psi, double time, double timeend,
	double accuracy_leaf, const Tree& tree) {
	//Derivative for standard mctdh nodes
	function<void(LayerInterface&, const double, Tensorcd&, const Tensorcd&)> ddt =
		&LayerInterface::Derivative;

	//Error function for standard MCTDH nodes
	function<double(const LayerInterface&, const Tensorcd&, const Tensorcd&)> Delta =
		&LayerInterface::Error;

	// Integrate on every layer_ with constant matrices
	for (const Node& node : tree) {
//		if (node.isToplayer() || node.isBottomlayer()) {
		if (node.isToplayer() || eom_spf) {
			LayerInterface I = interfaces_[node.address()];
			bs_integrator& layer_bs = bs_integrators_[node.address()];
			Tensorcd& Phi = Psi[node];
			double layertime = time;
			layer_bs.Integrate(Phi, layertime, timeend, dt_bs_[node.address()],
				accuracy_leaf, ddt, Delta, I);
		}
	}
//	orthogonal(Psi, tree); // @TODO: this should not be here?
}

double CMFIntegrator::Error(const Wavefunction& Psi, const Wavefunction& Chi,
	const MatrixTreecd& rho, const Tree& tree) const {
	function<double(const LayerInterface&, const Tensorcd&, const Tensorcd&)>
		Delta = &LayerInterface::Error;

	/// The node-local error is not invariant under transformation
	/// between nodes. Therefore, e.g. QR-Orthogonal does not work with CMF
	/// in its current form.
	/// @OTODO: build in a non-local error measure ||Psi-Chi||=|Psi|+|Chi|-2Re(<Psi|Chi>).

	double err = 0.0;
	for (const Node& node : tree) {
		const LayerInterface& layer = interfaces_[node.address()];
		if (node.isToplayer() || eom_spf) {
			err += pow(Delta(layer, Psi[node], Chi[node]), 2);
		}
	}
	return sqrt(err);
}

void CMFIntegrator::Output(double time, const Wavefunction& Psi,
	const Wavefunction& Psistart, const Hamiltonian& H,
	const Tree& tree, ostream& os) {
	// @TODO: Check if it should be calculated here (its not efficient)
//	matrices_.build(H, Psi, tree);
/*	cdvr_file_.open("grid/cdvr." + to_string(cdvr_nr_) + ".dat");
	matrices_.cdvr_.Update(Psi, H.V_, tree, 0, true, cdvr_file_);
	cdvr_file_.close();
	cdvr_nr_++;*/

	// Output for CMF Integrator
	cout << "Time: " << time << " a.u. /" << time / 41.362 << " fs\n";

	cout << "Energies:\n";
	auto h_matrix = Expectation(matrices_, Psi, H, tree);
	auto S = TreeFunctions::dotProduct(Psi, Psi, tree);
	auto s_top = S[tree.topNode()];
	for (size_t i = 0; i < h_matrix.dim1(); ++i) {
		double e = real(h_matrix(i, i) / s_top(i, i));
		cout << i << ":\t" << e * QM::cm << " 1/cm\t" << e << " a.u." << endl;
	}

	// Calculate Autocorrelation function
	auto autocorrelation = TreeFunctions::dotProduct(Psistart, Psi, tree);
	cout << "<Psi_0(0)|Psi_0(t)>=" << autocorrelation[tree.topNode()](0, 0) << endl;

	// Calculate Datout
	TreeIO::output(Psi, tree, os);
}
