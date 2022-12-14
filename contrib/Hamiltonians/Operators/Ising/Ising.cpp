#include "Ising.h"
#include <random>
#include "Pauli.h"


Ising::Ising(const mctdhBasis& basis, size_t N, size_t method) {
//	HamiltonianPaths(basis, N, method);
}

void Ising::SpinChain(const mctdhBasis& basis, double J, double H) {
    using namespace PauliMatrices;
    size_t N = basis.nLeaves();

    for (size_t i = 0; i < N; ++i) {
        MPO M;
		M.push_back(sigma_z, i);
        push_back(M, -H);
    }

    for (size_t i = 0; i < N; ++i) {
        MPO M;
        size_t j = (i + 1) % N;
		M.push_back(sigma_x, i);
		M.push_back(sigma_x, j);
        push_back(M, -J);
    }
}

Graph make_random_graph(size_t N, size_t seed) {
	// Initialize randomizer
//	using namespace std::chrono;
//	auto seed = system_clock::now().time_since_epoch().count();
	cout << "Seed used for graph generation: " << seed << endl;
	mt19937 gen(seed);

	// Build a (random) Graph
	Graph graph;
	for (size_t i = 0; i < N; ++i) {
		graph.push_back(Vertex(gen));
	}
	graph.all_distance_edges();

	return graph;
}

void Ising::HamiltonianPaths(const mctdhBasis& basis, size_t N, size_t seed) {
	using namespace PauliMatrices;

	size_t Nnod = N + 1;
	Graph graph = make_random_graph(Nnod, seed);
	cout << "randomly generated graph:\n";
	graph.print();
	size_t count_a = 0;
	size_t count_b = 0;
	size_t count_c = 0;
	double A = 1.;

	// Penalty term H_B
	for (size_t i = 0; i < Nnod; ++i) {
		shared_ptr<qbit_x> x_i(new qbit_x(i));
		for (size_t nu = 0; nu < Nnod; ++nu) {
			for (size_t mu = 0; mu < Nnod; ++mu) {

				MultiParticleOperator M;
//				if (nu < N) {
					M.push_back(x_i, nu);
//				}
//				if (mu < N) {
					M.push_back(x_i, mu);
//				}
//				mpos.push_back(M);
//				coeff.push_back(A);
				push_back(M, A);
				count_b++;
			}
		}
	}

	// -N part of H_B
//	{
//		MultiParticleOperator M;
//		mpos.push_back(M);
//		coeff.push_back(-Nnod);
//		count_b++;
//		
//	}
	{
		MultiParticleOperator M;
		shared_ptr<qbit_x> x_i(new qbit_x(N));
		M.push_back(x_i, N);
//		mpos.push_back(M);
//		coeff.push_back(-A);
		push_back(M, -A);
	}

	// Weight term
	for (size_t i = 0; i < graph.nEdges(); ++i) {
		Edge e = graph.edge(i);
		size_t u = e.src;
		size_t v = e.dest;
		double W = e.weight;
		for (size_t j = 0; j < N; ++j) {
			size_t jp = (j + 1) % Nnod;
			shared_ptr<qbit_x> x_j(new qbit_x(j));
			shared_ptr<qbit_x> x_jp(new qbit_x(jp));

			MultiParticleOperator M;
//			if ((u < N ) && (j < N)) {
				M.push_back(x_j, u);
//			}
//			if ((v < N) && (jp < N)) {
				M.push_back(x_jp, v);
//			}
//			mpos.push_back(M);
//			coeff.push_back(W);
			push_back(M, W);
			count_c++;
		}
	}

	size_t count = count_a + count_b + count_c;
	cout << "Fed " << count << " parts into the hamiltonian.\n";
	cout << "count in H_A = " << count_a << "\n";
	cout << "count in H_B = " << count_b << "\n";
	cout << "count in H_C = " << count_c << "\n";
}

/*void Ising::HamiltonianPaths(const mctdhBasis& basis, size_t N, size_t method) {
	using namespace PauliMatrices;

	Graph graph = make_random_graph(N + 1);
	cout << "randomly generated graph:\n";
	graph.print();
	size_t count_a = 0;
	size_t count_b = 0;
	size_t count_c = 0;

	// Penalty term H_A
	double A = 1.0;
	for (size_t nu = 0; nu < graph.nVertices() - 1; ++nu) {
		for (size_t i = 0; i < graph.nVertices() - 1; ++i) {
			for (size_t j = 0; j < graph.nVertices() - 1; ++j) {
				MultiParticleOperator M;
				M.push_back(bit_x, idx(nu, i, N));
				M.push_back(bit_x, idx(nu, j, N));
				mpos.push_back(M);
				coeff.push_back(A);
				count_a++;
			}
		}
	}

	// Penalty term H_B
	for (size_t i = 0; i < graph.nVertices() - 1; ++i) {
		for (size_t nu = 0; nu < graph.nVertices() - 1; ++nu) {
			for (size_t mu = 0; mu < graph.nVertices() - 1; ++mu) {
				MultiParticleOperator M;
				M.push_back(bit_x, idx(nu, i, N));
				M.push_back(bit_x, idx(mu, i, N));
				mpos.push_back(M);
				coeff.push_back(A);
				count_b++;
			}
		}
	}

	// Weight term
	for (size_t i = 0; i < graph.nEdges(); ++i) {
		Edge e = graph.edge(i);
		size_t u = e.src;
		size_t v = e.dest;
		double W = e.weight;
		for (size_t j = 0; j < N - 1; ++j) {
			size_t jp = (j + 1) % N;
			MultiParticleOperator M;
			if ((u < (N - 1)) && (j < (N - 1))) {
				M.push_back(bit_x, idx(u, j, N));
			}
			if ((v < (N - 1)) && (jp < (N - 1))) {
				M.push_back(bit_x, idx(v, jp, N));
			}
			mpos.push_back(M);
			coeff.push_back(W);
			count_c++;
		}
	}

	size_t count = count_a + count_b + count_c;
	cout << "Fed " << count << " parts into the hamiltonian.\n";
	cout << "count in H_A = " << count_a << "\n";
	cout << "count in H_B = " << count_b << "\n";
	cout << "count in H_C = " << count_c << "\n";

}
*/


