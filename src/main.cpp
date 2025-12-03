#include <mcpp_heisenberg/mc.hpp>

constexpr uint64_t MAX = 1000;
constexpr uint64_t N = 10;

int main() {
  mch::init_logging();

  uint64_t PROD_START = 0;

  arma::mat lattice = arma::eye(3, 3);
  lattice.at(0, 0) = 2;
  lattice.at(1, 1) = 2;
  lattice.at(2, 2) = 100;

  arma::mat positions = arma::mat(1, 3);

  auto geometry = mch::Geometry("square", lattice, {{"H", 1}}, positions).to_supercell(N, N, 1);
  auto square_hamiltonian = mch::Hamiltonian::from_geometry(geometry, {"H"}, {{"H", "H", 2.0, 1.0}});

  arma::vec temperature_range = arma::linspace(0.1, 3, 16);

  for (uint64_t iT = 0; iT < temperature_range.n_rows; ++iT) {
    auto runner = mch::MonteCarloRunner(square_hamiltonian);

    double T = temperature_range.at(iT);

    arma::mat stats(MAX - PROD_START, 2);

    for (uint64_t i = 0; i < MAX; ++i) {
      runner.cluster_update(T);

      if (i > PROD_START) {
        stats.row(i - PROD_START) = {runner.energy(), arma::sum(runner.spins())};
      }
    }

    LOGI << T
         << "," << arma::mean(stats.col(0)) / runner.N()
         << "," << arma::mean(stats.col(1)) / runner.N();
  }

  return EXIT_SUCCESS;
}
