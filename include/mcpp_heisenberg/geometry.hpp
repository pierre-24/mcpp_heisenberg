#ifndef  INCLUDE_MCPP_HEISENBERG_GEOMETRY_HPP_
#define  INCLUDE_MCPP_HEISENBERG_GEOMETRY_HPP_

#include <cassert>
#include <vector>
#include <utility>
#include <string>
#include <memory>

#include <fmt/format.h>

#include <mcpp_heisenberg/arma.hpp>
#include <mcpp_heisenberg/logging.hpp>

namespace mch {

/**
 * A word lexer, yield space-separated "words" (which can also be `\0` and `\n`).
 */
class WordLexer {
 protected:
  int _c{0};
  uint64_t _line{0};
  std::shared_ptr<std::istream> _stream;

  void _skip_space() {
    while (_c != 0 && (_c == ' ' || _c == '\t')) {
      _next();
    }
  }

  void _next() { _c = _stream->get(); }

 public:
  explicit WordLexer(std::shared_ptr<std::istream>& stream): _stream(stream) {
    _next();
  }

  /// Get current line
  [[nodiscard]] uint64_t line() const { return _line; }

  /// Get next word, `\n` or `\0`.
  std::string next_word() {
    _skip_space();

    if (_stream->eof()) {
      return "\0";
    }

    if (_c == '\n') {
      _line++;

      _next();
      return "\n";
    }

    std::string word;
    while (!_stream->eof() && !std::isspace(_c)) {
      word += _c;
      _next();
    }

    return word;
  }

  /// Skip if
  void skip_if(std::function<bool(std::string)> predicate, std::string msg) {
    if (!predicate(next_word())) {
      throw std::runtime_error(fmt::format("error at line {} :: {}", line() + 1, msg));
    }
  }

  /// Get next double, yield error if not
  double next_double_or_throw() {
    double dbl = .0;

    try {
      dbl = std::stod(next_word());
    } catch (std::exception& e) {
      throw std::runtime_error(fmt::format("error at line {}: cannot convert to double", line() + 1));
    }

    return dbl;
  }

  /// Get next integer, yield error if not
  int64_t next_integer_or_throw() {
    int64_t ix = 0;

    try {
      ix = std::stol(next_word());
    } catch (std::exception& e) {
      throw std::runtime_error(fmt::format("error at line {}: cannot convert to double", line() + 1));
    }

    return ix;
  }
};

using ion_type_t = std::pair<std::string, uint64_t>;

class Geometry {
 protected:
  /// Title
  std::string _title;
  /// Lattice vectors
  arma::mat _lattice_vectors;
  /// types & number of each ions
  std::vector<ion_type_t> _ions;
  /// positions (fractional coordinates)
  arma::mat _positions;

 public:
  Geometry() = delete;

  Geometry(
      const std::string& title, const arma::mat& lattice, const std::vector<ion_type_t>& ions, const arma::mat& pos):
  _title(title), _lattice_vectors(lattice), _ions(ions), _positions(pos) {
    uint64_t N = std::accumulate(
        ions.begin(), ions.end(), 0, [](uint64_t i, const ion_type_t& j) { return i + j.second; });
    assert(N == _positions.n_rows);
  }

  /// Get number of atoms
  [[nodiscard]] uint64_t number_of_atoms() const { return _positions.n_rows; }

  /// Get lattice
  [[nodiscard]] const arma::mat& lattice() const { return _lattice_vectors; }

  /// Get ions
  [[nodiscard]] const std::vector<ion_type_t>& ions() const { return _ions; }

  /// Get positions (fractional coordinates)
  [[nodiscard]] const arma::mat& positions() const { return _positions; }

  /// format as a POSCAR file
  [[nodiscard]] std::string to_poscar() const;

  /// Create a geometry out of POSCAR file
  static Geometry from_poscar(std::shared_ptr<std::istream> istream);

  /// Get a, b, c
  [[nodiscard]] arma::vec get_abc() const {
    arma::vec abc(3);
    for (uint64_t i = 0; i < 3; ++i) {
      abc.at(i) = arma::norm(_lattice_vectors.row(i));
    }

    return abc;
  }

  /// Get alpha beta gamma
  [[nodiscard]] arma::vec get_alpha_beta_gamma() const {
    arma::vec abg(3);

    auto angle_between = [](arma::rowvec v1, arma::rowvec v2) {
      v1 /= arma::norm(v1);
      v2 /= arma::norm(v2);

      return acos(arma::dot(v1, v2));
    };

    abg.at(0) = angle_between(_lattice_vectors.row(1), _lattice_vectors.row(2));
    abg.at(1) = angle_between(_lattice_vectors.row(0), _lattice_vectors.row(2));
    abg.at(2) = angle_between(_lattice_vectors.row(0), _lattice_vectors.row(1));

    return abg;
  }

  /// Create a `nx` x `ny` x `nz` supercell.
  [[nodiscard]] Geometry to_supercell(uint64_t nx, uint64_t ny, uint64_t nz, bool sort = true) const;

  /// Get a new geometry, but only with the atoms in `atoms`
  [[nodiscard]] Geometry filter_atoms(const std::vector<std::string>& atoms) const;
};

}  // namespace mch

#endif  // INCLUDE_MCPP_HEISENBERG_GEOMETRY_HPP_
