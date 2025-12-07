#include <sstream>
#include <string>
#include <vector>
#include <utility>
#include <algorithm>

#include <mcpp_heisenberg/logging.hpp>
#include <mcpp_heisenberg/geometry.hpp>

namespace mch {

std::string Geometry::to_poscar() const {
  std::stringstream ss;

  // title and multiplication factor
  ss << _title << "\n1.0\n";

  // lattice
  for (uint64_t il = 0; il < 3; ++il) {
    ss << fmt::format(" {: .8f} {: .8f} {: .8f}\n", _lattice_vectors.at(il, 0), _lattice_vectors.at(il, 1),
                      _lattice_vectors.at(il, 2));
  }

  // ions
  for (auto& it : _ions) {
    ss << " " << it.first;
  }

  ss << "\n";

  for (auto& it : _ions) {
    ss << " " << it.second;
  }

  ss << "\n";

  // positions
  ss << "Direct\n";

  uint64_t ni = 0;
  for (auto& it : _ions) {
    for (uint64_t ii = 0; ii < it.second; ++ii) {
      ss << fmt::format(" {: .8f} {: .8f} {: .8f} {}\n", _positions.at(ni + ii, 0), _positions.at(ni + ii, 1),
                        _positions.at(ni + ii, 2), it.first);
    }

    ni += it.second;
  }

  return ss.str();
}

void _skip_space(std::istream& istream) {
  char c = istream.get();

  while (c != 0 && std::isspace(c)) {
    c = istream.get();
  }
}

std::string _get_word(std::istream& istream) {
  std::stringstream ss;

  _skip_space(istream);

  char c = istream.get();

  if (c == 0) {
    throw std::runtime_error("unexpected end of file");
  }

  while (!std::isspace(c)) {
    ss << c;
    c = istream.get();
  }

  return ss.str();
}

Geometry Geometry::from_poscar(std::shared_ptr<std::istream> istream) {
  elapsed::Chrono chrono;
  LOGI << "*> Read geometry from POSCAR";

  // read the first line as is
  std::string title;

  char c = istream->get();

  while (c != 0 && c != '\n') {
    title += c;
    c = istream->get();
  }

  LOGV << "title is `" << title << "`";

  // then create lexer
  auto lexer = WordLexer(istream);

  // factor
  double factor = lexer.next_double_or_throw();
  LOGV << "factor is `" << factor << "`";

  auto is_EOL = [](std::string f) { return f == "\n"; };
  lexer.skip_if(is_EOL, "expected EOL");

  // lattice
  arma::mat lattice(3, 3);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      lattice.at(i, j) = factor * lexer.next_double_or_throw();
    }

    lexer.skip_if(is_EOL, "expected EOL");
  }

  // ions
  std::string next_word;
  std::vector<std::string> ion_types;
  next_word = lexer.next_word();

  while (next_word != "\0" && next_word != "\n") {
    ion_types.push_back(next_word);
    next_word = lexer.next_word();
  }

  std::vector<ion_type_t> ions;
  uint64_t natoms = 0;
  for (uint64_t in = 0; in < ion_types.size(); ++in) {
    ions.push_back(std::make_pair(ion_types[in], static_cast<uint64_t>(lexer.next_integer_or_throw())));
    natoms += ions.at(in).second;
  }

  LOGV << "ions are " << ions;

  lexer.skip_if(is_EOL, "expected EOL");

  // type
  bool is_direct = false;
  next_word = lexer.next_word();

  if (next_word.at(0) == 'S' || next_word.at(0) == 's') {  // skip selective dynamics
    next_word = lexer.next_word();

    while (next_word != "\n") {
      next_word = lexer.next_word();
    }

    next_word = lexer.next_word();
  }

  if (next_word.at(0) == 'D' || next_word.at(0) == 'd') {
    is_direct = true;
  }

  lexer.skip_if(is_EOL, "expected EOL");

  // positions
  arma::mat positions(natoms, 3);

  for (uint64_t iatm = 0; iatm < natoms; ++iatm) {
    for (int j = 0; j < 3; ++j) {
      positions.at(iatm, j) = lexer.next_double_or_throw();
    }

    if (iatm < natoms - 1) {  // skip whatever is remaining on that line
      next_word = lexer.next_word();
      while (next_word != "\n") {
        next_word = lexer.next_word();
      }
    }
  }

  LOGD << "done reading POSCAR";

  if (!is_direct) {
    LOGV << "convert to fractional coordinates";

    arma::mat cartesian_positions(positions);
    arma::mat inv = arma::inv(lattice);

    positions.zeros();

    for (uint64_t i = 0; i < natoms; ++i) {
      for (uint64_t k = 0; k < 3; ++k) {
        for (uint64_t j = 0; j < 3; ++j) {
          positions.at(i, k) += cartesian_positions.at(i, j) * inv.at(j, k);
        }
      }
    }
  }

  LOGI << "<* Done (took " << chrono.format() << ")";

  return {title, lattice, ions, positions};
}

Geometry Geometry::to_supercell(uint64_t nx, uint64_t ny, uint64_t nz, bool sort) const {
  elapsed::Chrono chrono;
  LOGI << "*> Make a " << nx << "x" << ny << "x" << nz << " supercell";

  // make new positions and ions
  std::vector<ion_type_t> ions;
  uint64_t n = _positions.n_rows;
  uint64_t N = nx * ny * nz;

  auto positions = arma::mat(N * n, 3);

  uint64_t nc = 0;

  for (uint64_t ix = 0; ix < nx; ++ix) {
    for (uint64_t jy = 0; jy < ny; ++jy) {
      for (uint64_t kz = 0; kz < nz; ++kz) {
        auto shift = arma::vec(3);
        shift.at(0) = ix;
        shift.at(1) = jy;
        shift.at(2) = kz;

        positions.rows(nc * n, (nc + 1) * n - 1) = _positions;
        positions.rows(nc * n, (nc + 1) * n - 1).each_row() += shift.t();
        nc++;

        ions.insert(ions.end(), _ions.begin(), _ions.end());
      }
    }
  }

  positions.col(0) /= nx;
  positions.col(1) /= ny;
  positions.col(2) /= nz;

  auto diagvec = arma::vec(3);
  diagvec.at(0) = nx;
  diagvec.at(1) = ny;
  diagvec.at(2) = nz;
  arma::mat lattice = (_lattice_vectors.t() * arma::diagmat(diagvec)).t();

  if (sort) {
    std::vector<ion_type_t> sorted_ions;
    auto sorted_positions = arma::mat(N * n, 3);

    uint64_t ni = 0;
    for (auto& it : _ions) {
      uint64_t nj = 0, nions = 0;
      for (auto& jt : ions) {
        uint64_t nk = jt.second;
        if (jt.first == it.first) {
          sorted_positions.rows(ni, ni + nk - 1) = positions.rows(nj, nj + nk - 1);
          ni += nk;
          nions += nk;
        }

        nj += nk;
      }

      sorted_ions.push_back(std::make_pair(it.first, nions));
    }

    ions = sorted_ions;
    positions = sorted_positions;
  }

  LOGI << "<* Done (took " << chrono.format() << ")";

  return {_title, lattice, ions, positions};
}

Geometry Geometry::filter_atoms(const std::vector<std::string>& atoms) const {
  LOGD << "Filter atoms: " << atoms;

  // count
  uint64_t N = 0;

  for (auto& it : _ions) {
    for (auto& jt : atoms) {
      if (it.first == jt) {
        N += it.second;
      }
    }
  }

  // filter
  auto positions = arma::mat(N, 3);
  auto ions = std::vector<ion_type_t>();

  uint64_t ni = 0, nj = 0;
  for (auto& it : _ions) {
    for (auto& jt : atoms) {
      if (it.first == jt) {
        ions.push_back(it);
        positions.rows(ni, ni + it.second - 1) = _positions.rows(nj, nj + it.second - 1);
        ni += it.second;
      }
    }

    nj += it.second;
  }

  return {_title, _lattice_vectors, ions, positions};
}

void Geometry::to_h5_group(HighFive::Group& group) const {
  // Save ion number & types
  std::vector<uint64_t> ion_numbers(_ions.size());
  std::transform(_ions.cbegin(), _ions.cend(), ion_numbers.begin(), [](auto& t) { return t.second; });
  auto dset_ion_numbers = group.createDataSet("ion_numbers", ion_numbers);
  dset_ion_numbers.write(ion_numbers);

  std::vector<std::string> ion_types(_ions.size());
  std::transform(_ions.cbegin(), _ions.cend(), ion_types.begin(), [](auto& t) { return t.first; });
  auto dset_ion_types = group.createDataSet("ion_types", ion_types);
  dset_ion_types.write(ion_types);

  // save lattice
  auto dset_lattice = group.createDataSet<double>("lattice_vectors", HighFive::DataSpace({3, 3}));
  dset_lattice.write_raw(_lattice_vectors.memptr());

  // save positions
  auto dset_positions = group.createDataSet<double>("positions", HighFive::DataSpace({3, _positions.n_rows}));
  dset_positions.write_raw(_positions.memptr());
}

}  // namespace mch
