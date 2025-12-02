#include <sstream>
#include <string>
#include <vector>
#include <utility>

#include <fmt/format.h>

#include <mcpp_heisenberg/logging.hpp>
#include <mcpp_heisenberg/geometry.hpp>

std::string mch::Geometry::to_poscar() const {
  std::stringstream ss;

  // title and multiplication factor
  ss << _title << "\n1.0\n";

  // lattice
  for (uint64_t il = 0; il < 3; ++il) {
    ss << fmt::format(
        " {: .8f} {: .8f} {: .8f}\n",
        _lattice_vectors.at(il, 0), _lattice_vectors.at(il, 1), _lattice_vectors.at(il, 2));
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
      ss << fmt::format(
          " {: .8f} {: .8f} {: .8f} {}\n",
          _positions.at(ni + ii, 0),
          _positions.at(ni + ii, 1),
          _positions.at(ni + ii, 2),
          it.first);
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

mch::Geometry mch::Geometry::from_poscar(std::shared_ptr<std::istream> istream) {
  // read the first line as is
  std::string title;

  char c = istream->get();

  while (c != 0 && c != '\n') {
    title += c;
    c = istream->get();
  }

  LOGD << "title is `" << title << "`";

  // then create lexer
  auto lexer = WordLexer(istream);

  // factor
  double factor = lexer.next_double();
  LOGD << "factor is `" << factor << "`";

  auto is_EOL = [](std::string f) { return f == "\n"; };
  lexer.skip_if(is_EOL, "expected EOL");

  // lattice
  arma::mat lattice(3, 3);

  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
        lattice.at(i, j) = factor * lexer.next_double();
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

  std::vector<mch::ion_type_t> ions;
  uint64_t natoms = 0;
  for (uint64_t in = 0; in < ion_types.size(); ++in) {
    ions.push_back(std::make_pair(ion_types[in], static_cast<uint64_t>(lexer.next_integer())));
    natoms += ions.at(in).second;
  }

  LOGD << "ions are " << ions;

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
      positions.at(iatm, j) = lexer.next_double();
    }

    if (iatm < natoms - 1) {  // skip whatever is remaining on that line
      next_word = lexer.next_word();
      while (next_word != "\n") {
        next_word = lexer.next_word();
      }
    }
  }

  if (!is_direct) {
    LOGD << "convert to fractional coordinates";

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

  return Geometry(title, lattice, ions, positions);
}

