# `mcpp_heisenberg`

Yet another Monte Carlo code for the evolution of the Heisenberg/Ising hamiltonian.
See [Wikipedia on Ising model](https://en.wikipedia.org/wiki/Ising_model) for a brief introduction.

## Installation

```bash
# clone
git clone git@github.com:pierre-24/mcpp_heisenberg.git
cd mcpp_heisenberg

# compile
meson setup _build
meson compile -C _build

# test
meson test --print-errorlog -C _build

# install
meson install -C _build
```

## Usage

```bash
mcpp_heisenberg_run -g POSCAR -i simulation.toml -o output.h5
```

Where `POSCAR` is a valid geometry in the [VASP POSCAR format](https://www.vasp.at/wiki/index.php/POSCAR), and `simulation.toml` contains the parameters:

```toml
# system
supercell = [5, 5, 1]  # number of repetitions for each side
magnetic_sites = ['H']  # list of magnetic atoms

# definition of coupling between two atoms:
pair_defs = [
  # 'atom1', 'atom2', distance (in Å), J value
  ['H', 'H', 2.0, 1.0]
]

# simulation
T = 0.1  # Temperature
H = 0   # magnetic field (not implemented yet)
N = 20000  # number of MC steps
step_type = 'sweep'  # either 'sweep' (update sites one by ones) or 'cluster' (switch whole cluster) 
```