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

Where `POSCAR` is a valid _geometry in the [VASP POSCAR format](https://www.vasp.at/wiki/index.php/POSCAR), and `simulation.toml` contains the parameters:

```toml
# system
supercell = [5, 5, 1]  # number of repetitions for each side
magnetic_sites = ['H']  # list of magnetic atoms

# definition of coupling between two atoms:
pair_defs = [
  # 'atom1', 'atom2', distance (in Å), J value
  ['H', 'H', 2.0, 1.0]
]

# for each spin, you can change the default, which is 1
spin_values = { H = 0.5 }

# for each spin, you can change the default, which is 0
magnetic_anisotropies = { H = 0.1 } 

# for each spin, you can change the default, which is "allup"
initial_configs = { H = 'random' }

# simulation
kB = 1.0  # boltzmann constant in unit of J/T
muB = 1.0  # bohr magneton, in units of J/H
T = 0.1  # Temperature
H = 0   # magnetic field
N = 20000  # number of MC steps

# save data frame
save_interval = 1000 # saving interval, thus requireing `save_interval*n_mag_sites*sizeof(double)` of cache memory
deflate_level = 6  # compression level for the `results/configs` dataset
chunk_size = 1024  # chunk size for the `results/configs` dataset
```

The value for `kB` and `muB` will depend on the unit for temperature (`T`), magnetic field (`H`), and energy (`J`).
In particular, `kB = 8.61733326e-5` and `muB = 5.788381e-5` are relevant when using Kelvin, Tesla, and eV, while `kB = 1` and `muB = 1` in reduced units.

## Results

The energy for a given configuration $\{S_i\}$ is computed as:

$E = -\sum_{(i,j)} J_{ij} S_i S_j - \mu_B H \sum_i S_i - \sum_i \Delta_i S_i^2,$

where $\sum_{(i,j)}$ runs over all pairs $(i, j)$, $J_{ij}$ is the coupling interaction, and is the external magnetic field (applied in the $z$ direction) and $\Delta$ is the magnetic anisotropy.

A heat bath Monte-Carlo is performed, based on $P=\max\left(1, \exp\left[-\frac{\Delta E_i}{k_B T}\right]\right)$, where $\Delta E_i$ is the change in energy due to flipping spin $i$.

### H5 file

After a run, the H5 file contains the following datasets:

Group `geometry/`:

+ `geometry/lattice_vectors` (`double (3, 3)`): lattice vectors;
+ `geometry/ion_types` (`uint64_t (Ntypes, )`): type (chemical symbol) of each magnetic ion/spin;
+ `geometry/ion_numbers` (`str (Ntypes, )`): number of each magnetic ion/spin type;
+ `geometry/positions` (`uint64_t (3, Nspins)`): position of each magnetic ion/spin;
+ `geometry/spin_values` (`double (Npsins, )`): for each magnetic/ion type, the initial spin value ($|S_i|$) ;
+ `geometry/magnetic_anisotropies` (`double (Npsins, )`): for each magnetic/ion type, the magnetic anisotropy ($\Delta_i$) ;

Group `hamiltonian/`:

+ `hamiltonian/pairs` (`uint64_t (Npairs, 2)`): list of pairs;
+ `hamiltonian/J` (`double (Npairs, )`): for each pair, the value of the magnetic coupling, $J_{ij}$ ;
+ `hamiltonian/parameters` (`double (4,)`): value of the [Boltzmann constant](https://en.wikipedia.org/wiki/Boltzmann_constant) ($k_B$), of $T$, of the [Bohr magneton](https://en.wikipedia.org/wiki/Bohr_magneton) ($\mu_B$), and of $H$ used during simulation;

Group `results/`:

+ `results/configs` (`int8_t (Nsteps, Nspins)`): for each step, `sign(spins)`;
+ `results/aggregated_data` (`double (Nsteps, 2)`): for each step, the (Ising) energy (col 0) and the sum of spins (col 1).