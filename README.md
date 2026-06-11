# phonon-lifetimes

Compute phonon power spectra, frequencies and lifetimes directly from molecular-dynamics (MD) atomic positions, for cubic cells — no perturbation theory required. Method and results are described in the accompanying publication: [Phys. Rev. Lett. 123, 235501 (2019)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.123.235501).

## Why

Phonon lifetimes are usually derived from perturbative force-constant expansions. This script instead extracts them non-perturbatively from the velocity/position autocorrelation of a plain MD trajectory: it Fourier-transforms the atomic motion per q-point, fits the resulting power-spectrum peaks, and reports frequencies and linewidths (inverse lifetimes) — including all anharmonic effects sampled by the MD.

## Requirements

- Python 3 (Python 2 also works — the script predates the split)
- `numpy`, `scipy`, `pandas` (and `matplotlib` for optional plotting)

```bash
pip install numpy scipy pandas matplotlib
```

## Get it

```bash
git clone https://github.com/glensk/phonon-lifetimes
cd phonon-lifetimes
```

## Usage (TiN example included)

```bash
cd example_TiN
tar -xf POSITIONs.tar.bzip2
../phonon_lifetimes.py -fftpy -atoms 216 -ps -v --peaks 2 -dt 1 -sc 3 -a 4.24 -pd -write_ps_fitted -write_full_ps -write_smooth_in
```

Key options (full list: `./phonon_lifetimes.py -h`):

| Option   | Meaning                                  |
| -------- | ----------------------------------------- |
| `-sc`    | supercell repetitions                     |
| `-atoms` | number of atoms in the cell               |
| `-dt`    | MD timestep in fs                         |
| `-a`     | lattice constant in Å                     |
| `-ps`    | compute the power spectrum                |
| `--peaks`| number of peaks to fit per q-point        |

The run writes the (smoothed/fitted) power spectra per q-point, e.g. `example_TiN/ps_smooth/ps_space_fft_1_1_1_12042.dat`:

![TiN power spectrum along [1 1 1]](example_TiN/images/ps_space_fft_1_1_1_12042.png "TiN [1 1 1]")

## Citation

If you use this code in your research, please cite:

A. Glensk et al., *Phys. Rev. Lett.* **123**, 235501 (2019). [doi:10.1103/PhysRevLett.123.235501](https://doi.org/10.1103/PhysRevLett.123.235501)

## License

[MIT](LICENSE)
