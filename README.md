<p align="center">
  <img src="\inst\app\www\LogoPGAcmdr.png" alt="PGA-cmdr" width="150"/>
</p>

<p align="center">
  Interactive R/Shiny application for population genetics analyses<br/>
  <sub>Intertryp · Universite de Montpellier · Cirad · IRD · UCAD</sub>
</p>

PGA-cmdr covers the full workflow for multilocus, individual-based
genotype datasets (microsatellites and other codominant markers): data
import and formatting, allele frequencies, general statistics, panmixia
and subdivision testing (Weir & Cockerham 1984), genetic diversities,
linkage disequilibrium, null allele screening (FreeNA), and isolation by
distance — all computed locally, with no data leaving your machine.

## Installation

All dependencies install automatically:

```r
install.packages("remotes")

remotes::install_github("KadNaf/PGA1")

pgacmdr::run_app()
```

### macOS prerequisites (compile from source)

PGA-cmdr contains C++ code compiled with OpenMP. Apple clang does **not**
include OpenMP or gfortran by default. Install these **before** running
`remotes::install_github()`:

**1. gfortran** — download from <https://mac.r-project.org/tools/> and
install `gfortran-14.2-universal.pkg` (or the current version).

**2. libomp** — OpenMP runtime:

```bash
brew install libomp
```

**3. Apple Silicon only** — add to `~/.R/Makevars` (create the file if
absent):

```
LDFLAGS += -L/opt/homebrew/opt/libomp/lib -lomp
CPPFLAGS += -I/opt/homebrew/opt/libomp/include -Xclang -fopenmp
```

## Docker

```bash
git clone https://github.com/KadNaf/PGA1.git
cd PGA1
docker compose up
# open http://localhost:3838
```

## Modules

Every module follows the same pattern: set parameters, click **Run**, and
a single `.txt` file (Allele Freq, Local/Global Panmixia, Subdivision, Diversities, LD) or a `.zip` of several files (General Stats, Null Alleles, IBD, Mantel) downloads
automatically with every result file needed — nothing is shown on screen
first. Column numbering, loci detection (single- or paired-column allele
format), and output file names are handled automatically wherever possible.

| Module | What it does |
|---|---|
| **Import Data** | CSV/TXT import (auto-detected comma/semicolon/tab separator, loaded as soon as a file is chosen). Population, Latitude/Longitude column assignment; loci range entered as "number of loci" + "first locus column" (auto-suggested, and automatically adjusted for single- or paired-column allele encoding); a numbered column-reference table is shown so column numbers are never ambiguous. |
| **Allele Freq** | Per-locus, per-population allele frequencies (plus a global column weighted by sample size), for every population and marker; a single `.txt` file. |
| **General Stats** | Per-locus Ho, Hs, Ht and FIS/FST/FIT (Weir & Cockerham 1984), optional Fst-max and Fst' (Meirmans), GST and GST''; per-population Ho, Hs and FIS with Ho/Hs averaged over loci **weighted by the number of individuals typed at each locus (Genepop convention)**; per-locus Ho of every population; per-allele F-statistics. Hs is the unbiased Nei & Chesser (1983) estimator. Output: a `.zip` of 4 files. |
| **Local Panmixia** | FIS (Weir & Cockerham 1984) by locus and by population in a **single `.txt` file** that states each type of bootstrap and its number of replicates: over individuals (`_indiv`), over sub-samples (`_subs`, by locus) and over loci (`_loci`, Overall row and each population), with two-sided permutation p-values. A second output gives FIS per locus and population with the one-sided heterozygote deficit/excess p-values (FSTAT convention). Fewer than 5 sub-samples or individuals gives NA. |
| **Global Panmixia** | Multilocus FIT across all populations in a single `.txt` file: per-locus bootstrap CI over sub-samples, and two Overall rows, `Overall (BS/Ss)` (sub-samples) and `Overall (BS/Loci)` (loci); two-sided permutation p-values. |
| **Subdivision** | A single `.txt` file with (1) FST per locus (bootstrap over sub-samples, one-sided permutation p-value based on the G statistic) and two Overall rows (over sub-samples, over loci); (2) FIS, FST, FIT, HS, HT multilocus with bootstrap CI over loci and over sub-samples; (3) the FSTAT-style G-based test of genotypic differentiation (`p_>=`, `p_>`); optionally the bootstrap replicate values. |
| **Diversities** | One file, `Diversities_confidence_intervals_<date>.txt`: HS per locus (bootstrap over individuals and over sub-samples), HS per population (over individuals and over loci), HT per locus, and the multilocus estimators (over loci and over sub-samples); Nei & Chesser (1983) estimators. Intervals from the bootstrap over individuals are shifted by the known downward bias of resampling individuals. Optional bootstrap replicate values. |
| **LD** | Pairwise linkage disequilibrium between every locus pair (loci, pairs and populations in the order of the data file), per population and combined ("All"), via a permutation-based G-test; a single `.txt` file. |
| **Null Alleles** | Null allele frequency estimation per locus × sub-sample (FreeNA EM algorithm, Chapuis & Estoup 2007), raw and ENA-corrected FST and DCSE with bootstrap CI (over loci and over sub-samples), and a full pairwise table (FST, FST-ENA, DCSE, DCSE-INA, linearised FST, geographic distance). |
| **IBD** | Isolation by distance: Rousset's (1997) regression of linearised genetic distance against geographic distance, and an independent Mantel permutation test (Pearson, Spearman, and Rousset 1D/2D statistics). |

## Getting started

See the package vignette for a full worked example on the bundled
*Boophilus* tick dataset:

```r
vignette("pgacmdr", package = "pgacmdr")
```

## Statistical methods

F-statistics follow the unbiased moment estimators of **Weir & Cockerham
(1984)**. Gene diversity (Hs, Ht) follows the unbiased estimator of **Nei
& Chesser (1983)**. Confidence intervals are obtained by non-parametric
bootstrap and p-values by Monte Carlo permutation (10,000 replicates by
default in every module, except 5,000 bootstrap replicates in Null Alleles), parallelised
in C++ via Rcpp and OpenMP. See the in-app **Help** tab for full
references.

## Citation

> PGA-cmdr: an interactive Shiny application for population genetics
> data import, exploration, and descriptive analyses. Intertryp / Universite de Montpellier /
> Cirad / IRD/ UCAD.

## Credits

**Programming:** Naffiou Kadiri and Vincent Manzanilla

**Conception:** Thierry de Meeûs

## License

MIT — see [LICENSE.md](LICENSE.md).

## Bugs & support

<https://github.com/KadNaf/PGA1/issues>
