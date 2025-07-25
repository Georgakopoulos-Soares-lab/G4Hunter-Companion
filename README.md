# G4Hunter-Companion

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This is a reworked version of G4Hunter, originally developed by Dr. Bedrat et al., licensed under GPLv3 [G4Hunter](https://github.com/AnimaTardeb/G4Hunter).

G4Hunter is a tool for identifying potential G-quadruplex (G4) forming sequences in DNA using sliding window scoring. This is a rework of G4Hunter, which uses Python3 and is designed to be more user-friendly and flexible. It includes features for consensus motif analysis, and customizable parameters for G4 identification. Other improvements include better output formatting, enhanced documentation, better memory utilization by using generators and a more modular code structure. The core logic of G4Hunter remains the same, but the implementation has been updated to use modern Python features and libraries.

This repository is developed as part of research conducted by the Georgakopoulos-Soares lab. The official page can be found at [G4Hunter-Companion](https://github.com/Georgakopoulos-Soares-lab/G4Hunter-Companion).

## Features

- **G4Hunter Algorithm**: Sliding window-based scoring for G-quadruplex identification
- **Consensus Motif Analysis**: Regex-based detection of G4 consensus sequences
- **Flexible Parameters**: Customizable window sizes, scoring thresholds, and motif parameters
- **Multiple Output Formats**: TSV output with optional sequence merging
- **Modern Python**: Clean, well-documented Python 3.8+ code with type hints
- **Command-line Interface**: Easy-to-use CLI tool

## Installation

The easiest way to install G4Hunter-Companion is through [Anaconda](https://www.anaconda.com/):

```
export PYTHON_VERSION=3.10
conda create -n g4hunter python=$PYTHON_VERSION
conda activate g4hunter
conda install -c bioconda bedtools
```

### From Source

```bash
git clone git@github.com:Georgakopoulos-Soares-lab/G4Hunter-Companion.git
cd g4hunter
pip install -e .
```

### Development Installation

```bash
git clone git@github.com:Georgakopoulos-Soares-lab/G4Hunter-Companion.git
cd g4hunter
pip install -e ".[dev]"
```

## Quick Start

### Command Line Usage

```bash
# Basic usage
export INPUT_FASTA=data/GCA_000005845.2_ASM584v2_genomic.fna.gz
g4hunter $INPUT_FASTA
```

The default parameters are window size of 25 and a minimum score of 1.5. The output will be saved in the `G4_results` directory.

```bash
# Custom parameters
g4hunter $INPUT_FASTA --window-size 25 --min-score 1.5 --outdir G4_results
```

You can also specify additional parameters for consensus motif analysis and merging sequences:

```bash
# With consensus motif analysis
g4hunter $INPUT_FASTA --parse-consensus 1 --merge-sequences
```

Or use the following command to specify all parameters:

```bash
# Advanced usage with custom parameters
g4hunter $INPUT_FASTA \
    --window-size 25 \
    --min-score 1.5 \
    --min-grun-length 3 \
    --max-loop-length 7 \
    --parse-consensus 1 \
    --merge-sequences 1 \
    --outdir my_results
```

### Python API Usage

```python
from g4hunter import G4Hunter

# Initialize G4Hunter
hunter = G4Hunter(
    outdir="G4_results",
    window_size=25,
    min_score=1.5,
    min_grun_length=3,
    max_loop_length=7
)

# Run analysis
scores = hunter.run(
    infile="input.fasta",
    parse_consensus=True,
    merge_sequences=True
)

print(f"Found {len(scores)} G4 regions")
```

## Parameters

### G4Hunter Parameters

- `--window-size`: Sliding window size for scoring (default: 25)
- `--min-score`: Minimum score threshold for G4 detection (default: 1.5)

### Consensus Motif Parameters

- `--min-grun-length`: Minimum G-run length (default: 3)
- `--min-loop-length`: Minimum loop length (default: 1)
- `--max-loop-length`: Maximum loop length (default: 7)
- `--regex-multiplicity`: Number of G-run repeats required (default: 3)

### Direct Repeat Parameters

G4Hunter-Companion also supports direct repeat extraction. This was required as a part of one of our projects but it is only a side utility and not the main focus of G4Hunter-Companion tool.

- `--min-arm-length`: Minimum arm length for direct repeats (default: 10)
- `--min-spacer-length`: Minimum spacer length (default: 0)
- `--max-spacer-length`: Maximum spacer length (default: 8)

### Analysis Options

- `--parse-consensus`: Enable consensus motif analysis
- `--merge-sequences`: Merge overlapping G4 regions
- `--outdir`: Output directory (default: G4_results)

## Output Files

G4Hunter generates several output files:

1. **`{accession}_pG4s.g4_hunter.tsv`**: Main G4Hunter results
2. **`{accession}_pG4s.consensus.tsv`**: Consensus motif results (if `--parse-consensus`)
3. **`{accession}_pG4s.merged.g4_hunter.tsv`**: Merged overlapping regions (if `--merge-sequences`)

### Output Format

Each output file contains the following columns:

- `seqID`: Sequence identifier
- `start`: Start position (0-based)
- `end`: End position (0-based)
- `sequence`: G4 sequence
- `length`: Sequence length
- `score`: G4Hunter score
- `strand`: Strand (+/-)
- `type`: Type (pG4)
- `method`: Detection method (G4Hunter, Consensus Motif, etc.)

### Consensus Motif Detection

Uses regex patterns to identify canonical G4 motifs:

- **Watson strand**: `G{n,}N{1-7}G{n,}N{1-7}G{n,}N{1-7}G{n,}`
- **Crick strand**: `C{n,}N{1-7}C{n,}N{1-7}C{n,}N{1-7}C{n,}`

Where `n` is the minimum G-run length and `N{1-7}` represents loops of 1-7 nucleotides. All the parameters can be customized via command-line options.


### Snakemake Workflow

Additionally, we provide a Snakemake workflow for running G4Hunter on multiple sequences in parallel. The workflow is defined in the `Snakefile` and can be executed with:

```bash
export CORES=4
bash submit_snake.sh $CORES
```

To use the Snakemake workflow, you need to create a `design.csv` file using `create_design.py` script.

Make sure to customize the `config.yaml` and `cluster.yaml` files to suit your environment and input data.

## Requirements

- Python 3.8+
- BioPython
- NumPy
- Pandas
- pybedtools
- attrs
- BEDTools

## License

This project is licensed under the GPL-3.0 License, consistent with the original G4Hunter implementation.

## How to Cite

If you use G4Hunter-Companion in your research, please cite the original work:

### Original G4Hunter Algorithm

**Bedrat, A., Lacroix, L., & Mergny, J. L.** (2016). Re-evaluation of G-quadruplex propensity with G4Hunter. *Nucleic Acids Research*, 44(4), 1746-1759. [https://doi.org/10.1093/nar/gkw006](https://doi.org/10.1093/nar/gkw006)

### G4Hunter-Companion Publications

You can also cite the following publications related to G4Hunter-Companion:

**Chantzi, N., Liew, S. W., Wijaya, A., Chan, C., Mouratidis, I., Amaral, E. O. S., ... & Georgakopoulos-Soares, I.** (2025). Landscape and mutational dynamics of G-quadruplexes in the complete human genome and in haplotypes of diverse ancestry. *bioRxiv*, 2025.06.17.660256. [https://doi.org/10.1101/2025.06.17.660256](https://doi.org/10.1101/2025.06.17.660256)

or

**Chantzi, N., Nayak, A., Baltoumas, F. A., Aplakidou, E., Liew, S. W., Galuh, J. E., Patsakis, M., Moeckel, C., Mouratidis, I., Sazed, S. A., Guiblet, W., Montgomery, A., Karmiris-Obratański, P., Wang, G., Zaravinos, A., Vasquez, K. M., Kwok, C. K., Pavlopoulos, G. A., & Georgakopoulos-Soares, I.** (2024). Quadrupia: Derivation of G-quadruplexes for organismal genomes across the tree of life. *bioRxiv*. https://doi.org/10.1101/2024.07.09.602008

### BibTeX Format

```bibtex
@article{bedrat2016g4hunter,
  title={Re-evaluation of G-quadruplex propensity with G4Hunter},
  author={Bedrat, Amina and Lacroix, Laurent and Mergny, Jean-Louis},
  journal={Nucleic acids research},
  volume={44},
  number={4},
  pages={1746--1759},
  year={2016},
  publisher={Oxford University Press},
  doi={10.1093/nar/gkw006}
}

@article{chantzi2025landscape,
  title={Landscape and mutational dynamics of G-quadruplexes in the complete human genome and in haplotypes of diverse ancestry},
  author={Chantzi, Nikol and Liew, Shia Wei and Wijaya, Aurell and Chan, Candace and Mouratidis, Ioannis and Amaral, Emilyane de Oliveira Santana and Uzun, Yasin and Hemberg, Martin and Vasquez, Karen M and Kwok, Chun Kit and others},
  journal={bioRxiv},
  pages={2025--06},
  year={2025},
  publisher={Cold Spring Harbor Laboratory},
  doi={10.1101/2025.06.17.660256}
}

@article{chantzi2024quadrupia,
  title={Quadrupia: Derivation of G-quadruplexes for organismal genomes across the tree of life},
  author={Chantzi, Nikol and Nayak, Akshatha and Baltoumas, Fotis A and Aplakidou, Eleni and Liew, Shiau Wei and Galuh, Jesslyn Elvaretta and Patsakis, Michail and Moeckel, Camille and Mouratidis, Ioannis and Sazed, Saiful Arefeen and others},
  journal={bioRxiv},
  pages={2024--07},
  year={2024},
  publisher={Cold Spring Harbor Laboratory},
  doi={10.1101/2024.07.09.602008}
}
```

## Contact

For any questions or issues, please contact Dr. Georgakopoulos-Soares or Nikol Chantzi at one of the following emails:

```
nicolechantzi@gmail.com
izg5139@psu.edu
georgakopoulos.soares@gmail.com
```
