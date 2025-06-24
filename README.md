# G4Hunter-Companion

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

G4Hunter is a tool for identifying potential G-quadruplex (G4) forming sequences in DNA using sliding window scoring. This is a modernized Python 3 implementation with improved functionality and command-line interface.

The original G4Hunter tool was developed by Dr. Bedrat et al. [G4Hunter](https://github.com/AnimaTardeb/G4Hunter).

This is a rework of G4Hunter, which uses Python3 and is designed to be more user-friendly and flexible. It includes features for consensus motif analysis, and customizable parameters for G4 identification. Other improvements include better output formatting, enhanced documentation, better memory utilization by using generators and a more modular code structure. 
The core logic of G4Hunter remains the same, but the implementation has been updated to use modern Python features and libraries.

## Features

- **G4Hunter Algorithm**: Sliding window-based scoring for G-quadruplex identification
- **Consensus Motif Analysis**: Regex-based detection of G4 consensus sequences
- **Direct Repeat Detection**: Identification of direct repeat sequences
- **Flexible Parameters**: Customizable window sizes, scoring thresholds, and motif parameters
- **Multiple Output Formats**: TSV output with optional sequence merging
- **Modern Python**: Clean, well-documented Python 3.8+ code with type hints
- **Command-line Interface**: Easy-to-use CLI tool

## Installation

You can create a new environment as follows:

```
export PYTHON_VERSION=3.10
micromamba create -n g4hunter python=$PYTHON_VERSION
micromamba activate g4hunter
micromamba install -c bioconda bedtools
```

### From Source

```bash
git clone https://github.com/yourusername/g4hunter.git
cd g4hunter
pip install -e .
```

### Development Installation

```bash
git clone https://github.com/yourusername/g4hunter.git
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
    --parse-consensus \
    --merge-sequences \
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

## Requirements

- Python 3.8+
- BioPython
- NumPy
- Pandas
- pybedtools
- attrs
- BEDTools

## License

This project is licensed under the MIT License.

## How to Cite

If you use G4Hunter in your research, please cite the original publication:

```
Amina Bedrat, Laurent Lacroix, Jean-Louis Mergny, Re-evaluation of G-quadruplex propensity with G4Hunter, Nucleic Acids Research, Volume 44, Issue 4, 29 February 2016, Pages 1746–1759, https://doi.org/10.1093/nar/gkw006
```

and either of our publications:

```
Landscape and mutational dynamics of G-quadruplexes in the complete human genome and in haplotypes of diverse ancestry
Nikol Chantzi, Shia Wei Liew, Aurell Wijaya, Candace Chan, Ioannis Mouratidis, Emilyane de Oliveira Santana Amaral, Yasin Uzun, Martin Hemberg, Karen M Vasquez, Chun Kit Kwok, Ilias Georgakopoulos Soares
bioRxiv 2025.06.17.660256; doi: https://doi.org/10.1101/2025.06.17.660256
```

or

```
Quadrupia: Derivation of G-quadruplexes for organismal genomes across the tree of life
Nikol Chantzi, Akshatha Nayak, Fotis A. Baltoumas, Eleni Aplakidou, Shiau Wei Liew, Jesslyn Elvaretta Galuh, Michail Patsakis, Camille Moeckel, Ioannis Mouratidis, Saiful Arefeen Sazed, Wilfried Guiblet, Austin Montgomery, Panagiotis Karmiris-Obratański, Guliang Wang, Apostolos Zaravinos, Karen M. Vasquez, Chun Kit Kwok, Georgios A. Pavlopoulos, Ilias Georgakopoulos-Soares
bioRxiv 2024.07.09.602008; doi: https://doi.org/10.1101/2024.07.09.602008
```

## Contact

For any questions or issues, please contact Dr. Georgakopoulos-Soares or Nikol Chantzi in one of the following emails:

```
nicolechantzi@gmail.com
izg5139@psu.edu
georgakopoulos.soares@gmail.com
```