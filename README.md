# Blast101 — Bioinformatics Sequence Search Tool

A Python implementation of a BLAST-like heuristic protein sequence search algorithm,
extended with a command line interface and automated test suite as part of a
Bioinformatics Algorithms course assessment.

## Overview

Blast101 implements the core components of the BLAST algorithm from first principles:
- Word-based heuristic search with neighbourhood expansion
- Smith-Waterman local alignment with BLOSUM scoring
- Gumbel distribution-based E-value and bit score calculation
- Searches against UniProt Swiss-Prot and custom FASTA databases

## Features Added

### Command Line Interface (`blast101_cli.py`)
- Single entry point replacing three separate PyCharm run configurations
- Four modes: `blast`, `sw` (Smith-Waterman), `stats`, `test`
- Input validation: DNA detection, invalid residue checking, database file verification
- FASTA file input via `--query-file` for long sequences
- Switchable BLOSUM matrix (45/50/62/80/90) and gap penalty
- Optional results file output via `--output`
- Auto-detects and switches to the correct virtual environment Python
- Settings summary printed before each run

### Automated Test Suite (`test_blast101.py`)
- 55 tests across 7 test classes using Python `unittest`
- Manually verified Smith-Waterman scores as numerical ground truth
- Mathematical invariant testing (score symmetry, E-value bounds)
- Real bug detection: module-level `bestscore` variable never reset between calls
- Clean formatted output grouped by class when run via the CLI

## Installation
```bash
# Clone the repository
git clone https://github.com/B291900-2025/Bioinformatics-Algorithms.git
cd Bioinformatics-Algorithms/Blast101_code

# Create virtual environment and install dependencies
python -m venv BA_ICA
BA_ICA/bin/pip install blosum pandas scipy seaborn matplotlib
```

## Usage
```bash
# Run BLAST search with default query and database
python blast101_cli.py --mode blast --database your_database.fasta

# Run with custom query sequence
python blast101_cli.py --mode blast \
    --query PWNAAPLHNFGEDFLQPYVQLQQNFSASDLEVNLEATRESHAHFSTPQALELFLNYSVTP \
    --database your_database.fasta

# Run with query from FASTA file
python blast101_cli.py --mode blast \
    --query-file query.fasta \
    --database your_database.fasta

# Change BLOSUM matrix and gap penalty
python blast101_cli.py --mode blast \
    --database your_database.fasta \
    --blosum 80 --gap -10

# Save results to file
python blast101_cli.py --mode blast \
    --database your_database.fasta \
    --output results.txt

# Run Smith-Waterman exhaustive search
python blast101_cli.py --mode sw --database your_database.fasta

# Run statistical analysis (Gumbel fit)
python blast101_cli.py --mode stats

# Run automated test suite
python blast101_cli.py --mode test

# Show all options
python blast101_cli.py --help
```

## Project Structure

```
Blast101_code/
├── blast101_cli.py          # Command line interface (new)
├── test_blast101.py         # Automated test suite (new)
├── blast_101_search.py      # BLAST heuristic search engine
├── smith_waterman_search.py # Smith-Waterman exhaustive search
├── smith_waterman_p.py      # SW algorithm implementation
├── process_fasta_file.py    # FASTA database parser
├── create_seq_word_dict.py  # BLAST word dictionary builder
├── calc_bit_and_evalues.py  # E-value and bit score calculations
├── build_expect_scores.py   # Gumbel distribution fitting
├── programme_settings.py    # Settings file reader/writer
├── print_logger.py          # Dual screen/file logger
└── settings.ini             # Persistent configuration
```

## Algorithm Overview

The BLAST search proceeds in three stages. First, the query sequence is broken into
overlapping words of length `word_size` (default 4), and a neighbourhood dictionary
is built containing all similar words exceeding a threshold score T against the
BLOSUM matrix. Second, the database is scanned for word matches, and matching word
pairs on the same diagonal within `max_extension_length` residues are extended to
form High-Scoring Segment Pairs (HSPs). Third, the top-scoring HSPs are re-scored
using full Smith-Waterman alignment, and E-values are calculated using a Gumbel
distribution fitted to random alignment score simulations.

## Configuration

Key parameters are stored in `settings.ini` and can be overridden at the command
line for the most important settings. Advanced parameters (word size, T-score
threshold, extension length) are managed via the settings file directly.

## Dependencies

- `blosum` — BLOSUM substitution matrices
- `pandas` — data handling for statistical analysis
- `scipy` — Gumbel distribution fitting
- `seaborn` / `matplotlib` — score distribution visualisation

## References

- Altschul et al. (1990) Basic local alignment search tool. *J Mol Biol* 215(3)
- Altschul et al. (1997) Gapped BLAST and PSI-BLAST. *Nucleic Acids Res* 25(17)
- Smith & Waterman (1981) Identification of common molecular subsequences. *J Mol Biol* 147(1)
- Henikoff & Henikoff (1992) Amino acid substitution matrices from protein blocks. *PNAS* 89(22)
