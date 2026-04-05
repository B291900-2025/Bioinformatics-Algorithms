#!/home/s2793337/Bioinformatics_Algorithms/ICA/Blast101_code/BA_ICA/bin/python3
import sys
VENV_PYTHON = "/home/s2793337/Bioinformatics_Algorithms/ICA/Blast101_code/BA_ICA/bin/python3"
if sys.executable != VENV_PYTHON:
    import subprocess
    result = subprocess.run([VENV_PYTHON] + sys.argv, stdout=sys.stdout, stderr=sys.stderr)
    sys.exit(result.returncode)

################################################################################
# blast101_cli.py
# Command Line Interface for the Blast101 application
# Provides a single entry point for BLAST search, SW search, stats, and tests
# Validates inputs before running to catch errors early
#
# Usage:
#   python blast101_cli.py --mode blast --query <seq> --database <file>
#   python blast101_cli.py --mode blast --query <seq> --database <file> --output results.txt
#   python blast101_cli.py --mode sw --query <seq>
#   python blast101_cli.py --mode stats
#   python blast101_cli.py --mode test
#
################################################################################

import argparse
import os
import sys

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
VALID_RESIDUES = set("GAVLITSMCPFYWHKRDENQ")
DNA_ONLY       = set("ACGT")
MODES          = ["blast", "sw", "stats", "test"]

BANNER = """
╔══════════════════════════════════════════════════════════════════╗
║              Blast101 - Bioinformatics Search Tool               ║
║                         Beta Version                             ║
║   Modes: blast | sw | stats | test                               ║
╚══════════════════════════════════════════════════════════════════╝
"""

# ---------------------------------------------------------------------------
# Sequence Validation
# ---------------------------------------------------------------------------
def validate_sequence(seq):
    """
    Validates a protein query sequence.
    Checks for: empty input, DNA-only composition, invalid characters.
    Returns cleaned sequence string or raises SystemExit with a clear message.
    """
    if not seq or len(seq.strip()) == 0:
        print("[ERROR] Query sequence is empty.")
        sys.exit(1)

    seq = seq.strip().upper()

    # Check for DNA before residue check — A/C/G/T are valid amino acids
    # so DNA would otherwise pass residue validation
    if all(c in DNA_ONLY for c in seq):
        print("[ERROR] Query sequence looks like DNA (only A/C/G/T found).")
        print("        Blast101 requires a PROTEIN sequence.")
        print(f"        Sequence received: {seq[:60]}{'...' if len(seq)>60 else ''}")
        sys.exit(1)

    # Check for invalid characters
    invalid = [c for c in seq if c not in VALID_RESIDUES]
    if invalid:
        print(f"[ERROR] Invalid residues found in query sequence: {set(invalid)}")
        print(f"        Valid amino acid codes are: {''.join(sorted(VALID_RESIDUES))}")
        print(f"        Sequence received: {seq[:60]}{'...' if len(seq)>60 else ''}")
        sys.exit(1)

    # Warn if sequence is very short
    import programme_settings
    word_size = int(programme_settings.settings["BLAST"]["word_size"])
    if len(seq) < word_size:
        print(f"[WARNING] Query sequence is shorter than word size ({word_size}). "
              f"BLAST search may return no results.")

    return seq


def validate_database(db):
    """
    Checks the database file exists and is readable.
    """
    if not os.path.isfile(db):
        print(f"[ERROR] Database file not found: '{db}'")
        print("        Check the filename and that you are running from the correct directory.")
        sys.exit(1)
    if not db.endswith(".fasta"):
        print(f"[WARNING] Database file '{db}' does not have a .fasta extension.")
    return db


def get_output_stream(output_path):
    """
    Returns an open file handle if output path is given, otherwise None.
    Results are always printed to screen; file is written in addition if specified.
    """
    if output_path is None:
        return None
    try:
        fh = open(output_path, 'w')
        print(f"[INFO] Results will also be saved to: {output_path}")
        return fh
    except IOError as e:
        print(f"[ERROR] Cannot open output file '{output_path}': {e}")
        sys.exit(1)


def tee_print(*args_p, fh=None, **kwargs):
    """Prints to screen and optionally to a file simultaneously."""
    print(*args_p, **kwargs)
    if fh:
        print(*args_p, file=fh, **kwargs)


# ---------------------------------------------------------------------------
# Settings summary
# ---------------------------------------------------------------------------
def print_settings_summary(args, settings):
    """Prints a clear summary of all parameters before the search runs."""
    print("\n" + "─" * 66)
    print("  Run Parameters")
    print("─" * 66)
    print(f"  Mode              : {args.mode.upper()}")
    print(f"  Query sequence    : {args.query[:55]}{'...' if len(args.query)>55 else ''}")
    print(f"  Query length      : {len(args.query)} residues")
    print(f"  Database          : {args.database}")
    print(f"  Max scores        : {args.max_scores}")
    print(f"  BLOSUM matrix     : {settings['DEFAULT']['blosum']}")
    print(f"  Gap penalty       : {settings['DEFAULT']['seq_gap']}")
    print(f"  Word size         : {settings['BLAST']['word_size']}")
    print(f"  T-score threshold : {settings['BLAST']['tscore']}")
    print(f"  Max extension len : {settings['BLAST']['max_extension_length']}")
    print(f"  Output file       : {args.output if args.output else 'stdout (screen only)'}")
    print("─" * 66 + "\n")


# ---------------------------------------------------------------------------
# Mode runners
# ---------------------------------------------------------------------------
def run_blast(args, settings, output_fh=None):
    """Runs the BLAST 101 search."""
    settings["DEFAULT"]["query_sequence"] = args.query
    settings["DEFAULT"]["database"]       = args.database
    settings["BLAST"]["max_scores"]       = str(args.max_scores)

    import programme_settings
    programme_settings.settings = settings

    import blast_101_search as b
    b.qsequence      = args.query
    b.query_sequence = __import__("create_seq_word_dict").create_word_dict(args.query)
    b.max_scores     = args.max_scores

    tee_print("[INFO] Starting BLAST search...", fh=output_fh)

    # Run the search phase
    res = b.process_fasta_file()

    # Print results to screen and optionally to file
    if output_fh:
        # Redirect stdout temporarily to capture output to file as well
        old_stdout = sys.stdout
        sys.stdout = output_fh
        b.print_final_results(res)
        sys.stdout = old_stdout
        # Also print to screen
        b.print_final_results(res)
    else:
        b.print_final_results(res)


def run_sw(args, settings, output_fh=None):
    """Runs the Smith-Waterman exhaustive search."""
    settings["DEFAULT"]["query_sequence"] = args.query
    settings["DEFAULT"]["database"]       = args.database
    settings["SWSEARCH"]["max_sw_scores"] = str(args.max_scores)

    import programme_settings
    programme_settings.settings = settings

    tee_print("[INFO] Starting Smith-Waterman search (this may be slow for large databases)...",
              fh=output_fh)

    import smith_waterman_p as SW
    import process_fasta_file as pff
    from operator import itemgetter
    import calc_bit_and_evalues as cbe
    import time

    query     = args.query
    aligntime = 0

    def processSW(myline_database):
        nonlocal aligntime
        t3  = time.time()
        res = SW.perform_smith_waterman(query, myline_database.upper(), False, False)
        t4  = time.time()
        aligntime += (t4 - t3)
        return res

    t0 = time.time()
    pff.res        = []
    pff.bestscore  = 0
    res = pff.process_fasta_file(args.database, processSW, args.max_scores)
    t1  = time.time()

    res.sort(key=itemgetter(2), reverse=True)

    k   = float(settings["DEFAULT"]["k"])
    lam = float(settings["DEFAULT"]["lam"])

    tee_print("\n###############################################################", fh=output_fh)
    tee_print("SW Alignment Results", fh=output_fh)
    tee_print(f"Query: {query}", fh=output_fh)
    tee_print(f"Database: {args.database}", fh=output_fh)
    tee_print("###############################################################\n", fh=output_fh)

    for i in res:
        line = (f'"{str(i[0])[1:40]}..."\t'
                f'Score={i[2]}\t'
                f'BitScore={cbe.get_bit_score_s(i[2], k, lam)}\t'
                f'E={cbe.get_expect_s(i[2], k, lam)}')
        tee_print(line, fh=output_fh)

    elapsed = t1 - t0
    tee_print(f"\nTime taken: {elapsed:.2f} secs", fh=output_fh)
    if elapsed > 0:
        tee_print(f"SW time: {aligntime / elapsed * 100:.1f}%", fh=output_fh)
    tee_print("~~~~~~~~~Finished~~~~~~~~~", fh=output_fh)


def run_stats(args, settings, output_fh=None):
    """Runs the Gumbel distribution fit on random SW score data."""
    tee_print("[INFO] Running statistical analysis (Gumbel distribution fit)...", fh=output_fh)
    import build_expect_scores
    build_expect_scores.run_trial()


def run_tests():
    """Runs the full automated test suite using the venv Python."""
    print("[INFO] Running automated test suite...\n")
    import subprocess
    test_file   = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_blast101.py")
    venv_python = "/home/s2793337/Bioinformatics_Algorithms/ICA/Blast101_code/BA_ICA/bin/python3"
    result      = subprocess.run([venv_python, test_file],
                                 cwd=os.path.dirname(os.path.abspath(__file__)))
    sys.exit(result.returncode)


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------
def build_parser(settings):
    """Builds and returns the argparse argument parser."""

    default_query = settings["DEFAULT"]["query_sequence"]
    default_db    = settings["DEFAULT"]["database"]
    default_max   = int(settings["BLAST"]["max_scores"])

    parser = argparse.ArgumentParser(
        prog="blast101_cli.py",
        description=(
            "Blast101 Command Line Interface\n"
            "Runs BLAST-like search, Smith-Waterman search, "
            "statistical analysis, or the automated test suite."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python blast101_cli.py --mode blast\n"
            "  python blast101_cli.py --mode blast --query PWNAAPLHNFGEDFLQ "
            "--database uniprot_bit2.fasta\n"
            "  python blast101_cli.py --mode blast --database uniprot_bit2.fasta "
            "--output results.txt\n"
            "  python blast101_cli.py --mode sw --query PWNAAPLHNFGEDFLQ\n"
            "  python blast101_cli.py --mode stats\n"
            "  python blast101_cli.py --mode test\n"
        )
    )

    parser.add_argument(
        "--mode", "-m",
        choices=MODES,
        required=True,
        help="Which mode to run: blast | sw | stats | test"
    )
    parser.add_argument(
        "--query", "-q",
        default=default_query,
        help="Protein query sequence. Default: sequence from settings.ini"
    )
    parser.add_argument(
        "--database", "-d",
        default=default_db,
        help=f"FASTA database file to search. Default: {default_db}"
    )
    parser.add_argument(
        "--max-scores", "-n",
        type=int,
        default=default_max,
        dest="max_scores",
        help=f"Maximum number of results to return. Default: {default_max}"
    )
    parser.add_argument(
        "--output", "-o",
        default=None,
        dest="output",
        help="Optional output file to save results (e.g. results.txt). "
             "Results are always printed to screen; this saves a copy to file."
    )

    return parser


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print(BANNER)

    # Read settings first — needed for defaults in the parser
    import programme_settings
    programme_settings.read()
    settings = programme_settings.settings

    parser = build_parser(settings)
    args   = parser.parse_args()

    # Test mode skips sequence/database validation
    if args.mode == "test":
        run_tests()
        return

    # --- Validate inputs ---
    args.query    = validate_sequence(args.query)
    args.database = validate_database(args.database)

    if args.max_scores < 1:
        print("[ERROR] --max-scores must be at least 1.")
        sys.exit(1)

    # --- Open output file if requested ---
    output_fh = get_output_stream(args.output)

    # --- Print settings summary ---
    print_settings_summary(args, settings)

    # --- Dispatch to correct mode ---
    if args.mode == "blast":
        run_blast(args, settings, output_fh)
    elif args.mode == "sw":
        run_sw(args, settings, output_fh)
    elif args.mode == "stats":
        run_stats(args, settings, output_fh)

    # --- Close output file ---
    if output_fh:
        output_fh.close()
        print(f"\n[INFO] Results saved to: {args.output}")


if __name__ == "__main__":
    main()
