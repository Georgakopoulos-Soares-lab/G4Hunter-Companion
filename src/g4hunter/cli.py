"""
Command-line interface for G4Hunter.

This module provides the main entry point for the g4hunter command-line tool.
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional
from .pg4_extraction import G4Hunter


def setup_logging(verbose: bool = False) -> None:
    """Setup logging configuration."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )


def create_parser() -> argparse.ArgumentParser:
    """Create and configure the argument parser."""
    parser = argparse.ArgumentParser(
        description="G4Hunter - Extract G-quadruplex forming sequences from DNA using sliding window scoring.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage
  g4hunter input.fasta

  # Custom parameters
  g4hunter input.fasta --window-size 30 --min-score 2.0 --outdir results

  # With consensus motif analysis
  g4hunter input.fasta --parse-consensus --merge-sequences

  # Direct repeat analysis parameters
  g4hunter input.fasta --min-arm-length 15 --max-spacer-length 10

For more information, visit: https://github.com/yourusername/g4hunter
        """
    )
    
    # Required arguments
    parser.add_argument(
        "infile",
        type=str,
        help="Input FASTA file containing DNA sequences (.fasta, .fa, .fna, optionally .gz compressed)"
    )
    
    # G4Hunter parameters
    g4_group = parser.add_argument_group("G4Hunter Parameters")
    g4_group.add_argument(
        "--window-size", "--window_size",
        type=int,
        default=25,
        help="Sliding window size for G4Hunter scoring (default: 25)"
    )
    g4_group.add_argument(
        "--min-score", "--min_score",
        type=float,
        default=1.5,
        help="Minimum G4Hunter score threshold for sequence extraction (default: 1.5)"
    )
    
    # Consensus motif parameters
    consensus_group = parser.add_argument_group("Consensus Motif Parameters")
    consensus_group.add_argument(
        "--min-grun-length", "--min_grun_length",
        type=int,
        default=3,
        help="Minimum length of G-runs for consensus motif detection (default: 3)"
    )
    consensus_group.add_argument(
        "--min-loop-length", "--min_loop_length",
        type=int,
        default=1,
        help="Minimum loop length between G-runs for consensus motif (default: 1)"
    )
    consensus_group.add_argument(
        "--max-loop-length", "--max_loop_length",
        type=int,
        default=7,
        help="Maximum loop length between G-runs for consensus motif (default: 7)"
    )
    consensus_group.add_argument(
        "--regex-multiplicity", "--regex_multiplicity",
        type=int,
        default=3,
        help="Number of G-run repeats required for consensus motif (default: 3)"
    )
    
    # Direct repeat parameters
    repeat_group = parser.add_argument_group("Direct Repeat Parameters")
    repeat_group.add_argument(
        "--min-arm-length", "--min_arm_length",
        type=int,
        default=10,
        help="Minimum arm length for direct repeat detection (default: 10)"
    )
    repeat_group.add_argument(
        "--min-spacer-length", "--min_spacer_length",
        type=int,
        default=0,
        help="Minimum spacer length for direct repeats (default: 0)"
    )
    repeat_group.add_argument(
        "--max-spacer-length", "--max_spacer_length",
        type=int,
        default=8,
        help="Maximum spacer length for direct repeats (default: 8)"
    )
    
    # Analysis options
    analysis_group = parser.add_argument_group("Analysis Options")
    analysis_group.add_argument(
        "--parse-consensus", "--parse_consensus",
        action="store_true",
        help="Parse consensus pG4 motifs using regex patterns"
    )
    analysis_group.add_argument(
        "--merge-sequences", "--merge_sequences",
        action="store_true",
        help="Merge overlapping sequences in output"
    )
    
    # Output options
    output_group = parser.add_argument_group("Output Options")
    output_group.add_argument(
        "--outdir",
        type=str,
        default="G4_results",
        help="Output directory for results files (default: G4_results)"
    )
    
    # General options
    general_group = parser.add_argument_group("General Options")
    general_group.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Enable verbose logging"
    )
    general_group.add_argument(
        "--version",
        action="version",
        version="%(prog)s 2.0.0"
    )
    
    return parser


def validate_args(args: argparse.Namespace) -> None:
    """Validate command-line arguments."""
    # Check input file exists
    infile = Path(args.infile)
    if not infile.exists():
        raise FileNotFoundError(f"Input file not found: {args.infile}")
    
    # Validate numeric parameters
    if args.window_size < 1:
        raise ValueError("Window size must be >= 1")
    
    if args.min_score < 0:
        raise ValueError("Minimum score must be >= 0")
    
    if args.min_grun_length < 1:
        raise ValueError("Minimum G-run length must be >= 1")
    
    if args.regex_multiplicity < 1:
        raise ValueError("Regex multiplicity must be >= 1")
    
    if args.min_loop_length < 0:
        raise ValueError("Minimum loop length must be >= 0")
    
    if args.max_loop_length < args.min_loop_length:
        raise ValueError("Maximum loop length must be >= minimum loop length")
    
    if args.min_arm_length < 1:
        raise ValueError("Minimum arm length must be >= 1")
    
    if args.min_spacer_length < 0:
        raise ValueError("Minimum spacer length must be >= 0")
    
    if args.max_spacer_length < args.min_spacer_length:
        raise ValueError("Maximum spacer length must be >= minimum spacer length")


def main(argv: Optional[list[str]] = None) -> int:
    """
    Main entry point for the g4hunter command-line tool.
    
    Args:
        argv: Command-line arguments (for testing)
        
    Returns:
        int: Exit code (0 for success, 1 for error)
    """
    parser = create_parser()
    args = parser.parse_args(argv)
    
    # Setup logging
    setup_logging(args.verbose)
    logger = logging.getLogger(__name__)
    
    try:
        # Validate arguments
        validate_args(args)
        
        # Create output directory
        outdir = Path(args.outdir).resolve()
        outdir.mkdir(exist_ok=True, parents=True)
        
        # Initialize G4Hunter
        hunter = G4Hunter(
            outdir=outdir,
            window_size=args.window_size,
            min_score=args.min_score,
            regex_multiplicity=args.regex_multiplicity,
            min_grun_length=args.min_grun_length,
            min_loop_length=args.min_loop_length,
            max_loop_length=args.max_loop_length,
            min_arm_length=args.min_arm_length,
            min_spacer_length=args.min_spacer_length,
            max_spacer_length=args.max_spacer_length,
        )
        
        # Run analysis
        logger.info(f"Starting G4Hunter analysis on {args.infile}")
        scores = hunter.run(
            infile=args.infile,
            parse_consensus=args.parse_consensus,
            merge_sequences=args.merge_sequences
        )
        
        # Log summary
        if scores:
            logger.info(f"Analysis complete. Found {len(scores)} G4 regions.")
            logger.info(f"Score range: {min(scores):.2f} - {max(scores):.2f}")
            logger.info(f"Mean score: {sum(scores)/len(scores):.2f}")
        else:
            logger.info("Analysis complete. No G4 regions found.")
        
        logger.info(f"Results written to {outdir}")
        return 0
        
    except KeyboardInterrupt:
        logger.error("Analysis interrupted by user")
        return 1
    except Exception as e:
        logger.error(f"Error during analysis: {e}")
        if args.verbose:
            logger.exception("Full traceback:")
        return 1


if __name__ == "__main__":
    sys.exit(main())