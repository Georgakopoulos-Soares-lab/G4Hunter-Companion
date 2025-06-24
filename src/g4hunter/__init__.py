"""
G4Hunter - Extract G-quadruplex forming sequences from DNA using sliding window scoring.

This package provides tools for identifying potential G-quadruplex (G4) forming sequences
in DNA using the G4Hunter algorithm and consensus motif analysis.

Author: Nikol Chantzi
Version: 1.0.0
"""

__version__ = "1.0.0"
__author__ = "Nikol Chantzi"
__email__ = "nicolechantzi@gmail.com"

from .pg4_extraction import G4Hunter
from .cli import main

__all__ = ["G4Hunter", "main"]
