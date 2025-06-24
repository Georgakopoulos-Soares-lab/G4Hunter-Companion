import attr 
from attr import field
import attrs
from termcolor import colored
from Bio.SeqIO.FastaIO import SimpleFastaParser
from pathlib import Path
import gzip 
import numpy as np
import pandas as pd
import csv
from typing import Iterable, ClassVar
import re
from collections import defaultdict
from pybedtools import BedTool
import logging

logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(levelname)s - %(message)s')

@attr.s(slots=True)
class G4Hunter:

    outdir: str = field(default=None, init=True)
    outfile: str = field(default=None, init=False)
    # pG4s G4Hunter Parameters
    window_size: int = field(converter=int, validator=attrs.validators.ge(0), default=25, init=True, repr=True)
    min_score: float = field(converter=float, validator=attrs.validators.ge(0.0), default=1.5, init=True, repr=True)
    # pG4s Consensus Motif Parameters
    regex_multiplicity: int = field(converter=int, validator=attrs.validators.ge(0), default=3, init=True, repr=True)
    min_grun_length: int = field(converter=int, validator=attrs.validators.ge(0), default=3, init=True, repr=True)
    min_loop_length: int = field(converter=int, validator=attrs.validators.ge(0), default=1, init=True, repr=True)
    max_loop_length: int = field(converter=int, validator=attrs.validators.ge(0), default=7, init=True, repr=True)
    regex_motifs: dict[str, str] = field(init=False, repr=True)
    # pG4s DIRECT Repeats Parameters
    min_arm_length: int = field(converter=int, validator=attrs.validators.ge(0), default=10)
    min_spacer_length: int = field(converter=int, validator=attrs.validators.ge(0), default=0)
    max_spacer_length: int = field(converter=int, validator=attrs.validators.ge(0), default=8)
    direct_consensus_motif: str = field(init=False, repr=True)
    # FIELDNAMES for FOUT
    FIELDNAMES: ClassVar[list[str]] = ["seqID", 
                                       "start", 
                                       "end", 
                                       "sequence", 
                                       "length", 
                                       "score", 
                                       "strand", 
                                       # NBR,
                                       "type", 
                                       "method"]

    def __attrs_post_init__(self) -> None:
        self.outdir = Path(self.outdir).resolve()
        self.outdir.mkdir(exist_ok=True)
        ########################################################################################################
        # Validate and set parameters to ensure they meet minimum requirements
        self.min_arm_length = max(self.min_arm_length, 1)  # Ensure minimum arm length is at least 1
        self.min_spacer_length = max(self.min_spacer_length, 0)  # Ensure minimum spacer length is non-negative
        self.max_spacer_length = max(self.max_spacer_length, self.min_spacer_length)  # Ensure max spacer length is at least min spacer length
        self.window_size = max(self.window_size, 1)  # Ensure window size is at least 1
        self.min_grun_length = max(self.min_grun_length, 1)  # Ensure minimum G-run length is at least 1
        self.min_loop_length = max(self.min_loop_length, 0)  # Ensure minimum loop length is non-negative
        self.max_loop_length = max(self.max_loop_length, self.min_loop_length)  # Ensure max loop length is at least min loop length        
        ########################################################################################################
        # Initialize regex patterns for direct repeats and G4 motifs
        self.direct_consensus_motif: str = r"(\w{%s,})(\w{%s,%s})\1" % (self.min_arm_length,
                                                                 self.min_spacer_length,
                                                                 self.max_spacer_length)
        self.regex_motifs: dict[str, str] = {
                            "Crick": r"c{%s,}\w{%s,%s}" % (self.min_grun_length, 
                                                           self.min_loop_length, 
                                                           self.max_loop_length) * self.regex_multiplicity + r"c{%s,}" % self.min_grun_length,
                            "Watson": r"g{%s,}\w{%s,%s}" % (self.min_grun_length, 
                                                            self.min_loop_length, 
                                                            self.max_loop_length) * self.regex_multiplicity + r"g{%s,}" % self.min_grun_length
                            }
        ########################################################################################################
        print(colored(f"Loading Direct Repeat Consensus Motif: {self.direct_consensus_motif}", "green"))
        print(colored(f"Chosen Direct Repeat Minimum Arm Length: {self.min_arm_length}.", "blue"))
        print(colored(f"Chosen Direct Repeat Minimum Spacer Length: {self.min_spacer_length}.", "blue"))
        print(colored(f"Chosen Direct Repeat Maximum Spacer Length: {self.max_spacer_length}.", "blue"))
        ########################################################################################################
        print(colored(f"Loading G4 engine...", "green"))
        print(colored(f"Chosen Window Size: {self.window_size}.", "blue"))
        print(colored(f"Chosen Minimum Score Threshold: {self.min_score}.", "blue"))
        print(colored(f"Chosen Minimum G-Run Length: {self.min_grun_length}.", "blue"))
        print(colored(f"Chosen Minimum Loop Length: {self.min_loop_length}.", "blue"))
        print(colored(f"Chosen Maximum Loop Length: {self.max_loop_length}.", "blue"))
        for motif, color in zip(self.regex_motifs.keys(), ["green", "magenta"]):
            print(colored(f"{motif}: {self.regex_motifs[motif]}", color))

    @staticmethod 
    def get_strand_symbol(strand: str) -> str:
        """
        Convert strand name to strand symbol.
        
        Args:
            strand: Strand name, either "Crick" or "Watson"
            
        Returns:
            str: "+" for Watson strand, "-" for Crick strand
        """
        if strand == "Crick":
            return "-"
        if strand == "Watson":
            return "+"
        raise ValueError(f"Invalid strand name: `{strand}`.")

    def extract_consensus_pG4s(self, accession: str, save_output: bool = True) -> pd.DataFrame:
        results_df = []
        for seqID, seq in G4Hunter.parse_fasta(accession):
            results = self._extract_consensus_pG4s(seq)
            results.loc[:, "seqID"] = seqID
            results_df.append(results)
        results_df = pd.concat(results_df, ignore_index=True)
        if save_output:
            accession_id = G4Hunter.extract_accession_id(accession)
            consensus_outfile = self.outdir / f"{accession_id}_pG4s.consensus.tsv"
            print(colored(f"Outsourcing pG4s consensus motifs results to: ➡️  {self.outfile}", "green"))
            results_df.to_csv(consensus_outfile, 
                              sep="\t", 
                              mode="w", 
                              header=True, 
                              index=False)
            print(colored(f"Finished processing {accession}. Results written to {self.outfile}", "green"))
        return results_df
    
    def _extract_consensus_pG4s(self, seq: str) -> dict[str, list]:
        """
        Extract consensus putative G4 sequences from DNA sequence using regex patterns.
        
        Args:
            seq: DNA sequence string to search for G4 motifs
            
        Returns:
            dict[str, list]: Dictionary containing lists of start, end, sequence, length, score, and strand information
        """
        results = defaultdict(list)
        for strand, consensus_motif in self.regex_motifs.items():
            strand_symbol = G4Hunter.get_strand_symbol(strand)
            matches = re.finditer(consensus_motif, seq)
            for match in matches:
                start = int(match.start())
                end = int(match.end())
                sequence = match.group()
                sequence_length = len(sequence)
                results["start"].append(start)
                results["end"].append(end)
                results["sequence"].append(sequence)
                results["length"].append(sequence_length)
                results["score"].append(".")
                results["strand"].append(strand_symbol)
        
        # Sort results by start position for consistent output
        if results["start"]:
            sorted_indices = sorted(range(len(results["start"])), key=lambda i: results["start"][i])
            for key in results:
                results[key] = [results[key][i] for i in sorted_indices]
        return results
    
    def extract_direct_repeats(self, accession: str, save_output: bool = True) -> pd.DataFrame:
        results_df = []
        for seqID, seq in G4Hunter.parse_fasta(accession):
            results = self._extract_direct_repeats(seq, return_frame=True)
            results.loc[:, "seqID"] = seqID
            results_df.append(results)
        results_df = pd.concat(results_df, ignore_index=True)
        if save_output:
            outfile = self.outdir / f"{G4Hunter.extract_accession_id(accession)}_pG4s.direct_repeats.tsv"
            print(colored(f"Outsourcing direct repeat results to: ➡️  {self.outfile}", "green"))
            results_df.to_csv(outfile, 
                              sep="\t", 
                              mode="w", 
                              header=True, 
                              index=False)
            print(colored(f"Finished processing {accession}. Results written to {self.outfile}", "green"))
        return results_df

    def _extract_direct_repeats(self, seq: str, return_frame: bool = False) -> dict[str, list]:
        """
        Extract direct repeat sequences from DNA sequence using regex patterns.
        
        Args:
            seq: DNA sequence string to search for direct repeats
            
        Returns:
            dict[str, list]: Dictionary containing lists of start, end, sequence, arm sequence, spacer sequence, length, score, and strand information
        """
        results = defaultdict(list)
        matches = re.finditer(self.direct_consensus_motif, seq)
        for match in matches:
           start = int(match.start())
           end = int(match.end())
           sequence = match.group()
           arm_sequence = match.group(1)
           arm_length = len(arm_sequence)
           spacer_sequence = match.group(2)
           spacer_length = len(spacer_sequence)
           sequence_length = len(sequence)
           results["start"].append(start)
           results["end"].append(end)
           results["sequence"].append(sequence)
           results["sequence_of_arm"].append(arm_sequence)
           results["arm_length"].append(arm_length)
           results["sequence_of_spacer"].append(spacer_sequence)
           results["spacer_length"].append(spacer_length)
           results["length"].append(sequence_length)
           results["score"].append(".")
           results["strand"].append("+")
           results["type"].append("DR")
           results["method"].append("Direct Repeat")
        
        # Sort results by start position for consistent output
        # Because there is only one motif, these are already sorted
        # if results["start"]:
        #    sorted_indices = sorted(range(len(results["start"])), key=lambda i: results["start"][i])
        #    for key in results:
        #        results[key] = [results[key][i] for i in sorted_indices]
        return results

    @staticmethod
    def parse_fasta(fasta: str) -> Iterable[tuple[str, str]]:
        """
        Parse a FASTA file and yield tuples of (seqID, sequence).

        :param fasta: Path to the FASTA file.
        :return: An iterable of tuples containing seqID and sequence.
        """
        if not Path(fasta).is_file():
            raise FileNotFoundError(f"File not found: {fasta}")
        if Path(fasta).name.endswith(".gz"):
            f = gzip.open(fasta, "rt")
        else:
            f = open(fasta, mode="r", encoding="UTF-8")
        for record in SimpleFastaParser(f):
             seqID = record[0].split(" ")[0]
             seq = record[1].lower().strip()
             yield seqID, seq
        f.close()

    @staticmethod 
    def extract_accession_id(accession: str) -> str:
        """
        Extract accession ID from filename by removing file extension.
        
        Args:
            accession: Filename or path containing sequence data
            
        Returns:
            str: Accession ID without file extension
        """
        accession = Path(accession.strip()).name
        if ".fna" in accession:
            return accession.split(".fna")[0]
        if ".fa" in accession:
            return accession.split(".fa")[0]
        if ".fasta" in accession:
            return accession.split(".fasta")[0]
        raise ValueError(f"Invalid accession format: `{accession}`. Expected .fna, .fa, or .fasta suffix.")
    
    def merge_overlapping_regions(self, save_file: bool = True) -> pd.DataFrame:
        """
        Merge overlapping regions in the output file using BedTool.

        Returns:
            None: Merged regions are written to the output file.
        """
        if not self.outfile:
            raise ValueError("Output file not set. Please run the `run` method first.")
        logging.info(f"Merging overlapping regions in {self.outfile}")
        # Read the output file as a BedTool object        
        df_g4_merged = pd.read_table(
            BedTool.from_dataframe(
                pd.read_table(self.outfile,
                          usecols=["seqID", "start", "end", "sequence"]
                          )
                )
                .sort()
                .merge(c="2,4", 
                       delim=";", 
                       o="collapse,collapse")
                .fn,
            header=None,
            dtype={"start": np.int32,
                   "end": np.int32,},
            names=["seqID", "start", "end", "multiple_starts", "multiple_sequences"],
        )
        merged_g4_sequences = defaultdict(list)
        # Iterate through the merged DataFrame and construct the merged sequences
        for _, row in df_g4_merged.iterrows():
            seqID = row["seqID"]
            start = row["start"]
            end = row["end"]
            multiple_starts = list(map(int, row["multiple_starts"].split(";")))
            multiple_sequences = row["multiple_sequences"].split(";")
            merged_sequence = ""
            pstart = start
            for sstart, seq in zip(multiple_starts[1:], multiple_sequences[:-1]):
                merged_sequence += seq[:sstart - pstart]
                pstart = sstart
            merged_sequence += multiple_sequences[-1]
            # validate length 
            assert len(merged_sequence) == end - start, f"Length mismatch for {seqID}: expected {end - start}, got {len(merged_sequence)}"
            merged_g4_sequences["seqID"].append(seqID)
            merged_g4_sequences["start"].append(start)
            merged_g4_sequences["end"].append(end)
            merged_g4_sequences["sequence"].append(merged_sequence)
            merged_g4_sequences["length"].append(len(merged_sequence))
        merged_g4_sequences = pd.DataFrame(merged_g4_sequences)
        merged_g4_sequences.loc[:, "type"] = "pG4"
        merged_g4_sequences.loc[:, "method"] = "G4Hunter"
        if save_file:
            merged_outfile = self.outdir / f"{self.outfile.name.split('_pG4s')[0]}_pG4s.merged.g4_hunter.tsv"
            merged_g4_sequences.to_csv(merged_outfile, 
                                       sep="\t", 
                                       mode="w",
                                       header=True,
                                       index=False)
            logging.info(f"Merged regions saved to {self.outfile}")
        return merged_g4_sequences
    
    def run(self, infile: str, 
            parse_consensus: bool = False,
            merge_sequences: bool = False) -> list[float]:
        """
        Extract G4 (G-quadruplex) forming sequences from FASTA file using sliding window scoring.
        
        Args:
            infile: Path to input FASTA file
            parse_consensus: Whether to also parse consensus G4 motifs (currently unused)

        Outputs:
             Writes to disk the results of G4Hunter analysis in TSV format.
             If `parse_consensus` is True, also writes consensus pG4 motifs to a separate TSV file.
             If `merge_sequences` is True, merges overlapping G4 regions and writes to a new TSV file.
            
        Returns:
            list[float]: List of G4 scores for extracted sequences
        """
        accession_id = G4Hunter.extract_accession_id(infile)
        self.outfile = self.outdir / f"{accession_id}_pG4s.g4_hunter.tsv"
        consensus_outfile = self.outdir / f"{accession_id}_pG4s.consensus.tsv"
        # direct_outfile = self.outdir / f"{accession_id}_pG4s.direct_repeats.tsv"
        if parse_consensus:
            consensus_fout = open(consensus_outfile, "w", encoding="UTF-8")
            consensus_writer = csv.DictWriter(consensus_fout, 
                                              delimiter="\t", 
                                              fieldnames=G4Hunter.FIELDNAMES)
            consensus_writer.writeheader()
            print(colored(f"Outsourcing consensus results to: ➡️  {consensus_outfile}", "green"))
        else:
            consensus_fout = None
            consensus_writer = None
        print(colored(f"Outsourcing results to: ➡️  {self.outfile}", "green"))
        with self.outfile.open("w", encoding="UTF-8") as fout:
            writer = csv.DictWriter(fout, delimiter="\t", fieldnames=G4Hunter.FIELDNAMES)
            writer.writeheader()
            for seqID, seq in self.parse_fasta(infile):
                if parse_consensus and consensus_writer:
                    # Extract consensus pG4s
                    results = self._extract_consensus_pG4s(seq)
                    for i in range(len(results["start"])):
                        consensus_writer.writerow({
                            "seqID": seqID,
                            "start": results["start"][i],
                            "end": results["end"][i],
                            "sequence": results["sequence"][i],
                            "length": results["length"][i],
                            "score": results["score"][i],
                            "strand": results["strand"][i],
                            "type": "pG4",
                            "method": "Consensus Motif"
                        })
                score_list = self.get_base_score(seq)
                avg_scores = [
                    sum(score_list[i:i+self.window_size]) / self.window_size for i in range(len(score_list) - (self.window_size - 1))
                ]
                # get_g4 function in original code|START
                consecutive_high_scoring = [
                    i for i in range(len(avg_scores)) if abs(avg_scores[i]) >= self.min_score
                ]
                if len(consecutive_high_scoring) > 0:
                    MSCORE = self.write_output(writer,
                                               seq,
                                               seqID,
                                               avg_scores, 
                                               consecutive_high_scoring,
                                               )
                # < get_g4 function in original code|END
        # return LScoreSeq, LSeq , LNumber, LseqID
        if merge_sequences:
            self.merge_overlapping_regions(save_file=True)
        if consensus_fout:
            consensus_fout.close()
            print(colored(f"Finished processing {infile}. Consensus results written to {consensus_outfile}", "green"))
        print(colored(f"Finished processing {infile}. Results written to {self.outfile}", "green"))
        return MSCORE

    @staticmethod 
    def determine_strand(sequence: str) -> str:
        """
        Determine the strand of the sequence based on G vs C base prevalence.
        
        Args:
            sequence: DNA sequence string
            
        Returns:
            str: "+" if G is more prevalent (Watson strand), "-" if C is more prevalent (Crick strand)
        """
        g_count = sequence.count("g")
        c_count = sequence.count("c")
        return "+" if g_count >= c_count else "-"

    def write_output(self, writer: csv.DictWriter, 
                     sequence: str, 
                     seqID: str, 
                     score_list: list[int], 
                     consecutive_high_scoring: list[int]) -> list[float]:
        """
        Write G4 sequences to output file by merging consecutive high-scoring windows.
        
        Args:
            writer: CSV writer object for output file
            sequence: Original DNA sequence
            seqID: Sequence identifier
            score_list: List of sliding window average scores
            consecutive_high_scoring: List of indices where scores exceed threshold
            
        Returns:
            list[float]: List of scores for merged G4 regions
        """
        i, k, I = 0, 0, 0
        a = b = consecutive_high_scoring[0]
        MSCORE: list[float] = []
        if (len(consecutive_high_scoring) > 1):
            c = consecutive_high_scoring[i+1]
            while (i < len(consecutive_high_scoring)-2):
                if (c == b+1):
                    k += 1
                    i += 1
                else:
                    I += 1
                    seq = sequence[a:a + self.window_size + k]
                    new_score_list = self.get_base_score(seq)
                    avg_score_list = round(np.mean(new_score_list), 2)
                    writer.writerow({
                         "seqID": seqID,
                         "start": a,
                         "end": a + k + self.window_size,
                         "sequence": seq,
                         "length": len(seq),
                         "score": avg_score_list,
                         "strand": G4Hunter.determine_strand(seq),
                         # "NBR": I,
                         "type": "pG4",
                         "method": "G4Hunter"
                        }
                    )
                    MSCORE.append(abs(avg_score_list))
                    k = 0
                    i += 1
                    a = consecutive_high_scoring[i]
                b = consecutive_high_scoring[i]
                c = consecutive_high_scoring[i + 1]
            I += 1
            seq = sequence[a:a + self.window_size + k + 1]
            new_score_list = self.get_base_score(seq)
            avg_score_list = round(np.mean(new_score_list), 2)
            writer.writerow({
                    "seqID": seqID,
                    "start": a,
                    "end": a + k + self.window_size + 1,
                    "sequence": seq,
                    "length": len(seq),
                    "score": avg_score_list,
                    "strand": G4Hunter.determine_strand(seq),
                    # "NBR": I,
                    "type": "pG4",
                    "method": "G4Hunter"
               }
            )
            MSCORE.append(abs(avg_score_list))
        else:
            I += 1
            seq = sequence[a:a + self.window_size]
            writer.writerow({
                    "seqID": seqID,
                    "start": a,
                    "end": a + k + self.window_size,
                    "sequence": seq,
                    "length": len(seq),
                    "score": score_list[a],
                    "strand": G4Hunter.determine_strand(seq),
                    # "NBR": I,
                    "type": "pG4",
                    "method": "G4Hunter"
               }
            )
            MSCORE.append(abs(score_list[a]))
        return MSCORE
    
    def get_base_score(self, line: str) -> tuple[str, list[int]]:
        """
        Calculate G4Hunter score for each base in the DNA sequence.
        
        G bases get positive scores (1-4) based on consecutive G runs.
        C bases get negative scores (-1 to -4) based on consecutive C runs.
        Other bases get score 0.
        
        Args:
            line: DNA sequence string
            
        Returns:
            list[int]: List of scores for each base position
        """
        item, score_list = 0, []
        # calcule le item de chaque base et la stock dans score_list
        while (item < len(line)):
            if (item < len(line) and (line[item]=="G" or line[item]=="g")):
                score_list.append(1)
                if(item+1< len(line) and (line[item+1]=="G" or line[item+1]=="g")):
                    score_list[item]=2
                    score_list.append(2)
                    if (item+2< len(line) and (line[item+2]=="G" or line[item+2]=="g")):
                        score_list[item+1]=3
                        score_list[item]=3
                        score_list.append(3)
                        if (item+3< len(line) and (line[item+3]=="G" or line[item+3]=="g")):
                            score_list[item]=4
                            score_list[item+1]=4
                            score_list[item+2]=4
                            score_list.append(4)
                            item=item+1
                        item=item+1
                    item=item+1
                item=item+1
                while(item < len(line) and (line[item]=="G" or line[item]=="g")):
                        score_list.append(4)
                        item=item+1

            elif (item < len(line) and line[item]!="G" and line[item]!="g" and line[item]!= "C" and line[item]!="c" ):
                        score_list.append(0)
                        item=item+1
                
            elif(item < len(line) and (line[item]=="C" or line[item]=="c")):
                score_list.append(-1)
                if(item+1< len(line) and (line[item+1]=="C" or line[item+1]=="c" )):
                    score_list[item]=-2
                    score_list.append(-2)
                    if (item+2< len(line) and (line[item+2]=="C" or line[item+2]=="c" )):
                        score_list[item+1]=-3
                        score_list[item]=-3
                        score_list.append(-3)
                        if (item+3< len(line) and (line[item+3]=="C" or line[item+3]=="c"  )):
                            score_list[item]=-4
                            score_list[item+1]=-4
                            score_list[item+2]=-4
                            score_list.append(-4)
                            item=item+1
                        item=item+1   
                    item=item+1
                item=item+1
                while(item < len(line) and (line[item]=="C" or line[item]=="c")):
                    score_list.append(-4)
                    item=item+1
            else:
                    item=item+1 
        # return line, score_list
        return score_list

if __name__ == "__main__":
    import argparse 
    parser = argparse.ArgumentParser(description="G4Hunter - Extract G-quadruplex forming sequences from DNA using sliding window scoring.")
    parser.add_argument("infile", type=str, help="Input FASTA file containing DNA sequences (.fasta, .fa, .fna, optionally .gz compressed)")
    parser.add_argument("--min_grun_length", type=int, default=3, help="Minimum length of G-runs for consensus motif detection (default: 3)")
    parser.add_argument("--min_loop_length", type=int, default=1, help="Minimum loop length between G-runs for consensus motif (default: 1)")
    parser.add_argument("--max_loop_length", type=int, default=7, help="Maximum loop length between G-runs for consensus motif (default: 7)")
    parser.add_argument("--regex_multiplicity", type=int, default=3, help="Number of G-run repeats required for consensus motif (default: 3)")
    parser.add_argument("--min_arm_length", type=int, default=10, help="Minimum arm length for direct repeat detection (default: 10)")
    parser.add_argument("--min_spacer_length", type=int, default=0, help="Minimum spacer length for direct repeats (default: 0)")
    parser.add_argument("--max_spacer_length", type=int, default=8, help="Maximum spacer length for direct repeats (default: 8)")
    parser.add_argument("--window_size", type=int, default=25, help="Sliding window size for G4Hunter scoring (default: 25)")
    parser.add_argument("--min_score", type=float, default=1.5, help="Minimum G4Hunter score threshold for sequence extraction (default: 1.5)")
    parser.add_argument("--mode", type=str, default="G4", help="Analysis mode (currently only G4 supported)", choices=["G4"])
    parser.add_argument("--parse_consensus", type=int, default=0, choices=[0, 1],
                        help="Parse consensus pG4 motifs using regex patterns (0=no, 1=yes, default: 0)")
    parser.add_argument("--merge_sequences", type=int, default=0, choices=[0, 1],
                        help="Merge overlapping sequences in output (0=no, 1=yes, default:)")
    parser.add_argument("--outdir", type=str, default="G4_results", help="Output directory for results files")
    args = parser.parse_args()
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(exist_ok=True)

    hunter = G4Hunter(min_grun_length=args.min_grun_length,
                      min_loop_length=args.min_loop_length,
                      max_loop_length=args.max_loop_length,
                      min_arm_length=args.min_arm_length,
                      min_spacer_length=args.min_spacer_length,
                      max_spacer_length=args.max_spacer_length,
                      regex_multiplicity=args.regex_multiplicity,
                      window_size=args.window_size,
                      min_score=args.min_score,
                      outdir=outdir)
    results_df = hunter.run(args.infile, 
                            parse_consensus=bool(args.parse_consensus),
                            merge_sequences=bool(args.merge_sequences))