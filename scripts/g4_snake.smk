import re 
import pandas as pd
from pathlib import Path
from pg4_extraction import G4Hunter

configfile: "config.yaml"

suffix = config["suffix"]
outdir = Path(config["outdir"]).resolve()
outdir.mkdir(parents=False, exist_ok=True)
files = pd.read_table(config["design"], 
                        usecols=["accession_id", "accession_path"],
                        index_col="accession_id"
                    )["accession_path"].to_dict()

rule all:
    input:
        expand(["%s/{accession}_pG4s.g4_hunter.tsv" % outdir,
                "%s/{accession}_pG4s.consensus.tsv" % outdir,
                "%s/{accession}_pG4s.merged.g4_hunter.tsv" % outdir
                ], 
                accession=list(files.keys())
                )

rule extract_pG4s:
    input:
        lambda wc: files[wc.accession]
    output:
        "%s/{accession}_pG4s.g4_hunter.tsv" % outdir,
        "%s/{accession}_pG4s.consensus.tsv" % outdir,
        "%s/{accession}_pG4s.merged.g4_hunter.tsv" % outdir,
    run:
        hunter = G4Hunter(outdir=outdir)
        hunter.run(input[0],
                   parse_consensus=True,
                   merge_sequences=True)