#!/bin/bash

OUTDIR="PANGENOMES"
mkdir -p $OUTDIR

while read -r line; do
    wget -P $OUTDIR/ "$line"
done < ftp_links.txt