#!/usr/bin/env bash

OUTFILE="master_counts.tsv"

FILES=$(ls mehreen_*_counts.txt | sort)

if [ -z "$FILES" ]; then
    echo "No count files found."
    exit 1
fi

echo "Merging:"
echo "$FILES"
echo ""

# Extract gene IDs from first file (skip header)
FIRST=$(echo "$FILES" | head -n 1)
tail -n +2 "$FIRST" | cut -f1 > genes.tmp

# Create header
echo -n "GeneID" > header.tmp
for f in $FILES; do
    SAMPLE=$(basename "$f" _counts.txt)
    echo -ne "\t$SAMPLE" >> header.tmp
done
echo "" >> header.tmp

# Extract counts (skip header)
for f in $FILES; do
    tail -n +2 "$f" | cut -f2 > "$f.tmp"
done

# Merge
paste genes.tmp mehreen_*_counts.txt.tmp > body.tmp
cat header.tmp body.tmp > "$OUTFILE"

# Cleanup
rm -f genes.tmp header.tmp body.tmp mehreen_*_counts.txt.tmp

echo "Clean master matrix written to $OUTFILE"
