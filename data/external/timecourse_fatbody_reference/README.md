# Gregaria time-course fat body single-pair reference

This folder is reserved for the two non-`FULL` *S. gregaria* fat body samples
from the `locust-time-course-RNAseq` project.

These samples are useful only as a qualitative tissue-matched phase reference,
because there is one solitarious/isolated and one gregarious/crowded fat body
sample. They should not be used for formal DEG testing.

Expected featureCounts files:

- `readcounts/GREG_G_CCT_11_FAT_featurecounts.txt`
- `readcounts/GREG_S_ICT_9_FAT_featurecounts.txt`

The main time-course metadata file lists `GREG_S_ICT_9_FAT_MERGE` with
`Phase = gregarious`, but the sample name and derived time-course metadata
identify it as the isolated/solitarious control. The local metadata in this
folder uses the corrected biological label.
