# Visualisation threshold blast all vs all

Script to have a visual view of the distribution of the blast all vs all based on the cluster done.

This script could help to visualised if some groups exists inside the cluster form with previous threshold

Dependencies :
--------------

- Python >=3.9
   - tqdm 4.62.3
   - plotly 5.6.0
   - pandas 1.4.0
   - biopython 1.79 

Could be quickly install using conda and mamba

```bash
mamba create -n visual_threshold -c biopython -c conda-forge tqdm plotly pandas biopython
conda activate visual_threshold
```

Usage:
------

```
usage: visualisation_threshold_blast.py [-h] -f <folder> -i <file> -b <file> -css <file> [-o <OUTPUT_DIR>] [-filter <SIZE>]

Threshold helper

options:
  -h, --help            show this help message and exit

General input dataset options:
  -f <folder>, --fasta_folder <folder>
                        Path to the folder with all the fasta of each groups: one fasta per group
  -i <file>, --input_fasta <file>
                        The path to the fasta file that was used for the blast all vs all
  -b <file>, --blastfile <file>
                        Blast all vs all file
  -css <file>, --cssfile <file>
                        CSS file
  -o <OUTPUT_DIR>, --output <OUTPUT_DIR>
                        Using <OUTPUT_DIR> for output files (default: blastFile directory)
  -filter <SIZE>, --filter_length <SIZE>
                        Use a filter in the length of the alignment

```

