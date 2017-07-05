# Eval-insertfill

Evaluates insertion genotyping tools.

## Installing

0. Clone this repository. Note: contains git submodules, clone recursively.
1. Install all the tools (There are links provided in *configuration.sh*).
2. Make sure *configuration.sh* points to the binaries.
3. Get a reference genome and point to it in *configuration.sh*.

We used chr17 from GRCh38 and removed all Ns with *scripts/fix_chr17.py*.

## Running

The experiments are given as shell scripts, which run and analyze the results for the experiments.
