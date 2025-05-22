####
# make a taxa-specific sourmash database and pangenome
#
# run `watch 'ls -lathr dbs/'` to watch directsketch add to dbs
#
# snakemake -s directsketch.smk --profile slurm --resources allowed_jobs=100
####

import os
import csv

NCBI_API_KEY = os.environ.get("NCBI_API_KEY")

configfile: "config/config.yaml"

OUTDIR = config.get('output_directory')

# map various NAMEs to NCBI taxonomic ID
NAMES_TO_TAX_ID = config.get('taxon_name_and_id')

KSIZES = config.get('k_values')

# Dictionary for correct slurm batch allocations
PART_JOBS = {1: ['low2', 1], 2: ['low2', 1], 3: ['med2', 10], 4: ['med2', 10], 5: ['high2', 33]}
BIGMEM_JOBS = {1: ['bmm', 12], 2: ['bmm', 12], 3: ['bmm', 25], 4: ['bmh', 33], 5: ['bmh', 33]}

email = config.get('email')

if email:
    onsuccess:
        print("\nWorkflow finished without error\n")
        shell("mail -s 'Workflow finished without error' {email} < {log}")

    onerror:
        print("\nAn error occurred\n")
        shell("mail -s 'an error occurred' {email} < {log}")

if config.get('batch_size'):
    BATCHLIST_FILES = expand("{o}/data/{NAME}-all.batchlist.txt", o=OUTDIR, NAME=NAMES_TO_TAX_ID.keys())

    onsuccess:
        print("\nRemoving the database batch files\n")
        shell("cat {batchlists} | xargs rm -v", batchlists=" ".join(BATCHLIST_FILES))

wildcard_constraints:
    ksize = "\\d{2}",
    NAME = "\\w[^-.]+",

# psuedo rule to define the final output
rule all:
    input:
        #expand('{o}/data/{NAME}-all.csv', o=OUTDIR, NAME=set(NAMES_TO_TAX_ID)),
        #expand('{o}/data/{NAME}-ref.csv', o=OUTDIR, NAME=set(NAMES_TO_TAX_ID)),
        expand("{o}/dbs/{NAME}-ref-k{ksize}.zip", o=OUTDIR, NAME=set(NAMES_TO_TAX_ID), ksize = KSIZES),
        expand("{o}/dbs/{NAME}-k{ksize}.zip", o=OUTDIR, NAME=set(NAMES_TO_TAX_ID), ksize = KSIZES),

# target rule set
rule get_tax_info:
    input: 'scripts/1-get-by-tax.py',
    output:
        all = temporary('{o}/data/{NAME}-all.pickle'),
        ref = temporary('{o}/data/{NAME}-ref.pickle'),
    params:
        tax_id = lambda w: { **NAMES_TO_TAX_ID }.get(w.NAME)
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs= lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    conda:
        "envs/scripts.yaml",
    shell: """
        {input} --taxons {params.tax_id} -o {output.all} --all-genomes
        {input} --taxons {params.tax_id} -o {output.ref}
    """

checkpoint get_tax_csv:
    input:
        script = 'scripts/2-output-directsketch-csv.py',
        allpickle = '{o}/data/{NAME}-all.pickle',
        refpickle = '{o}/data/{NAME}-ref.pickle',
    output:
        all = '{o}/data/{NAME}-all.csv',
        ref = '{o}/data/{NAME}-ref.csv',
    params:
        tax_id = lambda w: { **NAMES_TO_TAX_ID }.get(w.NAME)
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs= lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    conda:
        "envs/scripts.yaml",
    shell: """
        {input.script} {input.allpickle} -o {output.all}
        {input.script} {input.refpickle} -o {output.ref}
    """

rule directsketch_reference:
    benchmark: '{o}/benchmark/directsketch-reference-{NAME}.tsv'
    input:
        ref = '{o}/data/{NAME}-ref.csv',
    output:
        failed = "{o}/data/directsketch-ref-{NAME}.failures.csv",
        ref = "{o}/dbs/{NAME}-ref.zip",
    conda: "envs/directsketch.yaml"
    threads: 16
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        disk_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    params:
        k_list = lambda wildcards: ",".join([f"k={ksize}" for ksize in KSIZES]),
        scale = config.get('scale_value'),
        api_key = NCBI_API_KEY,
        threads = lambda wildcards: int(10 if NCBI_API_KEY and NCBI_API_KEY.strip() else 3),
        log = "{o}/logs/directsketch-ref-{NAME}.log",
    shell:'''
        sourmash scripts gbsketch {input.ref} -o {output.ref} --failed {output.failed} \
            --param-string "dna,{params.k_list},scaled={params.scale},abund" \
            -a '{params.api_key}' -r 10 -n 30 -c {params.threads} -g 2> {params.log}
    '''

rule directsketch_all:
    benchmark: '{o}/benchmark/directsketch-all-{NAME}.tsv'
    input:
        all = '{o}/data/{NAME}-all.csv',
    output:
        all = "{o}/data/{NAME}-all.batchlist.txt",
    conda: "envs/directsketch.yaml"
    threads: 16
    resources:
        mem_mb = 96 * 1024,
        disk_mb = 96 * 1024,
        time = lambda wildcards, attempt: 72 * 60 * attempt,
        runtime = lambda wildcards, attempt: 72 * 60 * attempt,
        allowed_jobs=34,
        partition="bmh",
    params:
        k_list = lambda wildcards: ",".join([f"k={ksize}" for ksize in KSIZES]),
        scale = config.get('scale_value'),
        api_key = NCBI_API_KEY,
        threads = lambda wildcards: int(10 if NCBI_API_KEY and NCBI_API_KEY.strip() else 3),
        log = "{o}/logs/directsketch-all-{NAME}.log",
        failed = "{o}/data/directsketch-all-{NAME}.failures.csv",
        batch_size = config.get('batch_size'),
    shell:'''
        sourmash scripts gbsketch {input.all} -o {wildcards.o}/dbs/{wildcards.NAME}-all.zip --failed {params.failed} \
            --param-string "dna,{params.k_list},scaled={params.scale},abund" \
            -a '{params.api_key}' -r 10 -n 30 -c {params.threads} -g --batch-size {params.batch_size} \
            --allow-completed 2> {params.log}
        #mv {wildcards.o}/dbs/{wildcards.NAME}-all.zip.batches.txt {output.all}
    '''

#rule check_directsketch_all:
#    input: "{o}/data/{NAME}-all.done",
#    output:
#        sketch = "{o}/dbs/{NAME}-all.{n}.zip",
#    params:
#        k_list = KSIZES,
#        scale = config.get('scale_value'),
#    run:
#        from pathlib import Path
#        from sourmash import load_file_as_index, MinHash, signature
#        from sourmash.save_load import SaveSignaturesToLocation
#        
#        prev_file = Path(str(output.sketch))
#        
#        valid_sketch = False
#        try:
#            sigs = load_file_as_index(f"{output.sketch}")
#            print(sigs)
#            valid_sketch = True
#        except Exception as e:
#            print(f"Could not load sourmash signature: {e}")
#            valid_sketch = False
#        
#        if valid_sketch:
#            print(f"{prev_file} contains valid sourmash sketches. Skipping...")
#        else:
#            print(f"{prev_file} is missing or invalid. Creating empty sketch.")
#            with SaveSignaturesToLocation(output[0]) as save_sig:
#                for k in params.k_list:
#                    mh = MinHash(n=0, ksize=k, scaled=params.scale)
#                    ss = signature.SourmashSignature(mh, name=f'empty k{k} {output[0]}')
#                    save_sig.add(ss)
#
#
#class Checkpoint_MakePattern:
#    def __init__(self, pattern):
#        self.pattern = pattern
#
#    def get_batch_count(self, csv_file, batch_size):
#        with open(csv_file) as f:
#            reader = csv.reader(f)
#            next(reader)  # skip header
#            n_rows = sum(1 for row in reader)
#            n_batches = n_rows // batch_size
#            if n_rows % batch_size != 0:
#                n_batches += 1
#        return n_batches
#
#    def __call__(self, wildcards):
#        checkpoint_output = checkpoints.get_tax_csv.get(o=wildcards.o, NAME=wildcards.NAME)
#        csv_file = checkpoint_output.output.all
#
#        batch_size = config["batch_size"]
#        n_batches = self.get_batch_count(csv_file, batch_size)
#
#        return expand(self.pattern, n=range(1, n_batches + 1), **wildcards)
#
#rule cat_collect_all:
#    input:
#        #batches = Checkpoint_MakePattern("{o}/dbs/{NAME}-all.{n}.zip"),
#        batches = "{o}/data/{NAME}-all.batchlist.txt",
#    output:
#        collect = "{o}/dbs/{NAME}.zip",
#    conda: "envs/sourmash.yaml"
#    resources:
#        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
#        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
#        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
#        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
#        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
#    shell:'''
#        sourmash sig cat {input.batches} -o {output.collect}
#    '''

rule cat_sketch_ref:
    benchmark: '{o}/benchmark/cat_sketch-ref-{NAME}-{ksize}.tsv'
    input:
        ref = "{o}/dbs/{NAME}-ref.zip",
    output:
        ref = "{o}/dbs/{NAME}-ref-k{ksize}.zip",
    conda: "envs/sourmash.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        disk_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell: """
        sourmash sig cat {input.ref} -k {wildcards.ksize} -o {output.ref}
    """

rule cat_sketch_all:
    benchmark: '{o}/benchmark/cat_sketch-all-{NAME}-{ksize}.tsv'
    input:
        all = "{o}/data/{NAME}-all.batchlist.txt",
    output:
        all = "{o}/dbs/{NAME}-k{ksize}.zip",
    conda: "envs/sourmash.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: 42 * 1024 * attempt,
        disk_mb = lambda wildcards, attempt: 42 * 1024 * attempt,
        time = lambda wildcards, attempt: 12 * 60 * attempt,
        runtime = lambda wildcards, attempt: 12 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: 15, #BIGMEM_JOBS[attempt][1],
        partition=lambda wildcards, attempt: BIGMEM_JOBS[attempt][0],
    shell: """
        sourmash sig cat {input.all} -k {wildcards.ksize} -o {output.all}
    """
