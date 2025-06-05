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

    existing_batchlists = [fp for fp in BATCHLIST_FILES if os.path.exists(fp)]
    print(existing_batchlists)

    onsuccess:
        if not existing_batchlists:
            print("\nNo batchlist files found.\n")
        else:
            print("\nProcessing batchlist files...\n")

            for batchlist in existing_batchlists:
                print(f"\nReading batchlist: {batchlist}")
                with open(batchlist) as f:
                    for filepath in f:
                        filepath = filepath.strip()
                        if not filepath:
                            continue
                        if os.path.exists(filepath):
                            print(f"Deleting {filepath}")
                            os.remove(filepath)
                        else:
                            print(f"Already deleted {filepath}")

wildcard_constraints:
    ksize = "\\d{2}",
    NAME = "\\w[^-.]+",

def createStandaloneManifest(wildcards):
    file = dict()

    first_ksize = KSIZES[0] if isinstance(KSIZES, (list, tuple)) else KSIZES

    file["mf"] = f"{wildcards.o}/dbs/{wildcards.NAME}-k{first_ksize}.zip"
    print(file)

    return file

# psuedo rule to define the final output
rule all:
    input:
        expand("{o}/dbs/{NAME}-ref-k{ksize}.zip", o=OUTDIR, NAME=set(NAMES_TO_TAX_ID), ksize = KSIZES),
        expand("{o}/dbs/{NAME}-k{ksize}.zip", o=OUTDIR, NAME=set(NAMES_TO_TAX_ID), ksize = KSIZES),
        expand("{o}/dbs/{NAME}-k{ksize}.species.csv", o=OUTDIR, NAME=set(NAMES_TO_TAX_ID), ksize = KSIZES),


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
        failed = "{o}/data/directsketch-ref-{NAME}.download-failures.csv",
        check = "{o}/data/directsketch-ref-{NAME}.checksum-failures.csv",
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
        sourmash scripts gbsketch {input.ref} -o {output.ref} \
            --checksum-fail {output.check} --failed {output.failed} \
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
        failed = "{o}/data/directsketch-all-{NAME}.download-failures.csv",
        check = "{o}/data/directsketch-all-{NAME}.checksum-failures.csv",
        batch_size = config.get('batch_size'),
    shell:'''
        sourmash scripts gbsketch {input.all} -o {wildcards.o}/dbs/{wildcards.NAME}-all.zip \
            --checksum-fail {params.check} --failed {params.failed} \
            --param-string "dna,{params.k_list},scaled={params.scale},abund" \
            -a '{params.api_key}' -r 10 -n 30 -c {params.threads} -g --batch-size {params.batch_size} \
            --allow-completed 2> {params.log}
        mv {wildcards.o}/dbs/{wildcards.NAME}-all.zip.batchlist.txt {output.all}
    '''

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

rule collect_all:
    input:
        unpack(createStandaloneManifest), #unpack the databases and use only the first ksize to run this rule
    output:
        mf = "{o}/data/{NAME}.manifest.csv",
    conda: "envs/sourmash.yaml",
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell: """
        sourmash sig manifest --no-rebuild {input.mf} -o {output.mf}
    """

rule lineage_file:
    input:
        script = "scripts/make_lineage_csv.py",
        mf = "{o}/data/{NAME}.manifest.csv",
    output:
        taxa = "{o}/lineages.{NAME}.csv",
    conda: "envs/sourmash.yaml",
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    params:
        domain = config.get("lineage_domain"),
    shell: """
        {input.script} --download --domain {params.domain} -m {input.mf} -o {output.taxa}
    """



rule pangenome_database:
    input:
        db = "{o}/dbs/{NAME}-k{ksize}.zip",
        taxa = "{o}/lineages.{NAME}.csv",
    output:
        pan = "{o}/dbs/{NAME}-k{ksize}.species.zip",
    conda: "envs/pangenome.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 8 * 60 * attempt,
        runtime = lambda wildcards, attempt: 8 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    params:
        scale = config.get('scale_value'),
    shell: """
        sourmash scripts pangenome_createdb {input.db} -t {input.taxa} -a -k {wildcards.ksize} --scaled {params.scale} -o {output.pan}
    """

rule pangenome_ranktable:
    input:
        pan = "{o}/dbs/{NAME}-k{ksize}.species.zip",
    output:
        ranktable = "{o}/dbs/{NAME}-k{ksize}.species.csv",
    conda: "envs/pangenome.yaml"
    resources:
        mem_mb = lambda wildcards, attempt: 16 * 1024 * attempt,
        time = lambda wildcards, attempt: 8 * 60 * attempt,
        runtime = lambda wildcards, attempt: 8 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    params:
        scale = config.get('scale_value'),
    shell: """
        sourmash scripts pangenome_ranktable \
            {input.pan} \
            -o {output.ranktable} \
            -k {wildcards.ksize}
    """

