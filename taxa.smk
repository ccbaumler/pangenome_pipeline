configfile: "./config/config.yaml"

OUTDIR = config.get('output_directory')

# map various NAMEs to NCBI taxonomic ID
NAMES_TO_TAX_ID = config.get('taxon_name_and_id')

# Dictionary for correct slurm batch allocations
PART_JOBS = {1: ['low2', 1], 2: ['low2', 1], 3: ['med2', 10], 4: ['med2', 10], 5: ['high2', 33]}
#PART_JOBS = {1: ['bml', 1], 2: ['bml', 1], 3: ['bmm', 25], 4: ['bmm', 25], 5: ['bmh', 33]}


# psuedo rule to define the final output
rule all:
    input:
        expand('{o}/data/{NAME}-all.csv', o=OUTDIR, NAME=set(NAMES_TO_TAX_ID)),
        expand('{o}/data/{NAME}-ref.csv', o=OUTDIR, NAME=set(NAMES_TO_TAX_ID)),
 
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

rule get_tax_csv:
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

