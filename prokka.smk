####
# make sure to load the apptainer module and activate snakemake conda env to run
####

import os
import csv

API_KEY = os.environ.get('NCBI_API_KEY')

configfile: "./config/config.yaml"

OUTDIR = config['output_directory'] if config.get('output_directory') else './results'
print(OUTDIR)

# map various NAMEs to NCBI taxonomic ID
NAMES_TO_TAX_ID = config.get('taxon_name_and_id')

tax_dict = {}
directory = 'data'
types = ['all', 'ref']

for name in set(NAMES_TO_TAX_ID):
    tax_dict[name] = {}

    for type_ in types:
        filename = f'{name}-{type_}.csv'
        print(filename)
        file_path = os.path.join(OUTDIR, directory, filename)
        print(file_path)
        if os.path.exists(file_path):
            with open(file_path, mode='r', newline='') as fp:
                reader = csv.DictReader(fp)
                samples = []

                for row in reader:
                    samples.append(row['accession'])

                if config.get('test'):
                    tax_dict[name][type_] = samples[:config.get('test')]
                else:
                    tax_dict[name][type_] = samples

if config.get('test'):
    print('Only returning a max of', config.get('test'), 'samples:')
    print(tax_dict)

# List comprehension to generate the list
rule_all_list = [
    f'{OUTDIR}/{name}-{type_}/check/{sample}.prokka.done'
    for name, types in tax_dict.items()
    for type_, sample_list in types.items()
    for sample in sample_list
]

rule_fasta_list = [
    f'{OUTDIR}/{name}-{type_}/fasta/{sample}.fna'
    for name, types in tax_dict.items()
    for type_, sample_list in types.items()
    for sample in sample_list
]

print(NAMES_TO_TAX_ID, tax_dict.keys())

#print(rule_fasta_list)
#print(tax_dict)
# Dictionary for correct slurm batch allocations
#PART_JOBS = {1: ['low2', 1], 2: ['low2', 1], 3: ['med2', 10], 4: ['med2', 10], 5: ['high2', 33]}
PART_JOBS = {1: ['bml', 1], 2: ['bml', 1], 3: ['bmm', 25], 4: ['bmm', 25], 5: ['bmh', 33]}


# psuedo rule to define the final output 
rule all:
    input:
        rule_all_list

rule fasta:
    input:
        rule_fasta_list

# target rule set
rule download_genomes:
    input:
        '{o}/data/{NAME}-{type}.csv',
    output:
        acc = temporary('{o}/data/{NAME}-{type}-acc-list.txt'),
        zip = temporary('{o}/{NAME}-{type}.zip'),
    params:
        dir = '{o}/{NAME}-{type}',
    threads:
        10
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=10,
        partition='high2' #lambda wildcards, attempt: PART_JOBS[attempt][0],
    conda:
        'envs/prokka.yaml'
    shell: """
        # I should add a backup step that downloads from ENA if the NCBI fails for any reason...

        awk -F, 'NR>1 {{print $1}}' {input} > {output.acc}

        datasets download genome accession --dehydrated --inputfile {output.acc} --filename {output.zip} --api-key {API_KEY} #--no-progressbar

        unzip -o {output.zip} -d {params.dir}

        datasets rehydrate --api-key {API_KEY} --max-workers 30 --directory {params.dir}
       """

rule check_genomes:
    input:
        zip = '{o}/{NAME}-{type}.zip',
    output:
        fna = '{o}/{NAME}-{type}/fasta/{sample}.fna',
        #check = temporary('{o}/{NAME}-{type}/check/{sample}.md5sum'),
        #sum = temporary('{o}/{NAME}-{type}/check/{sample}-checksum.txt'),
        #log = '{o}/{NAME}-{type}/check/{sample}-checksum-log.txt',
    params:
        dir = '{o}/{NAME}-{type}',
    resources:
        mem_mb = lambda wildcards, attempt: 2 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=10,
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    shell: """ 

        ln -s {params.dir}/ncbi_dataset/data/{wildcards.sample}/{wildcards.sample}*genomic.fna {output.fna}
    """
        #grep '{wildcards.sample}' {params.dir}/md5sum.txt | awk '{{print $1}}' > {output.check}

        #md5sum {output.fna} | awk '{{print $1}}' > {output.sum}

        #if [[ -f "{output.check}" && -f "{output.sum}" ]]; then
        #    diff {output.check} {output.sum} &> {output.log}
        #    && echo
        #else
        #    echo "Missing input files for diff" > {output.log}
        #fi

# run prokka
rule prokka:
    input:
        fna = '{o}/{NAME}-{type}/fasta/{sample}.fna',
    output:
        check = "{o}/{NAME}-{type}/check/{sample}.prokka.done",
    conda: 
        "prokka"
    threads:
        8
    resources:
        mem_mb = lambda wildcards, attempt: 8 * 1024 * attempt,
        time = lambda wildcards, attempt: 1.5 * 60 * attempt,
        runtime = lambda wildcards, attempt: 1.5 * 60 * attempt,
        allowed_jobs=lambda wildcards, attempt: PART_JOBS[attempt][1],
        partition=lambda wildcards, attempt: PART_JOBS[attempt][0],
    params:
        output_folder="{o}/{NAME}-{type}/prokka/{sample}"
    shell:
        """ 
        prokka --kingdom Bacteria --outdir {params.output_folder} \
        --norrna --notrna --prefix {wildcards.sample} --force\
        --locustag {wildcards.sample} {input.fna} && touch {output.check}
        """

## run prokka on the gtdb sequences (can i put these rules together?)
#rule prokka_gtdb:
#    input:
#        fa_gtdb = f"{OUTPUT_DIR}/{pang_name_out}/MAGs/{{genome}}.fna.gz",
#    output:
#        check = f"{OUTPUT_DIR}/{pang_name_out}/check/{{genome}}.prokka.done",
#    conda: 
#        "prokka"
#    threads: 1
#    params:
#        output_folder=f"{OUTPUT_DIR}/{pang_name_out}/prokka/{{genome}}"
#    shell:
#        """ 
#        prokka --kingdom Bacteria --outdir {params.output_folder} \
#        --norrna --notrna --prefix {wildcards.genome} \
#        --locustag {wildcards.genome} {input.fa_gtdb} && touch {output.check}
#        """

