####
# make sure to load the apptainer module and activate snakemake conda env to run
#
# snakemake -s roary.smk --profile slurm --resources allowed_jobs=100 --use-singularity
####

import os
import csv

configfile: "./config/config.yaml"

# map various NAMEs to NCBI taxonomic ID
NAMES_TO_TAX_ID = config.get('taxon_name_and_id')

OUTDIR = [ config.get('output_directory') if config.get('output_directory') is not None else '../results' ]

tax_dict = {}
directory = './data/'
types = ['all']

for name in set(NAMES_TO_TAX_ID):
    tax_dict[name] = {}

    for type_ in types:
        filename = f'{name}-{type_}.csv'
        file_path = os.path.join(directory, filename)

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

for names, types in tax_dict.items():
    print(names, types)
    for type, value in types.items():
        print(names, type, value)


# List comprehension to generate the list
rule_all_list = [
    f'{name}-{type}/check/roary.done'
    for name, types in tax_dict.items()
    for type, sample_list in types.items()
    for sample in sample_list
]
print(rule_all_list)
# Dictionary for correct slurm batch allocations
PART_JOBS = {1: ['low2', 1], 2: ['low2', 1], 3: ['med2', 10], 4: ['med2', 10], 5: ['high2', 33]}
#PART_JOBS = {1: ['bml', 1], 2: ['bml', 1], 3: ['bmm', 25], 4: ['bmm', 25], 5: ['bmh', 33]}


# psuedo rule to define the final output
rule all:
    input:
        rule_all_list

# run roary
rule roary:
    output:
        check = "{name}-{type}/check/roary.done",
    container:
        "envs/roary.sif"
    threads: 24
    params:
        prokka_folder = f"{name}-{type}/prokka",
        gff_folder = f"{name}-{type}/prokka/gff",
        output_folder=f"{name}-{type}/roary"
    shell:
        """
        mkdir -p {params.gff_folder}
        find "{params.prokka_folder}" -type f -iname "*.gff" | while read file; do
            # Check if the file already exists in the destination folder
            dest="$PWD/fusobacterium_nucleatum-all/prokka/gff/$(basename "$file")"
            if [ "$file" != "$dest" ]; then
                cp "$file" "$dest"
            fi
        done    

        rm -rf {params.output_folder}
        roary -p {threads} -f {params.output_folder} -e -n -v {params.gff_folder}/*.gff && touch {output.check}
        """
