#! /usr/bin/env python

from __future__ import print_function
import os
import sys
import argparse
import csv
import tarfile
from urllib.request import urlopen
from shutil import copyfileobj

import ncbi_taxdump_utils
from sourmash import manifest as sm_manifest


def download_file(url, output_path, chunk_size=10000):
    "Download the files needed to create a lineage file"
    print(f"Downloading from {url} to {output_path}")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with urlopen(url) as response:
        total_size = response.getheader('Content-Length')
        if total_size is None:
            with open(output_path, 'wb') as out_file:
                copyfileobj(response, out_file)
            print("Download complete (no content length available).")
            return

        total_size = int(total_size)
        downloaded = 0

        with open(output_path, 'wb') as out_file:
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                out_file.write(chunk)
                downloaded += len(chunk)
                percent = downloaded * 100 / total_size
                sys.stdout.write(f"\rDownloaded: {downloaded / 1e6:.2f} MB / {total_size / 1e6:.2f} MB ({percent:.1f}%)")
                sys.stdout.flush()
        print("\nDownload complete.")


def download_taxdump(output_path='taxdump', chunk_size=10000):
    """Download and extract the NCBI taxdump archive using urllib."""
    url = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
    tarball_path = os.path.join(output_path, 'taxdump.tar.gz')

    os.makedirs(output_path, exist_ok=True)
    print(f"Downloading taxdump from {url} to {tarball_path}")
    with urlopen(url) as response:
        total_size = response.getheader('Content-Length')
        if total_size is None:
            with open(tarball_path, 'wb') as out_file:
                copyfileobj(response, out_file)
            print("Download complete (no content length available).")
            return

        total_size = int(total_size)
        downloaded = 0

        with open(tarball_path, 'wb') as out_file:
            while True:
                chunk = response.read(chunk_size)
                if not chunk:
                    break
                out_file.write(chunk)
                downloaded += len(chunk)
                percent = downloaded * 100 / total_size
                sys.stdout.write(f"\rDownloaded: {downloaded / 1e6:.2f} MB / {total_size / 1e6:.2f} MB ({percent:.1f}%)")
                sys.stdout.flush()
        print("\nDownload complete.")

    print(f"Extracting {tarball_path} to {output_path}")
    with tarfile.open(tarball_path, 'r:gz') as tar:
        tar.extractall(path=output_path)

    os.remove(tarball_path)

def get_suffix(name):
    ident = name.split(' ')[0]
    assert ident.startswith('GC')
    return ident[3:]

def get_suffix_no_version(name):
    suffix = get_suffix(name)
    assert '.' in suffix
    return suffix.split('.')[0]

def run_lineage_file(nodes_dmp, names_dmp, assembly_summary_files, output, manifest_filter=None, ictv=False):
    want_taxonomy = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']
    ictv_taxonomy = ['superkingdom', 'clade', 'subrealm', 'kingdom', 'subkingdom', 'phylum', 'subphylum', 'class', 'subclass', 'order', 'suborder', 'family', 'subfamily', 'genus', 'subgenus', 'species']

    #global want_taxonomy

    taxfoo = ncbi_taxdump_utils.NCBI_TaxonomyFoo()

    print(f"loading nodes file '{nodes_dmp}'")
    taxfoo.load_nodes_dmp(nodes_dmp)
    print(f"loading names file '{names_dmp}'")
    taxfoo.load_names_dmp(names_dmp)

    mf_ident = set()
    mf_ident_no_version = set()
    if manifest_filter:
        mf = sm_manifest.BaseCollectionManifest.load_from_filename(manifest_filter)
        print(mf)
        print(f"Loaded {len(mf.rows)} idents from manifest")
        for row in mf.rows:
            name = row['name']
            mf_ident.add(get_suffix(name))
            mf_ident_no_version.add(get_suffix_no_version(name))
        print(mf_ident)


    want_taxonomy = want_taxonomy if not ictv else ictv_taxonomy

    w = csv.writer(output)
    w.writerow(['ident', 'taxid'] + want_taxonomy)

    for filename in assembly_summary_files:
        print(f"reading assembly summary file from '{filename}'")
        r = csv.reader(open(filename, newline=""), delimiter='\t')

        count = 0
        for row in r:
            if not row: continue
            if row[0][0] == '#':
                continue


            acc = row[0]
            if mf_ident is not None and get_suffix(acc) not in mf_ident:
                continue
            #if filter_idents is not None and acc not in filter_idents:
            #    continue

            taxid = row[5]
            taxid = int(taxid)

            lin_dict = taxfoo.get_lineage_as_dict(taxid, want_taxonomy)
            if not lin_dict:
                print(f"WARNING: taxid {taxid} not in taxdump files. Producing empty lineage.")

            row = [acc, taxid]
            for rank in want_taxonomy:
                name = lin_dict.get(rank, '')
                row.append(name)

            w.writerow(row)

            count += 1

        print(f'output {count} lineages')


def main():
    p = argparse.ArgumentParser()
    p.add_argument('nodes_dmp', nargs='?')
    p.add_argument('names_dmp', nargs='?')
    p.add_argument('assembly_summary_files', nargs='*')
    p.add_argument('-m', '--manifest', nargs='?')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'))
    p.add_argument('--ictv', action='store_true')
    p.add_argument('--download', action='store_true', help='Download required scripts and taxonomy files')
    p.add_argument('--domain', nargs='?')
    args = p.parse_args()

    if args.download:
        # Download utility scripts
        download_file(
            f"https://ftp.ncbi.nlm.nih.gov/genomes/genbank/{args.domain}/assembly_summary.txt",
            f"data/assembly_summary.{args.domain}.txt"
        )
        download_taxdump("taxdump")
        print("Download complete.")

        run_lineage_file("taxdump/nodes.dmp", "taxdump/names.dmp", [f"data/assembly_summary.{args.domain}.txt"], args.output, manifest_filter=args.manifest, ictv=args.ictv)
    else:
        if not (args.nodes_dmp and args.names_dmp and args.assembly_summary_files and args.output):
            raise("Missing required arguments for lineage generation.")

        run_lineage_file(args.nodes_dmp, args.names_dmp, args.assembly_summary_files, args.output, manifest_filter=args.manifest, ictv=args.ictv)


if __name__ == '__main__':
    sys.exit(main())
