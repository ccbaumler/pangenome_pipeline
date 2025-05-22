import pandas as pd
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Process CSV to rename columns.")
parser.add_argument("input_file", help="Path to the input CSV file")
parser.add_argument("output_file", help="Path to the output CSV file")
args = parser.parse_args()

# Load the input CSV
df = pd.read_csv(args.input_file)

# Create new columns 'accession' and 'name' from 'ident'
df['accession'] = df['ident']
df['name'] = df['ident']

# Drop the 'ident' column if not needed
df = df[['accession', 'name']]


# Save the modified DataFrame to the output file
df.to_csv(args.output_file, index=False)

print(f"Processed CSV saved to {args.output_file}")
