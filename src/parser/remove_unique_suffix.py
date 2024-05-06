import pandas as pd
import argparse
import click

@click.command()
@click.argument('input', type=click.Path(exists=True))
@click.option('--output', type=click.Path())
def remove_unique_suffix(input: str, output: str):
    if output is None:
        output = input.replace('.tsv', '_removed_suffix.tsv')

    # Load the TSV file into a DataFrame
    df = pd.read_csv(input, sep='\t')

    # Remove the suffix "-[0-9]+" in the columns 'entry1', 'entry2', and 'value'
    df['entry1'] = df['entry1'].str.replace('-[0-9]+$', '', regex=True)
    df['entry2'] = df['entry2'].str.replace('-[0-9]+$', '', regex=True)
    df['value'] = df['value'].str.replace('-[0-9]+$', '', regex=True)

    # Write the new DataFrame to a new TSV file
    df.to_csv(output, sep='\t', index=False)


if __name__ == '__main__':
    remove_unique_suffix()
