import pandas as pd
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse TSV file')
    parser.add_argument('input', type=str, help='Input TSV file')
    parser.add_argument('output', type=str, help='Output TSV file')
    args = parser.parse_args()

    # Load the TSV file into a DataFrame
    df = pd.read_csv(args.input, sep='\t')

    # Create a list to store the new rows
    new_rows = []

    # Iterate over each row in the DataFrame
    for _, row in df.iterrows():
        # If the "name" field is "compound", create two new rows
        if row['name'] == 'compound':
            new_row1 = row.copy()
            new_row1['entry2'] = row['value']
            new_rows.append(new_row1)

            new_row2 = row.copy()
            new_row2['entry1'] = row['entry2']
            new_row2['entry2'] = row['value']
            new_rows.append(new_row2)
        else:
            new_rows.append(row)

    # Create a new DataFrame from the list of new rows
    new_df = pd.concat(new_rows, axis=1).transpose()

    # Write the new DataFrame to a new TSV file
    new_df.to_csv(args.output, sep='\t', index=False)