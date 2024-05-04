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
        # Split the "type", "value" and "name" fields by comma
        types = row['type'].split(',')
        values = row['value'].split(',')
        names = row['name'].split(',')

        # Create new rows for each value-name pair
        for type, value, name in zip(types, values, names):
            new_row = row.copy()
            new_row['type'] = type.strip()
            new_row['value'] = value.strip()
            new_row['name'] = name.strip()

            # Append the new row to the list
            new_rows.append(new_row)

    # Create a new DataFrame from the list of new rows
    new_df = pd.concat(new_rows, axis=1).transpose()


    # Write the new DataFrame to a new TSV file
    new_df.to_csv(args.output, sep='\t', index=False)
