import pandas as pd

def process_feature_counts(feature_counts_path):
    # Read the feature counts matrix from the provided file path.
    feature_counts = pd.read_csv(feature_counts_path, sep="\t", comment="#")
    
    # Clean up the column names:
    # Split by "/" and take the last part, then split by "." and take the first 2 elements and join them seperated by a dot.
    feature_counts.columns = (feature_counts.columns.str.split("/").str[-1].str.split(".").str[:2].str.join('.'))
    # feature_counts.columns.str.split("/").str[-1].str.split(".").str[0]
    
    # Identify columns where all values are zero.
    columns_to_drop = feature_counts.columns[(feature_counts == 0).all()]
    
    # Drop those columns from the DataFrame.
    feature_counts = feature_counts.drop(columns=columns_to_drop)
    
    # Print the list of dropped columns.
    print("Samples Dropped(due to no reads):", columns_to_drop.tolist())
    
    feature_counts = feature_counts.drop(columns=['Chr','Start','End', "Strand","Length"])
    return feature_counts
