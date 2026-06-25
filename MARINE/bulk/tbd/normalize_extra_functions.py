# def filter_marine_counts_with_regions(marine_counts_path, positions_bed_file_path, edit_type, strandedness):
#     """
#     Filters marine_counts DataFrame to retain only C>T edits within given BED regions.

#     Parameters:
#         marine_counts (pd.DataFrame): Input DataFrame containing 'feature_type', 'strand_conversion', 'position', 'contig'.
#         regions (pr.PyRanges): BED-like PyRanges object with regions to filter by.

#     Returns:
#         pd.DataFrame: Filtered marine_counts DataFrame.
#     """

#     marine_counts = pd.read_csv(marine_counts_path, sep="\t")
#     bed_file = pd.read_csv(positions_bed_file_path, sep="\t", header=None)
#     bed_file.columns = ['Chromosome', 'Start', 'End', 'Name', 'Length', 'Strand']
#     edit_type_complement = complement_edit(edit_type)

#     # edit_type_map = {
#     # "A>G": "T>C",
#     # "T>C": "A>G",
#     # "C>T": "G>A",
#     # "G>A": "C>T",
#     # "A>C": "T>G",
#     # "T>G": "A>C",
#     # "C>G": "G>C",
#     # "G>C": "C>G",
#     # "A>T": "T>A",
#     # "T>A": "A>T",
#     # "C>A": "G>T",
#     # "G>T": "C>A"}
    
#     # edit_type_complement = edit_type_map.get(edit_type)
#     edit_type_complement = complement_edit(edit_type)
    
#     if strandedness == 0:
#         filtered_df = marine_counts.query(f"(feature_strand == '+' and conversion == '{edit_type}') or (feature_strand == '-' and conversion == '{edit_type_complement}')").copy()
  
#     elif strandedness == 1 or strandedness == 2: # Here since Marine already converts the strand information to the correct feature strand based on the strandedness of the data we are using just the strand conversion of the data
#         filtered_df = marine_counts.query("strand_conversion == 'C>T'").copy()  # Filter for valid feature types and C>T strand conversions
        

#     # filtered_df = marine_counts.query("feature_type != '-1' and strand_conversion == 'C>T'").copy()  # Filter for valid feature types and C>T strand conversions
    
#     filtered_df["End"] = filtered_df["position"] # Prepare PyRanges input
#     filtered_df["Start"] = filtered_df["position"] - 1

#     # filtered_df = filtered_df.rename(columns={"contig": "Chromosome"})
#     filtered_df = filtered_df.rename(columns={"contig": "Chromosome",
#                                           'feature_strand': 'Strand', # This feature strand is taken for downstream analyses for filtering since it will give the correct strand information regardless of the strand
#                                           "feature_name":"Name"})
    
#     positions = pr.PyRanges(filtered_df[["Chromosome", "Start", "End","Name","Strand"]])  # Intersect with regions
#     regions = pr.PyRanges(bed_file[["Chromosome", "Start", "End","Name","Strand"]])
    
#     matched = positions.join(regions, strandedness= "same",suffix = "_feature_counts").df # None join is inner join,only keeps those which are in both rows
#     # keep only those regions where the name is the equal to Name_b
#     matched = matched[matched["Name"] == matched["Name_feature_counts"]]
#     # adding Name to the tuple to filter the final dataframe
#     matched_positions = set(matched[["Chromosome", "End","Name"]].apply(tuple, axis=1))  # Build list of (chrom, start) positions for filtering
#     filtered_df = filtered_df[filtered_df[["Chromosome", "position","Name"]].apply(tuple, axis=1).isin(matched_positions)]

#     # changing back the filtered_df to the original names to avoid any confusion in downstream steps
#     filtered_df = filtered_df.rename(columns={"Name": "feature_name",
#                                           "Strand": "feature_strand"})

#     return filtered_df