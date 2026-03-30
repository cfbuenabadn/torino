# Functions for processing s-NMF factors and obtain isoform structures
# Originally written manually in R, translated into Python both manually and with the help of ChatGPT.
# Functions commented for clarity with the help of ChatGPT.

import numpy as np
import pandas as pd
import tabix
from pybedtools import BedTool
from scipy.stats import linregress

import rpy2
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

def run_tabix_on_junc(junc_file, coords, transform_to_df = True):
    chrom, start = coords[0].split(':')
    start = int(start)
    end = int(coords[-1].split(':')[1])
    
    tb = tabix.open(junc_file)
#     print(chrom, start, end)
    juncs = tb.query(chrom, start, end)

    if transform_to_df:

        junctions_bed = []#pd.DataFrame()
        junctions_bed_cols = []
        
        # Iterate over queried junctions and store them in the DataFrame
        for idx, record in enumerate(juncs):
            junctions_bed.append(pd.Series(record))
            junctions_bed_cols.append(f'record{str(idx)}')
            # junctions_bed[f'record{str(idx)}'] = record
        # print(record)
        # print(junctions_bed)
        if len(junctions_bed) == 0:
            return pd.DataFrame(columns = ['chrom', 'start', 'end', 'gene'])
        junctions_bed = pd.concat(junctions_bed, axis=1)
        junctions_bed.columns = junctions_bed_cols
        junctions_bed = junctions_bed.T
        return junctions_bed

    else:
        return juncs

def get_junctions_bed(coords, junc_file, gene_id, junc_thres=None):
    """
    Function: run_tabix
    
    Description:
    This function retrieves junction coordinates from a specified genomic region 
    using Tabix, and filters the junctions for a specific gene ID. It then formats 
    the retrieved junctions into a pandas DataFrame and converts them to a 
    BedTool object.
    
    Input:
    - coords: A list of genomic coordinates in the format 'chromosome:start-end'.
    - junc_file: A file path to a Tabix-indexed file containing junction information.
    - gene_id: A string representing the gene ID to filter the junctions.
    
    Output:
    - junctions_bed: A BedTool object containing junction coordinates for the specified gene ID, 
                     or an empty DataFrame if no junctions are found.
    """

    
    # Extract chromosome and start position from the first coordinate
    # chrom, start = coords[0].split(':')
    # start = int(start)
    # end = int(coords[-1].split(':')[1])
    
    # # Open Tabix file and query junctions within the specified genomic region
    # tb = tabix.open(junc_file)
    # juncs = tb.query(chrom, start, end)
    
    # # Create an empty DataFrame to store junctions
    # junctions_bed = pd.DataFrame()
    
    # # Iterate over queried junctions and store them in the DataFrame
    # for idx, record in enumerate(juncs):
    #     junctions_bed[f'record{str(idx)}'] = record
        
    # # Transpose the DataFrame and process gene IDs
    # junctions_bed = junctions_bed.T
    junctions_bed = run_tabix_on_junc(junc_file, coords)
    if junctions_bed.shape[0] == 0:
        junctions_bed = BedTool.from_dataframe(junctions_bed)
        return junctions_bed
    junctions_bed[3] = [x.split('.')[0] for x in junctions_bed[3]]

    junctions_bed_counts = junctions_bed.loc[junctions_bed[3] == gene_id, junctions_bed.columns[4:]].astype(int).quantile(0.999,axis=1)

    if junc_thres is None:
        junc_thres = junctions_bed_counts.max()/200
    
    # Filter junctions for the specified gene ID
    junctions_bed = junctions_bed.loc[junctions_bed[3] == gene_id, range(4)]
    junctions_bed.columns = ['chrom', 'start', 'end', 'gene']

    junctions_bed = junctions_bed[junctions_bed_counts >= junc_thres]
    
    # Return empty DataFrame if no junctions are found for the gene ID
    if junctions_bed.shape[0] == 0:
        junctions_bed = BedTool.from_dataframe(junctions_bed)
        return junctions_bed
    
    # Convert start and end positions to integers
    junctions_bed.start = junctions_bed.start.astype(int)
    junctions_bed.end = junctions_bed.end.astype(int)
    
    # Generate unique junction names
    junctions_bed['junc_names'] = ['junction_' + str(i) for i in range(1, junctions_bed.shape[0]+1)]
    
    # Rearrange DataFrame columns
    junctions_bed = junctions_bed[['chrom', 'start', 'end', 'junc_names']]
    
    # Convert DataFrame to BedTool object
    junctions_bed = BedTool.from_dataframe(junctions_bed)
    
    return junctions_bed

def factor_lm(factor, strand=None):
    """
    Function: factor_lm
    
    Description:
    This function performs a linear regression-based correction on a factor. 
    It scales the factor values by a quantile normalization to ensure that they 
    fall within a specific range. Then, it calculates the slope and intercept 
    of a linear regression model using the corrected factor values. 
    Based on the strand information provided, it applies different correction 
    methods. If the strand is not specified or if the slope of the regression 
    line matches the strand direction, it corrects the factor using the linear 
    regression model. Otherwise, it returns the original factor values.
    
    Input:
    - factor: An array-like object representing the factor to be corrected.
    - strand: A string indicating the strand direction ('plus' or 'minus'). 
              If None, the function corrects the factor regardless of strand.
    
    Output:
    - y_corrected: An array-like object containing the corrected factor values.
    
    Dependencies:
    - numpy: Numerical computing library for Python.
    - scipy.stats.linregress: Function for performing linear regression.
    """
    # Generate x-axis values
    x = np.arange(1, len(factor) + 1)
    # Convert factor to numpy array
    y = np.array(factor)
    # Normalize factor values using quantile normalization
    y /= np.quantile(y, 0.99)
    
    # Filter out extreme values
    y_ = y[y > 0.25]
    x_ = x[y > 0.25]
    
    # Perform linear regression
    slope, intercept, _, _, _ = linregress(x_, y_)
    
    # Correct factor based on strand information
    if strand is None:
        # Apply correction using linear regression model
        y_model = x * slope + intercept
        y_corrected = y / np.maximum(y_model, 0.1)
        y_corrected = np.where(y > 0.1, y_corrected, y)
    else:
        # Check if the strand direction matches the slope of the regression line
        if ((((strand == 'plus') or (strand == '+')) and slope > 0) or (((strand == 'minus') or (strand =='-')) and slope < 0)):
            # Apply correction using linear regression model
            y_model = x * slope + intercept
            y_corrected = y / np.maximum(y_model, 0.1)
            y_corrected = np.where(y > 0.1, y_corrected, y)
        else:
            # Return original factor values if strand direction does not match regression slope
            y_corrected = factor
    
    return y_corrected




def binarize_factor(factor, cutoff_=0.25, binary_fraction=0.5, step_fraction=10, exon_quant=0.99, cutoff_strict=0.1, pos_window=100):
    """
    Function: binarize_factor
    
    Description:
    This function binarizes a given factor, converting continuous values into binary 
    values based on specified parameters. It aims to identify peaks or significant 
    values in the factor. The binarization process involves scaling the factor values 
    by a quantile normalization to ensure they fall within a specific range. Then, 
    it iterates through the factor values, determining whether each value exceeds 
    a certain threshold and whether it represents a peak relative to neighboring values. 
    Based on these criteria, it assigns binary values to each position in the factor.
    
    Input:
    - factor: An array-like object representing the factor to be binarized.
    - cutoff_: A float representing the initial cutoff threshold for binarization.
    - binary_fraction: A float representing the fraction of the exon peak used as 
                       the threshold for binarization.
    - step_fraction: An integer representing the fraction of the factor length 
                     used for calculating the step size.
    - exon_quant: A float representing the quantile used for normalization.
    - cutoff_strict: A float representing the strict cutoff threshold for binarization.
    - pos_window: An integer representing the window size for considering neighboring 
                  positions during binarization.
    
    Output:
    - factor_binary: An array-like object containing the binarized factor values.
    
    Dependencies:
    - numpy: Numerical computing library for Python.
    """

    # Initialize parameters
    cutoff = cutoff_
    # Normalize factor values using quantile normalization
    factor = factor / np.quantile(factor, exon_quant)
    # Calculate step size
    step = len(factor) / step_fraction

    # Initialize array to store binarized factor values
    factor_binary = np.zeros(len(factor))

    # Initialize variable to store last peak position
    last_pos = None

    # Iterate through factor values
    for idx in range(len(factor)):
        f = factor[idx]

        # Apply strict cutoff
        if f <= cutoff_strict:
            factor_binary[idx] = 0
            last_pos = None
            continue

        # Determine binary value based on position and previous values
        if idx == 0:
            if f > cutoff:
                f_binary = 1
                last_pos = 0
            else:
                f_binary = 0
                last_pos = None
        else:
            if factor_binary[idx - 1] == 0:
                if f < cutoff:
                    f_binary = 0
                    last_pos = None
                else:
                    f_binary = 1
                    last_pos = idx
            else:
                # Update last peak position
                last_pos = max(last_pos, (idx - pos_window))
                # Calculate exon peak
                f_exon_peak = np.quantile(factor[last_pos:idx+1], exon_quant)

                # Determine binary value based on exon peak and cutoff
                if f < (f_exon_peak * binary_fraction) or (f < cutoff_strict):
                    f_binary = 0
                    cutoff = max(f, cutoff_, cutoff_strict)
                    last_pos = None
                else:
                    f_binary = 1

        factor_binary[idx] = f_binary

    # Store initial binary values
    first_binary = factor_binary.copy()

    # Iterate through factor values in reverse
    for idx in range(len(factor)-1, -1, -1):
        f = factor[idx]

        if idx == len(factor) - 1:
            # Check if the last position represents a peak
            if factor_binary[idx] == 1:
                last_pos = len(factor) - 1
            else:
                continue
        else:
            # Check neighboring values to identify peaks
            previous_f = factor_binary[idx + 1]
            current_f = factor_binary[idx]

            if previous_f == 0 and current_f == 1:
                last_pos = idx
                continue
            elif previous_f == 1:
                last_pos = min(last_pos, (idx + pos_window))

                if (current_f == 0) and ((last_pos - idx) < 3):
                    # Remove isolated peaks
                    factor_binary[idx:(last_pos + 1)] = 0
                    continue

                # Calculate exon peak
                f_exon_peak = np.quantile(factor[(idx + 1):(last_pos + 1)], exon_quant)

                # Determine binary value based on exon peak and cutoff
                if f > (f_exon_peak * binary_fraction) and (f > cutoff_strict):
                    factor_binary[idx] = 1
                else:
                    continue

    return factor_binary

def get_factor_gaps(df):
    """
    Function: get_factor_gaps
    
    Description:
    This function extracts the gaps between consecutive genomic regions 
    represented as DataFrame rows. It creates a new DataFrame containing 
    the intervals between the end of one region and the start of the next 
    region along with a unique identifier for each gap.
    
    Input:
    - df: A pandas DataFrame containing genomic regions represented by 
          'chrom', 'start', and 'end' columns.
    
    Output:
    - factor_junctions_bed: A pandas DataFrame representing the gaps 
                            between consecutive genomic regions, with 
                            columns 'chrom', 'start', 'end', and 'names'.
    
    Dependencies:
    - pandas: Data manipulation library for Python.
    """
    # Create a new DataFrame with one less row to represent gaps
    new_df = pd.DataFrame({
        'chrom': list(df['chrom'][1:]),  # Chromosome from the second row onwards
        'start': list(df['end'][:-1]),   # End of the previous region
        'end': list(df['start'][1:])     # Start of the next region
    })

    # Generate unique names for each gap
    new_df['names'] = ['factor_junction' + str(i) for i in range(1, new_df.shape[0] + 1)]
    # Rename the DataFrame for clarity
    factor_junctions_bed = new_df

    return factor_junctions_bed



def find_continuous_segments(coords, factor):
    """
    Function: find_continuous_segments
    
    Description:
    This function identifies continuous segments in genomic regions based on a binary factor 
    indicating the presence or absence of certain features. It iterates through the factor 
    values to detect transitions between 0 and 1, representing the start and end of segments. 
    It returns a DataFrame containing the start and end coordinates of each identified segment.
    
    Input:
    - coords: A list of genomic coordinates indicating the start and end of each segment.
    - factor: An array-like object representing genomic factors associated with each segment.
    
    Output:
    - segments_df: A pandas DataFrame containing the identified continuous segments with 
                   columns 'chrom', 'start', and 'end'.
    
    Dependencies:
    - pandas: Data manipulation library for Python.
    """
    # Initialize variables
    segments_df = pd.DataFrame(columns=['chrom', 'start', 'end'])
    in_segment = False
    start_idx = None
    current_chrom = None
    
    # Iterate through factor values
    for i in range(len(factor)):
        # Detect segment start
        if factor[i] == 1 and not in_segment:
            in_segment = True
            start_idx = int(coords[i].split(':')[1])
        # Detect segment end
        elif factor[i] == 0 and in_segment:
            in_segment = False
            end_idx = int(coords[i - 1].split(':')[1])
            # Add segment to DataFrame
            segment = pd.DataFrame({'chrom': [current_chrom], 'start': [start_idx], 'end': [end_idx]})
            segments_df = pd.concat([segments_df, segment], ignore_index=True)
        
        # Capture chromosome information
        if in_segment and current_chrom is None:
            current_chrom = coords[i].split(':')[0]
    
    # Check if the last segment is still continuing till the end
    if in_segment:
        end_idx = int(coords[-1].split(':')[1])
        segment = pd.DataFrame({'chrom': [current_chrom], 'start': [start_idx], 'end': [end_idx]})
        segments_df = pd.concat([segments_df, segment], ignore_index=True)
    
    return segments_df



# def get_corrected_exons(segments_df, factor_junctions_bed, bed_intersection, junctions_bed, coords):
#     """
#     Function: get_corrected_exons
    
#     Description:
#     This function corrects exon boundaries based on junction information and the gaps between 
#     continuous segments. It iterates through the 'factor_junctions_bed' DataFrame to find 
#     the best junction for each gap and corrects the exon boundaries accordingly. It returns 
#     a DataFrame containing the corrected exon coordinates.
    
#     Input:
#     - segments_df: A pandas DataFrame containing the identified continuous segments with 
#                    columns 'chrom', 'start', and 'end'.
#     - factor_junctions_bed: A pandas DataFrame containing junction information for the factors.
#     - bed_intersection: A pandas DataFrame representing the intersection between genomic 
#                         regions and junctions.
#     - junctions_bed: A pandas DataFrame containing information about genomic junctions.
#     - coords: A list of genomic coordinates indicating the start and end of each segment.
    
#     Output:
#     - corrected_exons: A pandas DataFrame containing the corrected exon coordinates with 
#                         columns 'chrom', 'start', and 'end'.
    
#     Dependencies:
#     - pandas: Data manipulation library for Python.
#     """
#     # Calculate gene size
#     gene_size = int(coords[-1].split(':')[1]) - int(coords[0].split(':')[1])
#     corrected_exons = pd.DataFrame(columns=['chrom', 'start', 'end'])
    
#     # Initialize lists to store exon boundaries
#     start_list = [segments_df['start'][0]]
#     end_list = []
#     # Iterate through factor junctions
#     for i in range(len(factor_junctions_bed['start'])):
#         factor_junc = factor_junctions_bed['names'][i]
#         # print(factor_junc)
#         # print(bed_intersection)
#         bed_intersection_slice = bed_intersection[bed_intersection['names'] == factor_junc]
#         # print(bed_intersection_slice)
#         # Find the best junction for the current factor junction
#         best_junction = find_best_junction(bed_intersection_slice)
    
#         if best_junction is None:
#             continue
#         elif best_junction == 'gap':
#             # Handle gap
#             end_list.append(int(factor_junctions_bed['start'][i]))
#             start_list.append(int(factor_junctions_bed['end'][i]))
#         else:
#             # Retrieve junction coordinates
#             end_list.append(junctions_bed[junctions_bed['junc_names'] == best_junction]['start'].values[0])
#             start_list.append(junctions_bed[junctions_bed['junc_names'] == best_junction]['end'].values[0])
    
#     # Add the last segment end coordinate
#     end_list.append(segments_df['end'].iloc[-1])

#     # Create DataFrame for corrected exons
#     for i in range(len(start_list)):
#         df_row = pd.DataFrame({'chrom': [segments_df['chrom'][0]], 'start': [start_list[i]], 'end': [end_list[i]]})
#         corrected_exons = pd.concat([corrected_exons, df_row], ignore_index=True)

#     ### Add correction for negative exons

#     corrected_exons.start += 1
#     corrected_exons.end -= 1

#     corrected_exons_len = corrected_exons.end - corrected_exons.start

#     corrected_exons = corrected_exons.loc[corrected_exons_len >= 3]
        
#     return corrected_exons



# def find_best_junction(bed_intersection_slice, gene_size=10000):
#     """
#     Function: find_best_junction
    
#     Description:
#     This function determines the best junction for a given gap between continuous segments based 
#     on intersection information and junction characteristics. It calculates various metrics such 
#     as distance and overlap percentage to select the most suitable junction. If no suitable junction 
#     is found, it returns 'gap' to indicate the presence of a gap.
    
#     Input:
#     - bed_intersection_slice: A DataFrame slice representing the intersection between genomic regions 
#                               and junctions for a specific gap.
#     - gene_size: An integer representing the size of the gene (default value: 10000).
    
#     Output:
#     - best_junc: A string representing the best junction found for the given gap.
    
#     Dependencies:
#     None.
#     """
#     # Initialize variables
#     gap_start = int(bed_intersection_slice['start'].iloc[0])
#     gap_end = int(bed_intersection_slice['end'].iloc[0])
#     gap_size = gap_end - gap_start
    
#     min_distance = gene_size/50
#     min_distance = np.min([min_distance, 200])
#     min_distance = np.max([min_distance, 30])
#     # print(min_distance)
    
#     best_junc = None
#     any_close = False
#     any_start_close = False
#     any_end_close = False

#     # Check if the slice is empty
#     if len(bed_intersection_slice) == 0:
#         if gap_size >= gene_size / 10:
#             best_junc = "gap"
#             return best_junc
#         else:
#             return None
#     # Check if the gap size indicates a gap
#     elif bed_intersection_slice['overlap'].iloc[0] == '0':
#         gap_min_size = gene_size / 10
#         if gap_size >= gap_min_size:
#             best_junc = "gap"
#     else:
#         # Initialize variables for best junction selection
#         best_distance = 1000
#         #best_percent = 0
        
#         # Extract unique junctions in the intersection slice
#         overlapping_junctions = bed_intersection_slice['junc_names'].unique()
#         for junc in overlapping_junctions:
#             bed_intersection_junc = bed_intersection_slice[bed_intersection_slice['junc_names'] == junc]
#             junc_start = int(bed_intersection_junc['junc_start'].iloc[0])
#             junc_end = int(bed_intersection_junc['junc_end'].iloc[0])
#             junc_size = int(bed_intersection_junc['overlap'].iloc[0])
#             # print(junc_size)

#             if junc_size < 1:
#                 continue
#                 # print(bed_intersection_junc)

#             # Calculate distances and percentages
#             start_close = abs(gap_start - junc_start)
#             end_close = abs(gap_end - junc_end)
#             gap_percent = junc_size / gap_size
#             junc_percent = gap_size / junc_size
#             current_distance = start_close + end_close

#             # Update flags for close distances
#             if start_close <= min_distance:
#                 any_start_close = True
#             if end_close <= min_distance:
#                 any_end_close = True

#             # Update flag for any close distance
#             if start_close <= min_distance or end_close <= min_distance:
#                 any_close = True

#             # Select best junction based on criteria
#             if (start_close <= min_distance) and (end_close <= min_distance) and (gap_percent >= 0.75) and (junc_percent > 0.75):
#                 if (current_distance < best_distance):# and (gap_percent > best_percent):
#                     best_junc = junc
#                     best_distance = current_distance
#                     #best_percent = gap_percent
            
#     # Handle cases where no suitable junction is found
#     if best_junc is None:
#         gap_min_size = gene_size / 10
#         if gap_size > 1000:
#             best_junc = 'gap'
#         elif any_start_close and any_end_close:
#             best_junc = 'gap'
            
#     return best_junc


def get_corrected_exons(segments_df, factor_junctions_bed, bed_intersection, junctions_bed, coords):
    """
    Function: get_corrected_exons
    
    Description:
    This function corrects exon boundaries based on junction information and the gaps between 
    continuous segments. It iterates through the 'factor_junctions_bed' DataFrame to find 
    the best junction for each gap and corrects the exon boundaries accordingly. It returns 
    a DataFrame containing the corrected exon coordinates.
    
    Input:
    - segments_df: A pandas DataFrame containing the identified continuous segments with 
                   columns 'chrom', 'start', and 'end'.
    - factor_junctions_bed: A pandas DataFrame containing junction information for the factors.
    - bed_intersection: A pandas DataFrame representing the intersection between genomic 
                        regions and junctions.
    - junctions_bed: A pandas DataFrame containing information about genomic junctions.
    - coords: A list of genomic coordinates indicating the start and end of each segment.
    
    Output:
    - corrected_exons: A pandas DataFrame containing the corrected exon coordinates with 
                        columns 'chrom', 'start', and 'end'.
    
    Dependencies:
    - pandas: Data manipulation library for Python.
    """
    # Calculate gene size
    gene_size = int(coords[-1].split(':')[1]) - int(coords[0].split(':')[1])
    corrected_exons = pd.DataFrame(columns=['chrom', 'start', 'end'])
    
    # Initialize lists to store exon boundaries
    start_list = [segments_df['start'].iloc[0]]
    end_list = []
    # Iterate through factor junctions
    for i in range(len(factor_junctions_bed['start'])):
        factor_junc = factor_junctions_bed['names'][i]
        if i > 0:
            start_limit = int(list(factor_junctions_bed.end)[i-1])
        else:
            start_limit = -1
        if i < (len(factor_junctions_bed['start']) - 1):
            end_limit = int(list(factor_junctions_bed.start)[i+1])
        else:
            end_limit = 1e100
        
        # print(factor_junc)
        # print(bed_intersection)
        bed_intersection_slice = bed_intersection[bed_intersection['names'] == factor_junc]
        # print(bed_intersection_slice)
        # Find the best junction for the current factor junction
        best_junction = find_best_junction(bed_intersection_slice, start_limit, end_limit, gene_size=gene_size)

        print(best_junction)

        # print(best_junction)
    
        if best_junction is None:
            continue
        elif best_junction[:3] == 'gap':
            gap, start_gap, end_gap = best_junction.split(':')
            if start_gap == '':
                end_list.append(int(factor_junctions_bed['start'][i]))
            else:
                end_list.append(int(start_gap))
            if end_gap == '':
                start_list.append(int(factor_junctions_bed['end'][i]))
            else:
                start_list.append(int(end_gap))
            # Handle gap
            
            
        else:
            # Retrieve junction coordinates
            end_list.append(junctions_bed[junctions_bed['junc_names'] == best_junction]['start'].values[0])
            start_list.append(junctions_bed[junctions_bed['junc_names'] == best_junction]['end'].values[0])
    
    # Add the last segment end coordinate
    end_list.append(segments_df['end'].iloc[-1])

    # Create DataFrame for corrected exons
    for i in range(len(start_list)):
        df_row = pd.DataFrame({'chrom': [segments_df['chrom'][0]], 'start': [start_list[i]], 'end': [end_list[i]]})
        corrected_exons = pd.concat([corrected_exons, df_row], ignore_index=True)

    ### Add correction for negative exons

    corrected_exons.start += 1
    corrected_exons.end -= 1

    corrected_exons_len = corrected_exons.end - corrected_exons.start

    corrected_exons = corrected_exons.loc[corrected_exons_len >= 3]
        
    return corrected_exons



def find_best_junction(bed_intersection_slice, start_limit, end_limit, gene_size=10000):
    """
    Function: find_best_junction
    
    Description:
    This function determines the best junction for a given gap between continuous segments based 
    on intersection information and junction characteristics. It calculates various metrics such 
    as distance and overlap percentage to select the most suitable junction. If no suitable junction 
    is found, it returns 'gap' to indicate the presence of a gap.
    
    Input:
    - bed_intersection_slice: A DataFrame slice representing the intersection between genomic regions 
                              and junctions for a specific gap.
    - gene_size: An integer representing the size of the gene (default value: 10000).
    
    Output:
    - best_junc: A string representing the best junction found for the given gap.
    
    Dependencies:
    None.
    """
    # Initialize variables
    gap_start = int(bed_intersection_slice['start'].iloc[0])
    gap_end = int(bed_intersection_slice['end'].iloc[0])
    gap_size = gap_end - gap_start
    
    min_distance = gene_size/50
    min_distance = np.min([min_distance, 200])
    min_distance = np.max([min_distance, 30])
    min_distance_start = min_distance
    min_distance_end = min_distance
    closest_start = 0
    closest_end = 0
    # print(min_distance)
    
    best_junc = None
    # any_close = False
    any_start_close = False
    any_end_close = False

    # Check if the slice is empty
    if len(bed_intersection_slice) == 0:
        if gap_size >= gene_size / 10:
            best_junc = "gap::"
            return best_junc
        else:
            return None
    # Check if the gap size indicates a gap
    elif bed_intersection_slice['overlap'].iloc[0] == '0':
        gap_min_size = gene_size / 10
        if gap_size >= gap_min_size:
            
            best_junc = "gap::"
    else:
        # Initialize variables for best junction selection
        best_distance = 1000
        #best_percent = 0
        
        # Extract unique junctions in the intersection slice
        overlapping_junctions = bed_intersection_slice['junc_names'].unique()
        for junc in overlapping_junctions:
            bed_intersection_junc = bed_intersection_slice[bed_intersection_slice['junc_names'] == junc]
            junc_start = int(bed_intersection_junc['junc_start'].iloc[0])
            junc_end = int(bed_intersection_junc['junc_end'].iloc[0])
            junc_size = int(bed_intersection_junc['overlap'].iloc[0])
            # print(junc_size)

            if junc_size < 1:
                continue
                # print(bed_intersection_junc)

            # Calculate distances and percentages
            start_close = abs(gap_start - junc_start)
            end_close = abs(gap_end - junc_end)
            gap_percent = junc_size / gap_size
            junc_percent = gap_size / junc_size
            current_distance = start_close + end_close

            # Update flags for close distances
            if (start_close <= min_distance_start) and (junc_start > start_limit):
                any_start_close = True
                min_distance_start = start_close
                closest_start = junc_start
            if (end_close <= min_distance_end) and (junc_end < end_limit):
                any_end_close = True
                min_distance_end = end_close
                closest_end = junc_end

            # # Update flag for any close distance
            # if (start_close <= min_distance) or (end_close <= min_distance):
            #     any_close = True

            # Select best junction based on criteria
            if (start_close <= min_distance) and (end_close <= min_distance) and (junc_start > start_limit) and (junc_end < end_limit):
                if (current_distance < best_distance):# and (gap_percent > best_percent):
                    best_junc = junc
                    best_distance = current_distance
                    #best_percent = gap_percent
            
    # Handle cases where no suitable junction is found
    if best_junc is None:
        if (((~any_start_close) or (~any_end_close)) and (gap_size < 200)):
            return None
        else:
            
            best_junc = 'gap:'
            if any_start_close:
                best_junc += str(closest_start)
            best_junc += ':'
            if any_end_close:
                best_junc += str(closest_end)
            if best_junc == 'gap::':
                best_junc = None

    else:
        min_overall_distance = min_distance_start + min_distance_end
        if (best_distance > 100) and (min_overall_distance < (best_distance - 50)) and (any_start_close and any_end_close):
            best_junc = 'gap:'
            if any_start_close:
                best_junc += str(closest_start)
            best_junc += ':'
            if any_end_close:
                best_junc += str(closest_end)

            if best_junc == 'gap::':
                best_junc = None


    # if best_junc == 'gap::':
    #     print(gap_min_size)
    return best_junc


def correct_factor(factor, junctions_bed, coords, coords_all, correct_bias=True, binary_fraction=0.25, strand=None, print_gene = None):
    """
    Function: correct_factor
    
    Description:
    This function corrects biases in a given factor representing genomic features. It employs 
    a series of processing steps including bias correction, smoothing, segmentation, and 
    exon correction based on junction information. The corrected genomic segments are returned 
    as a DataFrame.
    
    Input:
    - factor: An array-like object representing genomic features.
    - junctions_bed: A BedTool object containing genomic junction information.
    - coords: A list of genomic coordinates indicating the start and end of each segment.
    - correct_bias: A boolean indicating whether to perform bias correction (default: True).
    - smooth_fraction: A float specifying the fraction of smoothing for the factor (default: 0.25).
    - strand: A string indicating the strand information for bias correction (default: None).
    
    Output:
    - corrected_exons: A pandas DataFrame containing the corrected genomic segments with 
                       columns 'chrom', 'start', and 'end'.
    
    Dependencies:
    - numpy: Numerical computing library for Python.
    - pandas: Data manipulation library for Python.
    """

    coords_gap = int(coords[2].split(':')[1]) - int(coords[1].split(':')[1])

    factor, coords = correct_factor_coords(factor, coords, coords_all)

    
    # Perform bias correction if specified
    if correct_bias:
        y = factor_lm(factor, strand)
    else:
        y = factor
        
    # Normalize the factor
    y = y / np.quantile(y, 0.99)
    y = np.minimum(y, 1)
    
    # Smooth the factor
    binary_y = binarize_factor(y, binary_fraction=binary_fraction, exon_quant=1)
    
    # Identify continuous segments
    segments_df = find_continuous_segments(coords, binary_y)

    ##### CORRECTED FOR PKM 7/15/2024

    segment_min = np.min([coords_gap*5, 15])

    print(segment_min)

    segments_df = segments_df.loc[(segments_df.end - segments_df.start) >= segment_min].reset_index()
    #####
    
    # Check the number of segments
    n_segments = segments_df.shape[0]
    if (n_segments < 2) or (len(junctions_bed.to_dataframe()) == 0):
        return segments_df
    
    # Extract gaps between segments
    factor_junctions_bed = get_factor_gaps(segments_df)

    # if (len(factor_junctions_bed) == 0) or (len(junctions_bed.to_dataframe()) == 0):
    #     return segments_df
    
    # Intersect factor junctions with junctions_bed
    bed_intersection = BedTool.from_dataframe(factor_junctions_bed).intersect(junctions_bed, wao=True).to_dataframe()
    bed_intersection.columns = ['chrom', 'start', 'end', 'names', 'junc_chrom', 'junc_start', 'junc_end', 'junc_names', 'overlap']
    
    # Convert junctions_bed to DataFrame
    junctions_bed_df = junctions_bed.to_dataframe()
    junctions_bed_df.columns = ['chrom', 'start', 'end', 'junc_names']
    
    # Correct exons based on junction information
    corrected_exons = get_corrected_exons(segments_df, factor_junctions_bed, bed_intersection, junctions_bed_df, coords)

    corrected_exons = corrected_exons.reset_index(drop=True)
    
    return corrected_exons




def get_isoforms(EF, junctions_bed, coordinates, coordinates_all, correct_bias=True, binary_fraction=0.25, strand=None, print_gene=None):
    """
    Function: get_isoforms
    
    Description:
    This function generates isoforms based on given expression factors (EF) for different genes.
    It iterates through each expression factor column, corrects biases, and smoothes the factors 
    to generate isoforms. The corrected isoforms are returned as a dictionary.
    
    Input:
    - EF: A 2D numpy array containing expression factors (EF) for different genes.
    - junctions_bed: A BedTool object containing genomic junction information.
    - coordinates: A list of genomic coordinates indicating the start and end of each segment.
    - correct_bias: A boolean indicating whether to perform bias correction (default: True).
    - smooth_fraction: A float specifying the fraction of smoothing for the factors (default: 0.25).
    - strand: A string indicating the strand information for bias correction (default: None).
    
    Output:
    - isoforms: A dictionary where keys represent isoform names and values represent 
                corresponding corrected isoforms (pandas DataFrames).
    
    Dependencies:
    - numpy: Numerical computing library for Python.
    - correct_factor: A function to correct expression factors.
    """
    # Extract the number of isoforms
    K = EF.shape[1]

    # Exclude first and last rows of EF
    EF = EF[1:-1]

    # Check dimensions
    assert EF.shape[0] == len(coordinates), "Dimensions mismatch"

    # Initialize dictionary to store isoforms
    isoforms = {}
    
    # Iterate through each expression factor
    for k in range(K):
        print(k+1)
        # Extract factor for current isoform
        factor = np.array(EF[:, k])
        
        # Correct the factor to generate isoform
        isoform_k = correct_factor(factor, junctions_bed, coordinates, coordinates_all, correct_bias=correct_bias, binary_fraction=binary_fraction, strand=strand, print_gene = print_gene)
        
        # Store isoform in dictionary
        isoforms[f'factor_{str(k+1)}'] = isoform_k

    return isoforms

def correct_factor_coords(factor, coords_factor, coords_all):
    coords1 = int(coords_factor[0].split(':')[1])
    coords2 = int(coords_factor[1].split(':')[1])

    coords_diff = coords2 - coords1
    if coords_diff == 1:
        assert len(factor) == len(coords_all)
        return factor, coords_all
    else:
        new_factor = []
        for f in factor:
            new_factor.extend(([f]*coords_diff))
        if len(new_factor) > len(coords_all):
            new_factor = new_factor[:len(coords_all)]
        elif len(coords_all) > len(new_factor):
            coords_all = coords_all[:len(new_factor)]
        return new_factor, coords_all

def merge_isoforms(isoforms):
    """
    Function: merge_isoforms
    
    Description:
    This function merges isoforms that have similar exon structures based on certain criteria. 
    It iterates through each isoform, compares its exon structure with existing isoforms, 
    and merges isoforms if they meet specific conditions.
    
    Input:
    - isoforms: A dictionary containing isoforms where keys represent isoform names and 
                values represent corresponding corrected isoforms (pandas DataFrames).
    
    Output:
    - isoform_list: A dictionary containing merged isoforms where keys represent merged 
                    isoform names and values contain information about merged isoforms, 
                    including the merged DataFrame and the list of factors contributing 
                    to each merged isoform.
    """
    # Initialize dictionary to store merged isoforms
    isoform_list = {'isoform_1': {'df': isoforms['factor_1'], 'factors': []}}
    i = 2
    
    # Iterate through each isoform
    for name_iso in isoforms:
        iso = isoforms[name_iso]
        
        found_match = False
        
        # Compare with existing merged isoforms
        for name_iso2 in isoform_list:
            iso2 = isoform_list[name_iso2]['df']
            
            # Skip if already found a match
            if found_match:
                continue
            else:
                # Compare exon structures
                if iso.shape[0] == iso2.shape[0]:
                    length_iso = iso.shape[0]
                    exon_len_iso = iso['end'] - iso['start']
                    exon_len_iso2 = iso2['end'] - iso2['start']
                    
                    start_diff = abs(iso['start'] - iso2['start'])
                    end_diff = abs(iso['end'] - iso2['end'])
                    
                    # Check if exon structures are similar
                    if length_iso == 1:
                        suma_start = 0
                        suma_end = 0
                    else:
                        suma_start = sum(start_diff.iloc[1:])
                        suma_end = sum(end_diff.iloc[:-1])
                    
                    first_exon_diff = start_diff.iloc[0]
                    last_exon_diff = end_diff.iloc[-1]
                    
                    frac_start_1 = first_exon_diff / exon_len_iso.iloc[0]
                    frac_start_2 = first_exon_diff / exon_len_iso2.iloc[0]
        
                    frac_end_1 = last_exon_diff / exon_len_iso.iloc[-1]
                    frac_end_2 = last_exon_diff / exon_len_iso2.iloc[-1]
                    
                    first_exon_is_diff = (first_exon_diff > 200) and ((frac_start_1 > 0.25) or (frac_start_2 > 0.25))
                    last_exon_is_diff = (last_exon_diff > 200) and ((frac_end_1 > 0.25) or (frac_end_2 > 0.25))
                    
                    # Merge isoforms if exon structures are similar
                    if (suma_start == 0) and (suma_end == 0) and (not first_exon_is_diff) and (not last_exon_is_diff):
                        isoform_list[name_iso2]['factors'].append(name_iso)
                        found_match = True
                        
        # If no match found, create new isoform entry
        if not found_match:
            new_name_iso = f'isoform_{i}'
            isoform_list[new_name_iso] = {'df': iso, 'factors': [name_iso]}
            i += 1
            
    return isoform_list

def load_torino_rds(rds_file):
    readRDS = ro.r['readRDS']
    df = readRDS(rds_file)
    with (ro.default_converter + pandas2ri.converter).context():
        pd_from_r_df = ro.conversion.get_conversion().rpy2py(df)

    return pd_from_r_df