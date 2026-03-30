import numpy as np
import pandas as pd
import os

from matplotlib import pyplot as plt
import seaborn as sns
import tabix
from pybedtools import BedTool
import sys

from rpy2.rinterface_lib.sexp import NULLType

print('custom-made modules')

# sys.path.append('/project2/mstephens/cfbuenabadn/gtex-stm/code/scripts')
#from get_isoforms import *
#from prepare_counts import *


import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

def run_tabix_on_junc(junc_file, coords, transform_to_df = True):
    chrom, start = coords[0].split(':')
    start = int(start)
    end = int(coords[-1].split(':')[1])
    
    tb = tabix.open(junc_file)
    juncs = tb.query(chrom, start, end)

    if transform_to_df:

        junctions_bed = []
        junctions_bed_cols = []
        
        # Iterate over queried junctions and store them in the DataFrame
        for idx, record in enumerate(juncs):
            junctions_bed.append(pd.Series(record))
            junctions_bed_cols.append(f'record{str(idx)}')

        if len(junctions_bed) == 0:
            return pd.DataFrame(columns = ['chrom', 'start', 'end', 'gene'])
        junctions_bed = pd.concat(junctions_bed, axis=1)
        junctions_bed.columns = junctions_bed_cols
        junctions_bed = junctions_bed.T
        return junctions_bed

    else:
        return juncs


def binarize_factor(factor, cutoff_=0.25, binary_fraction=0.5, step_fraction=10, exon_quant=0.99, cutoff_strict=0.1, pos_window=100, min_gap = 200):
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
    last_pos = 0

    # Iterate through factor values
    # for idx in range(len(factor)):
    for idx, f in enumerate(factor):
        # f = factor[idx]

        # Apply strict cutoff
        if f <= cutoff_strict:
            factor_binary[idx] = 0
            # last_pos = None
            continue

        # Determine binary value based on position and previous values
        if idx == 0:
            if f > cutoff:
                f_binary = 1
                # last_pos = 0
            else:
                f_binary = 0
                # last_pos = 0
        else:
            if (factor_binary[idx - 1] == 0) and ((idx - last_pos) > min_gap):
                if f < cutoff:
                    f_binary = 0
                    # last_pos = None
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
                    # last_pos = None
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


def load_ebpmf_gene(gene, species):
    readRDS = ro.r['readRDS']
    if species in ['Mouse', 'Rat', 'Rabbit', 'Macaque', 'Chicken']:
        df = readRDS(f'results/torino/factorization/{species}/{gene}.rds')
    else:
        df = readRDS(f'results/torino/factorization_filtered/{species}/{gene}.rds')
    with (ro.default_converter + pandas2ri.converter).context():
        pd_from_r_df = ro.conversion.get_conversion().rpy2py(df)

    #coords = list(pd_from_r_df['torino_K10']['coords'])[1:-1]#pd_from_r_df['torino_K10']['coords']

    #junctions_bed = run_tabix_on_junc(junc_file, coords)

    output = {'rds':pd_from_r_df#, 
              #'junctions_bed':junctions_bed
             }

    return output


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

        # print(best_junction)

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

    return best_junc



def get_unmerged_isoforms(gene, junc_file, species, strand, gene_id = '.', use_low_juncs = False):
    rds = load_ebpmf_gene(gene, species)

    if type(rds['rds']['torino_K10']) is NULLType:
        return None
        
    #junc_file = '/project2/mstephens/cfbuenabadn/gtex-stm/code/junctions.tab.gz'
    coords = list(rds['rds']['torino_K10']['coords'])[1:-1]
    coords_all = list(rds['rds']['coords_all'])

    junctions_bed = get_junctions_bed(coords, junc_file, gene_id, species, use_low_juncs=use_low_juncs)

    isoforms = get_isoforms_junc_driven(rds['rds']['torino_K10']['EF_smooth'], junctions_bed, coords, coords_all, strand, gene)
    return isoforms

def get_isoforms_junc_driven(EF, junctions_bed, coordinates, coordinates_all, strand, gene):
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
    if EF.shape[0] == len(coordinates):

        # Initialize dictionary to store isoforms
        isoforms = {}
        
        # Iterate through each expression factor
        for k in range(K):
            # print(k+1)
            # Extract factor for current isoform
            factor = np.array(EF[:, k])
            
            # Correct the factor to generate isoform
            isoform_k = get_factor_isoform(factor, junctions_bed, coordinates, coordinates_all, strand)
            isoform_k['gene_id'] = [gene]*len(isoform_k)
            isoform_k['strand'] = [strand]*len(isoform_k)
    
            
    
            #isoform_k = isoform_k[['chrom', 'start', 'end', 'gene_id', 'transcript_id', 'strand', 'factors', 'exon_id']]
            
            # Store isoform in dictionary
            isoforms[f'factor_{str(k+1)}'] = isoform_k

            if len(isoform_k) == 0:
                print('Empty isoforms')
                return None

    else:
        print("Dimensions mismatch")
        isoforms = None

    return isoforms

def correct_intron_retention(corrected_exons, y, coords, juncs):
    chrom = coords[0].split(':')[0]
    coord_start = int(coords[0].split(':')[1])

    correct_ir = False
    found_correction = False

    new_exons = []
    
    for idx, row in corrected_exons.iterrows():
        start_idx = coords.index(chrom + ':' + str(row.start))
        end_idx = coords.index(chrom + ':' + str(row.end))

        if (np.quantile(y[start_idx:end_idx], 0.9) < 0.1) and ((row.end - row.start) > 200):

            correct_ir = True

            
            juncs_in_exon = juncs.loc[(juncs.start >= row.start) & (juncs.end <= row.end)]


            if len(juncs_in_exon) > 0:
                first_start = juncs_in_exon.start.min()
                last_end = juncs_in_exon.end.max()
                start_idx_junc = coords.index(chrom + ':' + str(first_start))
                end_idx_junc = coords.index(chrom + ':' + str(last_end))

                e1 = np.quantile(y[start_idx:start_idx_junc], 0.9)
                e2 = np.quantile(y[end_idx_junc:end_idx], 0.9)
                i = np.quantile(y[start_idx_junc:end_idx_junc], 0.9)


                if (e1 > 0.1) and (e2 > 0.1) and (i < 0.1):
                    new_exons.append([row.chrom, row.start, first_start-1])
                    new_exons.append([row.chrom, last_end+1, row.end])
                    found_correction = True
        else:
            new_exons.append([row.chrom, row.start, row.end])

    new_exons_df = pd.DataFrame(new_exons)
    new_exons_df.columns = ['chrom', 'start', 'end']

    return new_exons_df, correct_ir, found_correction
    
def correct_factor_coords(factor, coords_factor, coords_all):


    diffs = [int(coords_factor[i+1].split(':')[1]) - int(coords_factor[i].split(':')[1]) for i in range(len(coords_factor)-1)]

    coords_diff = int(np.median(diffs))
    
    # coords1 = int(coords_factor[0].split(':')[1])
    # coords2 = int(coords_factor[1].split(':')[1])

    # coords_diff = coords2 - coords1
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
                
                    


def get_factor_isoform(factor, junctions_bed, coords, coords_all, strand, binary_fraction=0.5, cutoff_strict = 0.1):
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
    # if correct_bias:
    #     y = factor_lm(factor, strand)
    # else:
    y = factor
        
    # Normalize the factor
    y = y / np.quantile(y, 0.99)
    y = np.minimum(y, 1)

    junctions_bed_df = junctions_bed.to_dataframe()
    if len(junctions_bed_df) == 0:
        junctions_bed_df = pd.DataFrame(columns = ['chrom', 'start', 'end', 'junc_names', 'junc_counts'])
        
    junctions_bed_df.columns = ['chrom', 'start', 'end', 'junc_names', 'junc_counts']


    factor_notes = ''
    
    compat_juncs = []

    for idx, junc in junctions_bed_df.iterrows():
        if get_junction_compatibility_score(y, coords, junc):
            compat_juncs.append(junc.junc_names)

    # print(compat_juncs)

    if len(compat_juncs) > 0:
        discard_junc_approach = False
        best_juncs = get_best_compatible_loop(junctions_bed_df.loc[junctions_bed_df.junc_names.isin(compat_juncs)], y, coords, k = 300)
        # best_juncs = get_best_compatible(junctions_bed_df.loc[junctions_bed_df.junc_names.isin(compat_juncs)], y, coords, k = 300)

        # print(best_juncs)
        # print('')

        if len(best_juncs) == 0:
            corrected_exons = get_exons(y, coords, junctions_bed_df.loc[junctions_bed_df.junc_names.isin(compat_juncs)], junctions_bed_df, strand)

        else:
            corrected_exons = get_exons(y, coords, junctions_bed_df.loc[junctions_bed_df.junc_names.isin(best_juncs)], junctions_bed_df, strand)
            


        #### See if most weight is outside of inferred transcript
        idx = int(coords[0].split(':')[1])
        transcript_start = corrected_exons.start.min() - idx
        transcript_end = corrected_exons.end.max() - idx

        if (np.sum(y[transcript_start:transcript_end])/np.sum(y)) <= 0.667:
            print('Junc approach misses 1/3 or more of transcript.')
            discard_junc_approach = True

        # Scan for fake IR from lack of good junctions and 3' bias

        for idx, row in corrected_exons.iterrows():
            chrom = row.chrom
            start = row.start
            end = row.end
            coord_start = chrom + ':' + str(start)
            coord_end = chrom + ':' + str(end)
            try:
                start_idx = coords.index(coord_start)
                end_idx = coords.index(coord_end)
                exon_body = y[start_idx:end_idx]/np.quantile(y[start_idx:end_idx], 0.99)
            except:
                coord_array = np.array([int(x.split(':')[1]) for x in coords])
                coord_start_array = (coord_array >= int(start))
                coord_end_array = (coord_array <= int(end))
                exon_body = np.array(y)[(coord_start_array & coord_end_array)]
            

            if np.median(exon_body) < 0.1:
                print('Spurious retained intron')
                discard_junc_approach = True
        
        
        

        #Now we correct for spurious intron retention (caused by overlapping transcripts or other issues)

        # Smooth the factor
        binary_y = binarize_factor(y, binary_fraction=binary_fraction, exon_quant=1, cutoff_strict = cutoff_strict)
        
        # Identify continuous segments
        segments_df = find_continuous_segments(coords, binary_y)
        
        segment_min = np.min([coords_gap*5, 15])
    
        # print(segment_min)
    
        if len(segments_df) == 0:
            print('odd segment issue')
            segment_df = pd.DataFrame(columns = ['chrom', 'start', 'end'])

        try:

            segments_df = segments_df.loc[(segments_df.end - segments_df.start) >= segment_min].reset_index(drop=True)

        except:
            print('Strange segments_df error')
            print(segments_df)
            segments_df = pd.DataFrame(columns = ['chrom', 'start', 'end'])
        #####
        
        # Check the number of segments
        n_segments = segments_df.shape[0]
        if n_segments >= 2:

            corrected_exons_bed = BedTool.from_dataframe(corrected_exons)

            ### Check for spurious intron retention events (this should be rare)
            
            factor_junctions_bed = get_factor_gaps(segments_df)

            factor_junctions_bed = factor_junctions_bed.loc[(factor_junctions_bed.end - factor_junctions_bed.start) > 200]

            if len(factor_junctions_bed) > 0:

                try:

                    new_exons_df, correct_ir, found_correction = correct_intron_retention(corrected_exons, y, coords, junctions_bed_df)
                except:
                    print('Missing junc coordinates.')
                    new_exons_df = None
                    correct_ir = False
                    found_correction = False

                if correct_ir and found_correction:
                    print('Unannotated exon skipping? This should be rare.')
                    factor_notes += 'Unannotated exon-skipping,'
                    corrected_exons = new_exons_df

                elif correct_ir and (not found_correction):
                    print("Couldn't correct intron retention. This should be rare.")
                    #discard_junc_approach = True

                else:

                    factor_junctions_bed = BedTool.from_dataframe(factor_junctions_bed)
                    
                    bad_overlaps = factor_junctions_bed.intersect(corrected_exons_bed, f=0.8).to_dataframe()
    
                    # if len(bad_overlaps) > 0:
                    #     discard_junc_approach = True
                    #     print('This should be rare')
    
                    ### Check for intronic transcripts
                    # print(segments_df)
                    segments_bed = BedTool.from_dataframe(segments_df).merge(d = 500)
    
                    junction_introns = BedTool.from_dataframe(get_factor_gaps(corrected_exons))
    
    
    
                    
                    outstanding_exons = segments_bed.intersect(corrected_exons_bed, v=True).intersect(junction_introns, f=1).to_dataframe()
                    if len(outstanding_exons) > 0:
                        for idx, row in outstanding_exons.iterrows():
                            chrom = row.chrom
                            start = row.start
                            end = row.end
                        
                            coord_start = chrom + ':' + str(start)
                            coord_end = chrom + ':' + str(end)
                        
                            idx_start = coords.index(coord_start)
                            idx_end = coords.index(coord_end)
        
                            outstanding_y = y[idx_start:idx_end]
                            median_outstanding_y = np.median(outstanding_y)
                            if median_outstanding_y >= 0.5:
                                outstanding_df = pd.DataFrame([[chrom, start, end]], columns=['chrom', 'start', 'end'])
        
                                corrected_exons = pd.concat([corrected_exons, outstanding_df], ignore_index=True).reset_index(drop=True).sort_values('start')
                                factor_notes += 'Intronic transcript,'
                                print('Intronic transcript found. This should be rare.')

    else:
        discard_junc_approach = True


    

    if discard_junc_approach:

        print("Let's try again")

        factor_notes += 'Old approach,'
    
        # Smooth the factor
        binary_y = binarize_factor(y, binary_fraction=binary_fraction, exon_quant=1, cutoff_strict = cutoff_strict)
        
        # Identify continuous segments
        segments_df = find_continuous_segments(coords, binary_y)
    
        ##### CORRECTED FOR PKM 7/15/2024
    
        segment_min = np.min([coords_gap*5, 15])
    
        # print(segment_min)
    
        if len(segments_df) == 0:
            print('odd segment issue')
            segment_df = pd.DataFrame(columns = ['chrom', 'start', 'end'])

        try:

            segments_df = segments_df.loc[(segments_df.end - segments_df.start) >= segment_min].reset_index(drop=True)

        except:
            print('Strange segments_df error')
            print(segments_df)
            segments_df = pd.DataFrame(columns = ['chrom', 'start', 'end'])
        #####
        
        # Check the number of segments
        n_segments = segments_df.shape[0]
        if (n_segments < 2) or (len(junctions_bed.to_dataframe()) == 0):
            corrected_exons = segments_df #return segments_df

        else:
        
            # Extract gaps between segments
            factor_junctions_bed = get_factor_gaps(segments_df)
        
            # if (len(factor_junctions_bed) == 0) or (len(junctions_bed.to_dataframe()) == 0):
            #     return segments_df
            
            # Intersect factor junctions with junctions_bed
            bed_intersection = BedTool.from_dataframe(factor_junctions_bed).intersect(junctions_bed, wao=True).to_dataframe()
            bed_intersection.columns = ['chrom', 'start', 'end', 'names', 'junc_chrom', 'junc_start', 'junc_end', 'junc_names', 'junc_counts', 'overlap']
            
            # Convert junctions_bed to DataFrame
            junctions_bed_df = junctions_bed.to_dataframe()
            junctions_bed_df.columns = ['chrom', 'start', 'end', 'junc_names', 'junc_counts']
            
            # Correct exons based on junction information
            corrected_exons = get_corrected_exons(segments_df, factor_junctions_bed, bed_intersection, junctions_bed_df, coords)
        
            corrected_exons = corrected_exons.reset_index(drop=True)


    corrected_exons['notes'] = [factor_notes] * len(corrected_exons)
    corrected_exons = corrected_exons.loc[(corrected_exons.start >= 1) & (corrected_exons.end >= 1) & (corrected_exons.end > corrected_exons.start)]

    # print(corrected_exons.columns)
    return corrected_exons



def get_splice_site_compatibility_score(y, splice_site, donor, k=30, min_thr = 0.05, min_ratio = 1):
    start = np.max([0, splice_site-k])
    end = np.min([splice_site+k, len(y)])
    if donor:
        exon = y[start:splice_site]
        intron = y[splice_site:end]
    else:
        exon = y[splice_site:end]
        intron = y[start:splice_site]

    if (len(exon) == 0) or (len(intron) == 0):
        return 0, 0, exon, intron

    # print(exon)
    exon_cov = np.quantile(exon, 0.95)
    intron_cov = np.quantile(intron, 0.95)

    if exon_cov < min_thr:
        compatibility_score = 0
        # print(exon_cov)

    else:

        total_diff = exon_cov - intron_cov
        compatibility_score = ((exon_cov/intron_cov)-1)*total_diff
    
    return compatibility_score, exon_cov, exon, intron

def get_coverage_score(y, donor, acceptor, donor_exon, acceptor_exon, q=0.9, intron_thr = 0.5):
    intron_q = np.quantile(y[donor:acceptor], q)
    donor_q = np.quantile(donor_exon, q)
    acceptor_q = np.quantile(acceptor_exon, q)

    exon_q = np.min([donor_q, acceptor_q])

    if intron_q > intron_thr:
        return 0

    return 0.75 - (intron_q/exon_q)
    

def get_junction_compatibility_score(y, coords, junc):
    donor_coord = junc.chrom + ':' + str(junc.start)
    acceptor_coord = junc.chrom + ':' + str(junc.end)

    if (donor_coord not in coords) or (acceptor_coord not in coords):
        return 0

    donor_idx = coords.index(donor_coord)
    acceptor_idx = coords.index(acceptor_coord)

    intron_len = acceptor_idx - donor_idx

    k = int(np.min([200, intron_len/5]))

    # print(donor_idx)
    # print(k)

    donor_score, donor_exon, e, i = get_splice_site_compatibility_score(y, donor_idx, True, k=k)
    acceptor_score, acceptor_exon, e, i = get_splice_site_compatibility_score(y, acceptor_idx, False, k=k)

    if (acceptor_idx - donor_idx) <= 10:
        print('Very short junction')
        print(acceptor_idx - donor_idx)
        coverage_score = 0

    else:
    
        coverage_score = get_coverage_score(y, donor_idx, acceptor_idx, donor_exon, acceptor_exon)

    all_compatible = int(np.all(
        [
        donor_score > 0,
        acceptor_score > 0,
        coverage_score > 0
    ]
    ))

    # print(donor_score, acceptor_score, coverage_score)

    return all_compatible


def get_best_compatible(compat_juncs, y, coords, k = 300):
    best_junc_list = []
    discarded_juncs = []
    discarded_dict = {}
    for junc in list(compat_juncs.junc_names):
        if (junc in discarded_juncs) or (junc in best_junc_list):
            continue
        junc_bed = BedTool.from_dataframe(pd.DataFrame(compat_juncs.loc[compat_juncs.junc_names == junc]))
        other_junc_bed = BedTool.from_dataframe(pd.DataFrame(compat_juncs.loc[compat_juncs.junc_names != junc]))

        overlap_bed =  junc_bed.intersect(other_junc_bed, wo=True).to_dataframe()
        
        if len(overlap_bed) == 0:
            best_junc_list.append(junc)

        else:
            overlap_bed.columns = ['chrom', 'start', 'end', 'junc', 'junc_counts', 'chr_', 'start_', 'end_', 'alt_junc', 'overlap']

            list_of_overlaps = [junc] + list(overlap_bed.alt_junc)

            start = np.min([overlap_bed.start.min(), overlap_bed.start_.min()])
            end = np.max([overlap_bed.end.max(), overlap_bed.end_.max()])



            start_coord = compat_juncs.chrom.iloc[0] + ':' + str(start)
            end_coord = compat_juncs.chrom.iloc[0] + ':' + str(end)
        
            start_idx = coords.index(start_coord)
            end_idx = coords.index(end_coord)

            # print(overlap_bed)

            # print(start_idx)
            # print(end_idx)

            start_idx = np.max([0, start_idx - 300])
            end_idx = np.min([end_idx + 300, len(y)-1])

            coords_slice = coords[start_idx:end_idx]
            y_slice = y[start_idx:end_idx]
            binary_y = binarize_factor(y_slice, binary_fraction=0.5, exon_quant=1, cutoff_strict = 0.1)


            segments_df = find_continuous_segments(coords_slice, binary_y)

            coords_gap = int(coords[2].split(':')[1]) - int(coords[1].split(':')[1])

            segment_min = np.min([coords_gap*5, 15])


            if len(segments_df) == 0:
                print('odd segment issue')
                segment_df = pd.DataFrame(columns = ['chrom', 'start', 'end'])
    
            try:

                segments_df = segments_df.loc[(segments_df.end - segments_df.start) >= segment_min].reset_index(drop=True)
    
            except:
                print('Strange segments_df error')
                print(segments_df)
                segments_df = pd.DataFrame(columns = ['chrom', 'start', 'end'])
    
            factor_junctions_bed = get_factor_gaps(segments_df)

            # plt.plot(y_slice)
            # plt.plot(binary_y)
            # plt.show()

            current_score = 1e10
            best_junc = ''

            # print(list_of_overlaps)

            for junc_test in list_of_overlaps:
                start_junc = int(compat_juncs.loc[compat_juncs.junc_names == junc_test].start.iloc[0])
                end_junc = int(compat_juncs.loc[compat_juncs.junc_names == junc_test].end.iloc[0])

                for idx, fjunc in factor_junctions_bed.iterrows():
                    junc_score = np.abs(fjunc.start - start_junc) + np.abs(fjunc.end - end_junc)
                    if junc_score < current_score:
                        best_junc = junc_test
                        current_score = junc_score

            discarded_juncs.extend([x for x in list_of_overlaps if x != best_junc])
            discarded_dict.update({best_junc:[x for x in list_of_overlaps if x != best_junc]})
            best_junc_list.append(best_junc)

    return best_junc_list


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





def get_exons(y, coords, juncs, all_juncs, strand):

    binary_y = binarize_factor(y, binary_fraction=0.5, exon_quant=1, cutoff_strict = 0.1)
    
    juncs = juncs.reset_index(drop=True)
    chrom = juncs.chrom.iloc[0]
    start = juncs.start.min()
    end = juncs.end.max()

    coord_start = chrom + ':' + str(start)
    coord_end = chrom + ':' + str(end)

    idx_start = coords.index(coord_start)
    idx_end = coords.index(coord_end)

    # print(idx_end)
    # print(len(y))

    if idx_start <= 30:
        exon_1 = int(coords[0].split(':')[1])

    else:
        y_slice = np.array(y[0:idx_start+30])
        y_slice = y_slice/np.max(y_slice)
        coords_slice = coords[0:idx_end+30]
        if (strand == '+') and (len(y_slice) < len(y)/5):
            
            binary_start = binarize_factor(y_slice, binary_fraction=0.25, exon_quant=1, cutoff_strict = 0.1)
        else:
            binary_start = binarize_factor(y_slice, binary_fraction=0.5, exon_quant=1, cutoff_strict = 0.25)

        segments_df = find_continuous_segments(coords_slice, binary_start)

        if len(segments_df) > 1:
            segments_df_other = segments_df.iloc[:-1]
            d = np.max([(segments_df_other.end - segments_df_other.start).max(), 50])

        else:

            d = int(len(y_slice)/10)

        
        
        segments_df = BedTool.from_dataframe(segments_df).merge(d=d).to_dataframe()

        # plt.plot(y_slice)
        # plt.plot(binary_start)
        # plt.show()



        segment_min = 15#np.min([coords_gap*5, 30])

        if len(segments_df) == 0:
            print('odd segment issue')
            segment_df = pd.DataFrame(columns = ['chrom', 'start', 'end'])

        try:

            segments_df = segments_df.loc[(segments_df.end - segments_df.start) >= segment_min].reset_index(drop=True)

        except:
            print('Strange segments_df error')
            print(segments_df)
            segments_df = pd.DataFrame(columns = ['chrom', 'start', 'end'])

        exon_1 = int(start)-31
        if len(segments_df) >= 1:
            # exon_1 = int(start)-31
        # else:
            # if strand == '-':
            current_median = -1e100
            coords_start = int(coords_slice[0].split(':')[1])
            for idx, row in segments_df[::-1].iterrows():
                row_end_idx =  (row.end) - coords_start
                row_start_idx =  (row.start) - coords_start
                y_segment = y_slice[row_start_idx:row_end_idx]
                peak_median = np.median(y_segment)
                if peak_median > current_median:
                    exon_1 = int(row.start)
                    current_median = peak_median
            # else:
            #     exon_1 = segments_df.start.max()

        junc_ends = (all_juncs.end - exon_1).abs().sort_values()
        junc_ends = all_juncs.loc[junc_ends < 300]

        if len(junc_ends) > 0:
            # junc_ends.sort_values('junc_counts').index[-1]
            junc_end = junc_ends.sort_values('junc_counts').index[-1]
            exon_1_alt = all_juncs.loc[junc_end].end

            if exon_1_alt < (int(juncs.iloc[0].start)-15):
                exon_1 = exon_1_alt

        # junc_ends = (all_juncs.end - exon_1).abs().sort_values()

        # if junc_ends.min() < 100:
        #     junc_end = junc_ends.index[0]
        #     exon_1 = all_juncs.loc[junc_end].end


    if (len(y)-idx_end) <= 30:
        exon_last = int(coords[-1].split(':')[1])
    else:
        y_slice = np.array(y[idx_end-30:])
        y_slice = y_slice/np.max(y_slice)
        coords_slice = coords[idx_end-30:]
        if (strand == '-') and (len(y_slice) < len(y)/5):
            binary_end = binarize_factor(y_slice, binary_fraction=0.25, exon_quant=1, cutoff_strict = 0.1)
        else:
            
            binary_end = binarize_factor(y_slice, binary_fraction=0.5, exon_quant=1, cutoff_strict = 0.25)

        # plt.plot(y_slice)
        # plt.plot(binary_end)
        # plt.show()

        segments_df = find_continuous_segments(coords_slice, binary_end)

        if len(segments_df) > 1:
            segments_df_other = segments_df.iloc[1:]
            d = np.max([(segments_df_other.end - segments_df_other.start).max(), 50])

        else:

            d = int(len(y_slice)/10)

        segments_df = BedTool.from_dataframe(segments_df).merge(d=d).to_dataframe()

        

        
        
        segment_min = 15#np.min([coords_gap*5, 30])

        if len(segments_df) == 0:
            print('odd segment issue')
            segment_df = pd.DataFrame(columns = ['chrom', 'start', 'end'])

        try:

            segments_df = segments_df.loc[(segments_df.end - segments_df.start) >= segment_min].reset_index(drop=True)

        except:
            print('Strange segments_df error')
            print(segments_df)
            segments_df = pd.DataFrame(columns = ['chrom', 'start', 'end'])

        exon_last = int(end)+31
        if len(segments_df) >= 1:
            # exon_last = int(end)+31
        # else:
            # if strand == '+':
            current_median = -1e100
            coords_start = int(coords_slice[0].split(':')[1])
            for idx, row in segments_df.iterrows():
                row_end_idx =  (row.end) - coords_start
                row_start_idx =  (row.start) - coords_start
                y_segment = y_slice[row_start_idx:row_end_idx]
                peak_median = np.median(y_segment)
                if peak_median > current_median:
                    exon_last = int(row.end)
                    current_median = peak_median
            # else:
            #     exon_last = segments_df.end.min()

        


        

        junc_starts = (all_juncs.start - exon_last).abs().sort_values()
        junc_starts = all_juncs.loc[junc_starts < 300]

        if len(junc_starts) > 0:
            # junc_ends.sort_values('junc_counts').index[-1]
            junc_start = junc_starts.sort_values('junc_counts').index[-1]
            exon_last_alt = all_juncs.loc[junc_start].start

            if exon_last_alt > (int(juncs.iloc[-1].end)+15):
                exon_last = exon_last_alt


            
        # junc_starts = (all_juncs.start - exon_last).abs().sort_values()
        
        # if junc_starts.min() < 100:
        #     junc_start = junc_starts.index[0]
        #     exon_last = all_juncs.loc[junc_start].start


    exon_chain = []
    # print(juncs.iloc[:-1])



    
    
    exon_chain.append([chrom, int(exon_1), int(juncs.iloc[0].start)])
    for idx, junc in juncs.iloc[:-1].iterrows():
            

        exon_chain.append([chrom, int(junc.end), int(juncs.iloc[idx+1].start)])

    exon_chain.append([chrom, int(juncs.iloc[-1].end), int(exon_last)])

    exon_df = pd.DataFrame(exon_chain)
    exon_df.columns = ['chrom', 'start', 'end']

    exon_df.start += 1
    exon_df.end -= 1
    corrected_exons_len = exon_df.end - exon_df.start
    exon_df = exon_df.loc[corrected_exons_len >= 3]

    return exon_df#y_slice, binary_end


def process_isoform(isoform_dict, isoform, gene_id='.', strand = '.'):
    df = isoform_dict['df']
    factors = isoform_dict['factors']
    n = len(df)
    
    df['gene_id'] = [gene_id]*n
    
    if gene_id == '.':
        transcript_id = isoform
    else:
        transcript_id = f'{gene_id}.{isoform}'
        
    df['transcript_id'] = [transcript_id]*n
    df['strand'] = [strand]*n
    df['factors'] = [','.join(factors)]*n
    
    exon_list = [f'exon_{str(j)}' for j in range(1, n+1)]
    df['exon_id'] = exon_list
    
    return df


def isoforms_to_bed(isoforms, gene_id = '.', strand = '.'):

    df = list()
    for isoform, isoform_dict in isoforms.items():
        df_iso = process_isoform(isoform_dict, isoform, gene_id, strand)
        df.append(df_iso)
        
    df = pd.concat(df, axis=0)
    return df


def get_best_compatible_loop(compat_juncs, y, coords, k = 300):

    current_list = list(compat_juncs.junc_names)

    counter = 0
    while True:
    
        best_juncs, discarded_dict, discarded_juncs = loop_compat_juncs(compat_juncs.loc[compat_juncs.junc_names.isin(current_list)], y, coords, k=300)

        if (sorted(best_juncs) == sorted(current_list)) or (counter > 30):
            break
        
        for djunc in discarded_juncs:
            if (djunc in discarded_dict.keys()) and (djunc in best_juncs):
                best_juncs.pop(best_juncs.index(djunc))
                best_juncs.extend(discarded_dict[djunc])

        current_list = best_juncs

        counter += 1

    return best_juncs



def loop_compat_juncs(compat_juncs, y, coords, k=300):
    best_junc_list = []
    discarded_juncs = []
    discarded_dict = {}
    
    for junc in list(compat_juncs.junc_names):
        # print('new junc')
        # print(junc)
        if (junc in discarded_juncs) or (junc in best_junc_list):
            continue
        junc_bed = BedTool.from_dataframe(pd.DataFrame(compat_juncs.loc[compat_juncs.junc_names == junc]))
        other_junc_bed = BedTool.from_dataframe(pd.DataFrame(compat_juncs.loc[compat_juncs.junc_names != junc]))

        overlap_bed =  junc_bed.intersect(other_junc_bed, wo=True).to_dataframe()
        
        if len(overlap_bed) == 0:
            best_junc_list.append(junc)

        else:
            overlap_bed.columns = ['chrom', 'start', 'end', 'junc', 'junc_counts', 'chr_', 'start_', 'end_', 'alt_junc', 'alt_counts', 'overlap']

            list_of_overlaps = [junc] + list(overlap_bed.alt_junc)

            start = np.min([overlap_bed.start.min(), overlap_bed.start_.min()])
            end = np.max([overlap_bed.end.max(), overlap_bed.end_.max()])



            start_coord = compat_juncs.chrom.iloc[0] + ':' + str(start)
            end_coord = compat_juncs.chrom.iloc[0] + ':' + str(end)
        
            start_idx = coords.index(start_coord)
            end_idx = coords.index(end_coord)

            # print(overlap_bed)

            # print(start_idx)
            # print(end_idx)

            start_idx = np.max([0, start_idx - 300])
            end_idx = np.min([end_idx + 300, len(y)-1])

            coords_slice = coords[start_idx:end_idx]
            y_slice = y[start_idx:end_idx]
            binary_y = binarize_factor(y_slice, binary_fraction=0.25, exon_quant=1, cutoff_strict = 0.25)


            segments_df = find_continuous_segments(coords_slice, binary_y)

            coords_gap = int(coords[2].split(':')[1]) - int(coords[1].split(':')[1])

            segment_min = np.min([coords_gap*5, 15])

            segments_df_ = segments_df.loc[(segments_df.end - segments_df.start) >= segment_min].reset_index(drop=True)

            factor_junctions_bed = get_factor_gaps(segments_df_)
    
            if len(factor_junctions_bed) < 1:
                segments_df_ = segments_df.loc[(segments_df.end - segments_df.start) >= 0].reset_index(drop=True)
    
                factor_junctions_bed = get_factor_gaps(segments_df_)
    
            if len(factor_junctions_bed) < 0:
                print('use score')

                best_junc = list_of_overlaps[0]
                discarded_juncs.extend([x for x in list_of_overlaps if x != best_junc])
                discarded_dict.update({best_junc:[x for x in list_of_overlaps if x != best_junc]})
                best_junc_list.append(best_junc)

            else:
    
                
    
                # print(list_of_overlaps)

                min_juncs = compat_juncs.loc[compat_juncs.junc_names.isin(list_of_overlaps)].junc_counts.max()/10
                refined_juncs = compat_juncs.loc[
                    compat_juncs.junc_names.isin(list_of_overlaps) & (compat_juncs.junc_counts >= min_juncs)
                    ].junc_names

                current_score = 300
                current_counts = 0
                best_junc = ''
    
                for junc_test in refined_juncs:
                    start_junc = int(compat_juncs.loc[compat_juncs.junc_names == junc_test].start.iloc[0])
                    end_junc = int(compat_juncs.loc[compat_juncs.junc_names == junc_test].end.iloc[0])

                    junc_counts = int(compat_juncs.loc[compat_juncs.junc_names == junc_test].junc_counts.iloc[0])
    
                    for idx, fjunc in factor_junctions_bed.iterrows():
                        junc_score = np.abs(fjunc.start - start_junc) + np.abs(fjunc.end - end_junc)
                        # print('junc_score')
                        # print(junc_score)
                        #junc_counts = fjunc.junc_counts
                        if (junc_score < current_score):# and (junc_counts > current_counts/10):
                            best_junc = junc_test
                            current_score = junc_score

                if best_junc == '':            

                    current_score = 1e10
                    current_counts = 0
                    best_junc = ''
        
                    for junc_test in list_of_overlaps:
                        start_junc = int(compat_juncs.loc[compat_juncs.junc_names == junc_test].start.iloc[0])
                        end_junc = int(compat_juncs.loc[compat_juncs.junc_names == junc_test].end.iloc[0])
    
                        junc_counts = int(compat_juncs.loc[compat_juncs.junc_names == junc_test].junc_counts.iloc[0])
        
                        for idx, fjunc in factor_junctions_bed.iterrows():
                            junc_score = np.abs(fjunc.start - start_junc) + np.abs(fjunc.end - end_junc)
                            # print('junc_score')
                            # print(junc_score)
                            #junc_counts = fjunc.junc_counts
                            if (junc_score < current_score):# and (junc_counts > current_counts/10):
                                best_junc = junc_test
                                current_score = junc_score
    
                discarded_juncs.extend([x for x in list_of_overlaps if x != best_junc])
                discarded_dict.update({best_junc:[x for x in list_of_overlaps if x != best_junc]})
                best_junc_list.append(best_junc)

    return best_junc_list, discarded_dict, discarded_juncs



def get_junctions_bed(coords, junc_file, gene_id, species, junc_thres=None, use_low_juncs = False):
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
    junctions_bed = run_tabix_on_junc(junc_file, coords)

    

    
    if junctions_bed.shape[0] == 0:
        junctions_bed = BedTool.from_dataframe(junctions_bed)
        return junctions_bed

    if species == 'Human':
        junctions_bed[3] = [x.split('.')[0] for x in junctions_bed[3]]

    junctions_bed = junctions_bed.loc[junctions_bed[3] == gene_id]
    junctions_bed_counts = junctions_bed.loc[junctions_bed[3] == gene_id, junctions_bed.columns[4:]].astype(int).quantile(0.999,axis=1)

    if junc_thres is None:
        junc_thres = junctions_bed_counts.max()/200

        if use_low_juncs and (junctions_bed_counts.max() < 50):
            junc_thres = -1

    # print(junc_thres)
    
    # Filter junctions for the specified gene ID
    junctions_bed = junctions_bed.loc[:, range(4)]
    junctions_bed.columns = ['chrom', 'start', 'end', 'gene']

    junctions_bed['junc_counts'] = list(junctions_bed_counts)

    

    # print(junctions_bed)


    # print(junctions_bed_counts)

    

    junctions_bed = junctions_bed[junctions_bed_counts >= junc_thres]

    # Convert start and end positions to integers
    junctions_bed.start = junctions_bed.start.astype(int)
    junctions_bed.end = junctions_bed.end.astype(int)

    ##### Make sure juncs are within gene body
    coords_start = int(coords[0].split(':')[1])
    coords_end = int(coords[-1].split(':')[1])
    junctions_bed = junctions_bed.loc[(junctions_bed.start >= coords_start) & (junctions_bed.end <= coords_end)]
    
    # Return empty DataFrame if no junctions are found for the gene ID
    if junctions_bed.shape[0] == 0:
        junctions_bed = BedTool.from_dataframe(junctions_bed)
        return junctions_bed
    
    


    
    
    # Generate unique junction names
    junctions_bed['junc_names'] = ['junction_' + str(i) for i in range(1, junctions_bed.shape[0]+1)]
    
    # Rearrange DataFrame columns
    junctions_bed = junctions_bed[['chrom', 'start', 'end', 'junc_names', 'junc_counts']]
    
    # Convert DataFrame to BedTool object
    junctions_bed = BedTool.from_dataframe(junctions_bed)
    
    return junctions_bed



def get_intron_chain(transcript_bed):
    if pd.DataFrame(transcript_bed).shape[1] >= 3:
        intron_chain = []
        for x in zip(transcript_bed.iloc[:-1].end, transcript_bed.iloc[1:].start):
            intron_chain.append(x)
    else:
        intron_chain = None
        
    return intron_chain

def get_exon_chain(transcript_bed):
    exon_chain = []
    for idx, row in transcript_bed.iterrows():
        exon_chain.append((int(row.start), int(row.end)))

    return exon_chain



def is_same_chain(chain_a, chain_b):
    '''
    Test if chain_a is the same as chain_b
    '''
    return chain_a == chain_b
    

def same_terminal_exons(exon_chain_1, exon_chain_2, strand,
                       three_prime_diff = 200,
                       three_prime_ratio = 0.75,
                       five_prime_diff = 500,
                       five_prime_ratio = 0.5):

    # if (len(exon_chain_1) == 0) or (len(exon_chain_2) == 0):
    #     return True
    if strand == '+':
        last_exon_1 = exon_chain_1[-1]
        last_exon_2 = exon_chain_2[-1]


        ###### THIS IS WHAT CHANGED

        if (len(exon_chain_1) > 1) or (len(exon_chain_2) > 1):

            if last_exon_1[0] != last_exon_2[0]:
                # print('Different chain')
                return False

        first_exon_1 = exon_chain_1[0]
        first_exon_2 = exon_chain_2[0]

        if (len(exon_chain_1) > 1) or (len(exon_chain_2) > 1):

            if first_exon_1[1] != first_exon_2[1]:
                # print('Different chain')
                return False

        ##############################

    else:
        last_exon_1 = exon_chain_1[0]
        last_exon_2 = exon_chain_2[0]

        if last_exon_1[1] != last_exon_2[1]:
            # print('Different chain')
            return False

        first_exon_1 = exon_chain_1[-1]
        first_exon_2 = exon_chain_2[-1]

        if first_exon_1[0] != first_exon_2[0]:
            # print('Different chain')
            return False
        

    last_exon_1_len = last_exon_1[1] - last_exon_1[0]
    last_exon_2_len = last_exon_2[1] - last_exon_2[0]

    long_exon = np.max([last_exon_1_len, last_exon_2_len])
    short_exon = np.min([last_exon_1_len, last_exon_2_len])

    exon_diff = long_exon - short_exon
    exon_ratio = short_exon/long_exon

    if (exon_ratio < three_prime_ratio) and (exon_diff > three_prime_diff):
        # print('Alt. UTR')
        return False

    

    first_exon_1_len = first_exon_1[1] - first_exon_1[0]
    first_exon_2_len = first_exon_2[1] - first_exon_2[0]

    long_exon = np.max([first_exon_1_len, first_exon_2_len])
    short_exon = np.min([first_exon_1_len, first_exon_2_len])

    exon_diff = long_exon - short_exon
    exon_ratio = short_exon/long_exon

    if (exon_ratio < five_prime_ratio) and (exon_diff > five_prime_diff):
        # print('Alt. initiation')
        return False

    return True

def is_same_isoform(isoform_1, isoform_2, strand,
                       three_prime_diff = 200,
                       three_prime_ratio = 0.75,
                       five_prime_diff = 500,
                       five_prime_ratio = 0.5):
    
    exon_chain_1 = get_exon_chain(isoform_1)
    exon_chain_2 = get_exon_chain(isoform_2)

    intron_chain_1 = get_intron_chain(isoform_1)
    intron_chain_2 = get_intron_chain(isoform_2)

    
    same_chain = is_same_chain(intron_chain_1, intron_chain_2)

    same_terminal = same_terminal_exons(exon_chain_1, exon_chain_2, strand,
                                       three_prime_diff,
                                        three_prime_ratio,
                                        five_prime_diff,
                                        five_prime_ratio
                                       )

    if same_chain and same_terminal:
        return True
    else:
        return False


# def is_subchain(isoform_1, isoform_2, strand):

#     if len(isoform_1) > len(isoform_2):
#         long_isoform = isoform_1
#         short_isoform = isoform_2

#     else:
#         long_isoform = isoform_2
#         short_isoform = isoform_1

#     if strand == '+':
#         long_isoform_sub = long_isoform.iloc[-len(short_isoform):]
#     else:
#         long_isoform_sub = long_isoform.iloc[:len(short_isoform)]

#     if len(short_isoform) > 1:

#         # print('Test if short isoform is sub chain of long isoform')

#         is_isoform_subchain = is_same_isoform(long_isoform_sub, short_isoform, strand,
#                                               five_prime_diff = 200,
#                                               five_prime_ratio = 0.75
#                                              )

#     else:
#         # print('Short isoform is single exon.')
#         is_isoform_subchain = last_exon_match(long_isoform, short_isoform, strand)

#     return is_isoform_subchain

###########################################


def is_subchain(isoform_1, isoform_2, strand):
    if len(isoform_1) > len(isoform_2):
        long_isoform = isoform_1
        short_isoform = isoform_2
    else:
        long_isoform = isoform_2
        short_isoform = isoform_1
    
    # Special case for single-row isoform_2 with "Old approach" in notes
    if (len(isoform_2) == 1) and ('notes' in isoform_2.columns):
        notes_value = isoform_2.iloc[0]['notes']
        if isinstance(notes_value, str) and "Old approach" in notes_value:
            # Compare with appropriate end of isoform_1 based on strand
            if strand == '+':
                target_entry = isoform_1.iloc[-1]  # last entry
                long_isoform_sub = long_isoform.iloc[-len(short_isoform):]
            else:
                target_entry = isoform_1.iloc[0]   # first entry
                long_isoform_sub = long_isoform.iloc[:len(short_isoform)]
            
            # Check if they match 75% (you'll need to implement this matching logic)
            # This is a placeholder - replace with your actual matching function

            match_75 = matches_75_percent(target_entry, isoform_2.iloc[0])


            is_isoform_subchain = last_exon_match(long_isoform, short_isoform, strand)
            
            return match_75 or is_isoform_subchain #matches_75_percent(target_entry, isoform_2.iloc[0])
    
    if strand == '+':
        long_isoform_sub = long_isoform.iloc[-len(short_isoform):]
    else:
        long_isoform_sub = long_isoform.iloc[:len(short_isoform)]
    if len(short_isoform) > 1:
        # print('Test if short isoform is sub chain of long isoform')
        is_isoform_subchain = is_same_isoform(long_isoform_sub, short_isoform, strand,
                                              five_prime_diff = 200,
                                              five_prime_ratio = 0.75
                                             )
    else:
        # print('Short isoform is single exon.')
        is_isoform_subchain = last_exon_match(long_isoform, short_isoform, strand)
    return is_isoform_subchain


def matches_75_percent(entry1, entry2):
    """
    Check if two entries match within 75% criteria.
    The absolute difference between start and end coordinates should be 
    less than 25% of the total length of entry2.
    """
    # Calculate the length of entry2 (isoform_2's single row)
    entry2_length = entry2['end'] - entry2['start']
    
    # Calculate 25% threshold
    threshold = 0.25 * entry2_length
    
    # Check if both start and end differences are within threshold
    start_diff = abs(entry1['start'] - entry2['start'])
    end_diff = abs(entry1['end'] - entry2['end'])
    
    return start_diff < threshold and end_diff < threshold


def last_exon_match(long_isoform, short_isoform, strand):
    exon_chain_1 = get_exon_chain(long_isoform)
    exon_chain_2 = get_exon_chain(short_isoform)
    
    if strand == '+':
        exon_long_iso = exon_chain_1[-1]
        exon_short_iso = exon_chain_2[0]

        len_diff = np.abs(exon_long_iso[1] - exon_short_iso[1])
        
        if len(long_isoform) == 1:
            
            if exon_long_iso[0] != exon_short_iso[0]:
                # ("Freestanding isoforms have to match at 5' end.")
                return False
        else:
            if exon_short_iso[0] < exon_long_iso[0]:
                # ("Freestanding isoforms have to match at 5' end of the last exon.")
                return False
            
    else:
        exon_long_iso = exon_chain_1[0]
        exon_short_iso = exon_chain_2[0]

        len_diff = np.abs(exon_long_iso[0] - exon_short_iso[0])
        if len(long_isoform) == 1:
            if exon_long_iso[1] != exon_short_iso[1]:
                # ("Freestanding isoforms have to match at 5' end.")
                return False
        else:
            if exon_short_iso[1] > exon_long_iso[1]:
                # ("Freestanding isoforms have to match at 5' end of the last exon.")
                return False
                    
    len_exon_1 = exon_long_iso[1] - exon_long_iso[0]
    len_exon_2 = exon_short_iso[1] - exon_short_iso[0]


    if len_exon_2 > len_exon_1:
        exon_ratio = len_diff/len_exon_1
    else:
        exon_ratio = len_diff/len_exon_2

    if (len_diff > 200) and (exon_ratio > 0.75):
        return False

    return True

# def merge_factors_into_isoforms(isoforms, strand):
#     factor_list = sorted(isoforms.keys())
#     merge_dict = {}
#     accounted_factor = []

#     # Check for isoforms that are identical
    
#     for i, factor_1 in enumerate(factor_list):
        
#         if factor_1 in accounted_factor:
#             continue
#         merge_dict.update({factor_1:[]})
#         for factor_2 in factor_list[i+1:]:
#             if factor_1 in accounted_factor:
#                 continue
            
#             same_isoform = is_same_isoform(isoforms[factor_1], isoforms[factor_2], strand)
#             # subchain_iso = is_subchain(isoforms[factor_1], isoforms[factor_2], strand)

#             if same_isoform:# or subchain_iso:
#                 merge_dict[factor_1].append(factor_2)
#                 accounted_factor.append(factor_2)
#         accounted_factor.append(factor_1)

#     # Now we account for isoforms that are sub-isoforms
#     # We do this separately so that we have a pre-defined hierarchy of isoforms to merge
    
#     candidate_isoforms = {}
#     accounted_factor = []
#     for i, factor_1 in enumerate(sorted(merge_dict.keys())):
#         if factor_1 in accounted_factor:
#             continue
#         candidate_isoforms.update({factor_1:merge_dict[factor_1]})
#         for factor_2 in sorted(merge_dict.keys()):
#             if factor_1 == factor_2:
#                 continue
#             if len(isoforms[factor_1]) > len(isoforms[factor_2]):
#                 # Keep proper order by only merging small to large
#                 is_sub = is_subchain(isoforms[factor_1], isoforms[factor_2], strand)
#                 if is_sub:
#                     candidate_isoforms[factor_1].append(factor_2)
#                     candidate_isoforms[factor_1].extend(merge_dict[factor_2])
#                     accounted_factor.append(factor_2)
                    

#     # for factor_1 in factor_list:
#     #     merge_dict[factor_1] = sorted(set(merge_dict[factor_1]))

#     # print(accounted_factor)
#     for factor in accounted_factor:
#         if factor in candidate_isoforms.keys():
#             candidate_isoforms.pop(factor)

#     merged_isoforms = deduplicate_by_longest_list(candidate_isoforms)
#     return merged_isoforms

# def deduplicate_by_longest_list(input_dict):
#     # Create a dictionary to track which lists each element appears in
#     element_occurrences = {}
#     for key, lst in input_dict.items():
#         for item in lst:
#             element_occurrences.setdefault(item, []).append(key)

#     # Sort keys by length of their list (descending)
#     sorted_keys = sorted(input_dict, key=lambda k: len(input_dict[k]), reverse=True)

#     # Create a copy to modify
#     output_dict = {k: [] for k in input_dict}

#     # Track which elements have already been assigned
#     assigned_elements = set()

#     for key in sorted_keys:
#         for item in input_dict[key]:
#             if item not in assigned_elements:
#                 output_dict[key].append(item)
#                 assigned_elements.add(item)

#     return output_dict





#############################################

def merge_factors_into_isoforms(isoforms, strand):
    factor_list = sorted(isoforms.keys())
    merge_dict = {}
    accounted_factor = []
    # Check for isoforms that are identical
    
    for i, factor_1 in enumerate(factor_list):
        
        if factor_1 in accounted_factor:
            continue
        merge_dict.update({factor_1:[]})
        for factor_2 in factor_list[i+1:]:
            if factor_1 in accounted_factor:
                continue
            
            same_isoform = is_same_isoform(isoforms[factor_1], isoforms[factor_2], strand)
            # subchain_iso = is_subchain(isoforms[factor_1], isoforms[factor_2], strand)
            if same_isoform:# or subchain_iso:
                merge_dict[factor_1].append(factor_2)
                accounted_factor.append(factor_2)
        accounted_factor.append(factor_1)
    # Now we account for isoforms that are sub-isoforms
    # We do this separately so that we have a pre-defined hierarchy of isoforms to merge
    
    candidate_isoforms = {}
    accounted_factor = []
    for i, factor_1 in enumerate(sorted(merge_dict.keys())):
        if factor_1 in accounted_factor:
            continue
        candidate_isoforms.update({factor_1:merge_dict[factor_1]})
        for factor_2 in sorted(merge_dict.keys()):
            if factor_1 == factor_2:
                continue
            if len(isoforms[factor_1]) > len(isoforms[factor_2]):
                # Keep proper order by only merging small to large
                is_sub = is_subchain(isoforms[factor_1], isoforms[factor_2], strand)
                if is_sub:
                    candidate_isoforms[factor_1].append(factor_2)
                    candidate_isoforms[factor_1].extend(merge_dict[factor_2])
                    accounted_factor.append(factor_2)
                    
    # for factor_1 in factor_list:
    #     merge_dict[factor_1] = sorted(set(merge_dict[factor_1]))
    # print(accounted_factor)
    for factor in accounted_factor:
        if factor in candidate_isoforms.keys():
            candidate_isoforms.pop(factor)
    merged_isoforms = deduplicate_by_longest_isoform(candidate_isoforms, isoforms)
    return merged_isoforms


def deduplicate_by_longest_isoform(input_dict, isoforms):
    # Create a dictionary to track which lists each element appears in
    element_occurrences = {}
    for key, lst in input_dict.items():
        for item in lst:
            element_occurrences.setdefault(item, []).append(key)
    
    # Sort keys by number of rows in their isoform dataframe (descending)
    sorted_keys = sorted(input_dict, key=lambda k: len(isoforms[k]), reverse=True)
    
    # Create a copy to modify
    output_dict = {k: [] for k in input_dict}
    # Track which elements have already been assigned
    assigned_elements = set()
    
    for key in sorted_keys:
        for item in input_dict[key]:
            if item not in assigned_elements:
                output_dict[key].append(item)
                assigned_elements.add(item)
    return output_dict


###################################################################

def factors_to_isoforms(isoforms, merged_isoforms_dict):
    isoform_list = []
    iso = 1
    for factor in sorted(merged_isoforms_dict.keys()):
        iso_name = f'isoform_{str(iso)}'
        iso += 1

        merge_list = [factor] + merged_isoforms_dict[factor]

        start = isoforms[factor].start.min()
        end = isoforms[factor].end.max()

        notes = ''

        for sub_factor in merge_list:
            notes += isoforms[sub_factor].notes.iloc[0]
            start_alt = isoforms[sub_factor].start.min()
            end_alt = isoforms[sub_factor].end.max()

            if start_alt < start:
                start = start_alt
            if end_alt > end:
                end = end_alt
            


        
        isoform = isoforms[factor].copy()
        isoform['start'] = [start] + list(isoform.start[1:])
        isoform['end'] = list(isoform.end[:-1]) + [end]
        isoform['transcript_id'] = iso_name
        factors_id = ':'.join(merge_list)
        isoform['factors'] = [factors_id]*len(isoform)
        isoform['exon_id'] = [f'exon_{str(x)}' for x in range(1, len(isoform)+1)]
        isoform['notes'] = [notes]*len(isoform)

        isoform = isoform[['chrom', 'start', 'end', 'gene_id', 'transcript_id', 'strand', 'factors', 'exon_id', 'notes']]

        isoform_list.append(isoform)

    return isoform_list


def process_isoforms(gene, species, genes_file, junc_file, use_low_juncs = False):
    strand = genes_file.loc[genes_file.gene_name == gene].strand.iloc[0]
    if species == 'Human':
        # junc_file = '/project2/mstephens/cfbuenabadn/gtex-stm/code/junctions.tab.gz'
        gene_id = genes_file.loc[genes_file.gene_name == gene].gene_id.iloc[0]
    else:
        gene_id = '.'
        # junc_file = f'results/leafcutter2/{species}/junctions.sorted.tab.gz'
        
    isoforms = get_unmerged_isoforms(gene, junc_file, species, gene_id = gene_id, strand = strand, use_low_juncs=use_low_juncs)
    if isoforms is None:
        return None, None
    merged_isoforms = pd.concat(factors_to_isoforms(isoforms, merge_factors_into_isoforms(isoforms, strand)))

    return merged_isoforms, isoforms


    # merged_isoforms
