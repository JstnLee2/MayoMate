import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import logging
import scipy.stats as stats

logging.basicConfig(
    handlers = [logging.FileHandler(
        filename="meiotic_frequencies.log",
        encoding = 'utf-8',
        mode = 'w'
        )],
    level=logging.INFO,
    format='%(levelname)s:.%(funcName)s: %(message)s',
    force = True
    )


def main():
    """
    Docstring for main
    """

    genotypes_to_plot = ["ung1∆ A3A", "ung1∆msh2∆ A3A", "ung1∆mlh1∆ A3A", "ung1∆SNPless A3A", "ung1∆pol3-01 A3A", "ung1∆pol3-5DV A3A", "ung1∆pol32∆exo1-nd,pif1-m2 A3A", "ung1∆msh2SNPless A3A"] # EV data: "ung1∆ EV", "ung1∆msh2∆ EV", "ung1∆mlh1∆ EV", "ung1∆SNPless EV", "ung1∆pol3-01 EV", "ung1∆pol3-5DV EV", "ung1∆pol32∆exo1-nd,pif1-m2 EV", "ung1∆msh2SNPless EV"       A3A data: "ung1∆ A3A", "ung1∆msh2∆ A3A", "ung1∆mlh1∆ A3A", "ung1∆SNPless A3A", "ung1∆pol3-01 A3A", "ung1∆pol3-5DV A3A", "ung1∆pol32∆exo1-nd,pif1-m2 A3A", "ung1∆msh2SNPless A3A"

    file_path = "meiotic_frequencies_mayomate.csv"

    
    plot_meiotic_frequencies(file_path, genotypes_to_plot)



    
def pull_df_data(file_name: str) -> pd.DataFrame:
    """
    Pulls data from a CSV file and returns a pandas Dataframe.

    Parameters:
    - file_name (str): The path of the CSV file to read.
    Returns:
    - pd.df: Dataframe containing the data from the CSV file.
    """
    data = pd.read_csv(file_name)
    logging.info(f"Data pulled from {file_name} successfully.")
    return data

def build_genotype_dict(file_path: str, genotypes_to_plot: list[str]) -> dict:
    """
    This function takes the dataframe from file_path and groups the data "pre_meiotic" and "post_meiotic" by genotype.

    Parameters:
    - data_df (pd.DataFrame): Dataframe containing the data.
    
    Returns:
    - dict: A dictionary of grouped data by genotype.
    """
    if not file_path.endswith('.csv'):
        raise ValueError("Error: Input file must be a CSV file to build genotype dictionary.")
    
    df = pull_df_data(file_path)    # Load the data from CSV

    
    genotype_dict: {str, list(int), list(int)} = {}  # Initialize an empty dictionary to hold the grouped data, with genotype as keys and lists of premeiotic and postmeiotic frequencies as values.
    for _, row in df.iterrows(): # iterate through each row of the dataframe. _ indicates that the index is not used in the loop, but iterates through all the rows.
        genotype = row['sample'] # extracts the genotype from the row in which the column is named 'sample'.
        if genotype not in genotypes_to_plot: # if the genotype is not found in the genotypes to plot, it is skipped and the loop continues without 
            continue                        # executing the rest of the loop.
        freq_data = {          # creates a dictionary to hold the two sets of data, premeiotic and postmeiotic. While iterating through the rows, it extracts the values for premeiotic and postmeiotic frequencies.
            'premeiotic': row['premeiotic_frequencies'],
            'postmeiotic': row['postmeiotic_frequencies']
        }

        if genotype not in genotype_dict: # if the genotype is not already a key in the genotype_dict, it initializes a new entry for that genotype.
            genotype_dict[genotype] = {'premeiotic': [], 'postmeiotic': []}

        for freq_type in ['premeiotic', 'postmeiotic']: # iterates throught the two frequency types.
            freq = freq_data[freq_type] 
            if freq == 'NA' or freq == '': # if the frequency is empty or NA, it gets set to NaN.
                freq = np.nan
                logging.info(f"{freq_type.capitalize()} frequency for genotype {genotype} set to NaN.")
            if not np.isnan(freq): # if the frequency is not NaN, it gets converted to an integer and multiplied by 10**10.
                value = int(freq * 10**10)
            else:
                value = np.nan 
            genotype_dict[genotype][freq_type].append(value) # after the value has been updated, it gets appended to the key [genotype] and the respective frequency type.
            logging.info(f"Appended premeiotic {freq_data['premeiotic']} Postmeiotic {freq_data['postmeiotic']} To genotype {genotype}.") # loggins after each data point is appended.

        

    return genotype_dict

def conduct_mann_whitney_test(genotype_dict: dict) -> tuple[dict, dict]:
    """
    Conducts Mann-Whitney U test between geneotypes and ung1∆ for premeiotic and postmeiotic frequencies.

    Parameters:
    - genotype_dict (dict): A dictionary of grouped data by genotype.

    Returns:
    - tuple[dict, dict]: A tuple containing the genotype_dict and results dictionary.
    """

    results = {}
    genotypes = list(genotype_dict.keys())
    
    # only compare each genotype to ung1, not all pairwise comparisons
    j = genotypes.index('ung1∆ A3A') # "ung1∆ EV"
    for i in range(len(genotypes)):
        genotype1 = genotypes[i]
        genotype2 = genotypes[j]

        pre_values1 = pd.to_numeric(genotype_dict[genotype1]['premeiotic'], errors='coerce')
        pre_values2 = pd.to_numeric(genotype_dict[genotype2]['premeiotic'], errors='coerce')
        post_values1 = pd.to_numeric(genotype_dict[genotype1]['postmeiotic'], errors='coerce')
        post_values2 = pd.to_numeric(genotype_dict[genotype2]['postmeiotic'], errors='coerce')

        pre_stat, pre_p_value = stats.mannwhitneyu(pre_values1[~np.isnan(pre_values1)], pre_values2[~np.isnan(pre_values2)], alternative='two-sided')
        post_stat, post_p_value = stats.mannwhitneyu(post_values1[~np.isnan(post_values1)], post_values2[~np.isnan(post_values2)], alternative='two-sided')

        results[f"{genotype1}"] = {f"{genotype1}": {
            'premeiotic': {
                'p_value': pre_p_value,
                f'median_genotype_{genotype1}': np.median(pre_values1),
                'median_genotype2': np.median(pre_values2),
                f'range_{genotype1}': (float(np.min(pre_values1)), float(np.max(pre_values1))) if len(pre_values1) > 0 else (np.nan, np.nan)
            },
            'postmeiotic': {
                'p_value': post_p_value,
                f'median_genotype_{genotype1}': np.median(post_values1),
                'median_genotype2': np.median(post_values2),
                f'range_{genotype1}': (float(np.min(post_values1)), float(np.max(post_values1))) if len(post_values1) > 0 else (np.nan, np.nan)
            }
        }}
        logging.info(f"Conducted Mann-Whitney U test between {genotype1} and {genotype2} for premeiotic and postmeiotic frequencies.")
        logging.info(f"Results: {results[f'{genotype1}']}")

    return results


# def meiotic_stats(file_path: str, genotypes_to_plot: list[str]) -> dict:
#     """
#     Calculates the median and range for each condition in the dataframe and saves it into a dictionary.

#     Parameters:
#     - data_df (pd.DataFrame): Dataframe containing the data.
    
#     Returns:
#     - dict: A dictionary of meiotic statistics for each genotype.
#     """
#     if not file_path.endswith('.csv'):
#         raise ValueError("Error: Input file must be a CSV file to build genotype dictionary.")
    

#     get_data = build_genotype_dict(file_path, genotypes_to_plot) # get the genotype dictionary
#     logging.info(f"Genotype data retrieved for meiotic stats calculation.")

#     meiotic_stats = {} # initialize an empty dictionary to hold the statistics
#     for genotype, freq_dict in get_data.items(): # iterates through each genotype and freq_dict represents the lists associated with that genotype
#         pre_values = pd.to_numeric(freq_dict['premeiotic'], errors='coerce') # first list in freq_dict is premeiotic frequencies. .dropna() Removes NaN values and associated data keeping complete data only.
#         post_values = pd.to_numeric(freq_dict['postmeiotic'], errors='coerce') # second list in freq_dict is postmeiotic frequencies. .dropna()Removes NaN values and associated data keeping complete data only.

#         pre_stats = {   # creates a dictionary to hold the all statistics for premeiotic frequencies
#             'median': float(np.median(pre_values)) if len(pre_values) > 0 else np.nan, # 1 line code calculating the median of premeiotic frequencies if there are more than 0 values, else NaN
#             'range': (float(np.min(pre_values)), float(np.max(pre_values))) if len(pre_values) > 0 else (np.nan, np.nan), # 1 line code calculating the range of premeiotic frequencies if there are more than 0 values, else (NaN, NaN)
#             'count': int(np.count_nonzero(~np.isnan(pre_values)))
#         }
#         post_stats = {  # creates a dictionary to hold the all statistics for postmeiotic frequencies
#             'median': float(np.median(post_values)) if len(post_values) > 0 else np.nan,
#             'range': (float(np.min(post_values)), float(np.max(post_values))) if len(post_values) > 0 else (np.nan, np.nan),
#             'count': int(np.count_nonzero(~np.isnan(post_values)))
#         }
#         meiotic_stats[genotype] = { # creates a dictionary to hold both premeiotic and postmeiotic statistics for each genotype
#             'premeiotic': pre_stats,
#             'postmeiotic': post_stats
#         }
#         logging.info(f"""Calculated stats for genotype {genotype}:
#                     Premeiotic={pre_stats} 
#                     Postmeiotic={post_stats}"""
#                     )
    
#     logging.info(f'meiotic_stats completed successfully.')
#     return meiotic_stats, genotype

def identify_significant_comparisons(mann_whitney_input: dict, genotype_dict: dict, alpha: float = 0.05) -> dict:
    """
    Identifies significant comparisons from the Mann-Whitney U test results.

    Parameters:
    - mann_whitney_results (dict): A dictionary containing the results of the Mann-Whitney U tests.
    - alpha (float): The significance level to use for determining significance. Default is 0.05.

    Returns:
    - tuple[dict, dict]: A tuple containing the significant comparisons and their p-values, and the genotype_dict.
    """

    significant_comparisons = {}
    
    for comparison, results in mann_whitney_input.items():
        for condition in ['premeiotic', 'postmeiotic']:
            p_value = results[condition]['p_value']
            if p_value < alpha:
                if comparison not in significant_comparisons:
                    significant_comparisons[comparison] = {}
                significant_comparisons[comparison][condition] = p_value
                logging.info(f"Significant comparison identified: {comparison} ({condition}) with p-value {p_value:.4f}")

    return significant_comparisons


def plot_meiotic_frequencies(file_path: str, 
                             genotypes_to_plot: list[str],
                             save = True,
                             save_path = "meiotic_frequencies_plot.png"
                             ) -> None:
    """
    Plots the meiotic frequencies from the CSV file.

    Parameters:
    - file_path (str): The path of the CSV file to read.
    - save (bool): Whether to save the plot as a PNG file. Default is True.
    - save_path (str): The path to save the PNG file. Default is "meiotic_frequencies_plot.png".

    Returns: 
    None
    """
    if not file_path.endswith('.csv'):
        raise ValueError("Error: Input file must be a CSV file to plot meiotic frequencies.")
    
    genotype_dict = build_genotype_dict(file_path, genotypes_to_plot) # build the genotype dictionary for plotting
    
    stats = conduct_mann_whitney_test(genotype_dict) # get the meiotic statistics
    logging.info(f"Meiotic statistics retrieved for plotting.")

  
    medians = {}
    ranges = {}

    # Extract medians and ranges from the statistics dictionary produced by conduct_mann_whitney_test
    for genotype in stats.keys():
        medians[genotype] = {
            'premeiotic': stats[genotype][genotype]['premeiotic'][f'median_genotype_{genotype}'],
            'postmeiotic': stats[genotype][genotype]['postmeiotic'][f'median_genotype_{genotype}']
        }
        ranges[genotype] = {
            'premeiotic': stats[genotype][genotype]['premeiotic'][f'range_{genotype}'],
            'postmeiotic': stats[genotype][genotype]['postmeiotic'][f'range_{genotype}']
        }

    logging.info(f"Medians and ranges extracted for plotting. \nMedians: {medians} \nRanges: {ranges}")

    # Set the colors (Blue for premeiotic, Red for meiotic)
    colors = ['#3274a1', '#cc3333', '#3274a1', '#cc3333']

    # --- 2. PLOTTING ---
    fig, ax = plt.subplots(figsize=(24, 14))

    x_pos = np.arange(len(labels))
    bars = ax.bar(x_pos, medians, color=colors, edgecolor='white', width=0.8)

    # Set Y-axis to Logarithmic
    ax.set_yscale('log')
    ax.set_ylim(10, max([r[1] for r in ranges]) * 2) # Adjust based on your data spread

    # Styling the axis and labels
    ax.set_ylabel('ura3-29 reversion frequency ($10^{10}$)', fontsize=14, fontweight='bold')
    ax.set_title('ura3-29 reversion frequencies', fontsize=18, fontweight='bold', pad=20)
    ax.set_xticks(x_pos)
    labels = [label.replace(" ", "\n") for label in labels]  # Replace spaces with newlines for better readability
    ax.set_xticklabels(labels, fontsize=10, rotation=60, fontweight='bold', fontstyle='italic')

    # Gridlines
    ax.yaxis.grid(True, which='major', linestyle='-', linewidth=0.5, alpha=0.7)
    ax.set_axisbelow(True)

    genotype_dict = build_genotype_dict(file_path, genotypes_to_plot)
    mann_whitney_input = conduct_mann_whitney_test(genotype_dict)
    significant_comparisons = identify_significant_comparisons(mann_whitney_input, genotype_dict)
    print("=" *80)
    print(significant_comparisons)
    print("=" *80)

    # --- 3. ADDING TEXT OVERLAYS ---
    for i, bar in enumerate(bars):
        height = bar.get_height()
        
        # Add Median Label (White box with black border)
        ax.text(bar.get_x() + bar.get_width()/2., height * 1.3,
                f'{int(medians[i])}',
                ha='center', va='center', fontweight='bold', fontsize=8,
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
        
        # Add Range Label (Top text)
        low, high = ranges[i]
        ax.text(bar.get_x() + bar.get_width()/2., height * .7,
                f'({int(low)} - {int(high)})',
                ha='center', va='bottom', fontsize=6, fontweight='bold')
        
        # If conduct_mann_whitney_test results are significant, add asterisk
        # Need to conduct the Mann-Whitney U test and check the p-values to determine significance
        
        # Example: If significant, add an asterisk above the bar
        # If genotype is in the comparison and condition matches, add asterisk
    
        if add_asterisk_to_significant_comparisons(significant_comparisons, genotype_dict):
            ax.text(bar.get_x() + bar.get_width()/2., height * 1.0,
                        '*',
                        ha='center', va='bottom', fontsize=50, color='black', fontweight='bold')
            

    plt.tight_layout()
    if save:
        plt.savefig(save_path, dpi=300)
        logging.info(f"Meiotic frequencies plot saved to {save_path}.")
    plt.show()

def add_asterisk_to_significant_comparisons(significant_comparisons: dict, genotype_dict: dict) -> tuple:
    """
    Adds an asterisk to the plot for significant comparisons.

    Parameters:
    - significant_comparisons (dict): A dictionary containing the significant comparisons and their p-values.

    Returns:
    - bool: True if an asterisk was added, False otherwise.
    """
    # This function is a placeholder for the logic to add asterisks to the plot based on significant comparisons.
    # The actual implementation would depend on how the significant comparisons are structured and how they relate to the bars in the plot.
    for genotype in genotype_dict.keys():
        if genotype in significant_comparisons:
            # Logic to add asterisk to the plot
            # Add asterisk to plot
            logging.info(f"Genotype {genotype} is statistically different to ung1∆")
            return genotype, True
    return False
    
    

if __name__ == "__main__":
    main()