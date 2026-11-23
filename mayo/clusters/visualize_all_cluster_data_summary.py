# This code was developed and authored by Jerzy Twarowski in Malkova Lab at the University of Iowa 
# Contact: jerzymateusz-twarowski@uiowa.edu, tvarovski1@gmail.com
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import plotly.express as px
import logging
from typing import Literal

from mayo.clusters import fill_in_missing_samples


def plot_pairplot_plt(df: pd.DataFrame, columns_to_plot: list, save_name:str = 'pairplot.png', save:bool = True) -> None:
    """
    Creates a pairplot of the given columns in the given DataFrame, grouped by genotype.

    Parameters:
        df (pandas.DataFrame): The DataFrame containing the data to plot.
        columns_to_plot (list): A list of column names to include in the pairplot.
        save_name (str): The name to use when saving the plot. Defaults to 'pairplot.png'.
        save (bool): Whether or not to save the plot. Defaults to True.

    Returns:
        None
    """
    sns.pairplot(df, hue='genotype', vars=columns_to_plot, diag_kind='kde', markers="+")
    if save:
        plt.savefig(f"figures/clusters/{save_name}", dpi=500)
    plt.show()

def conduct_ttest(set1: pd.DataFrame, set2: pd.DataFrame, column: str) -> float:
    """
    Conducts a two-sample t-test on the given column of two DataFrames.

    Parameters:
        set1 (pandas.DataFrame): The first DataFrame to use in the t-test.
        set2 (pandas.DataFrame): The second DataFrame to use in the t-test.
        column (str): The name of the column to use in the t-test.

    Returns:
        float: The p-value resulting from the t-test.
    """
    tstat, pval = (stats.ttest_ind(set1[column], set2[column]))
    print(f"t-test: p-value = {pval:.5e}, t-statistic = {tstat:.5e}")
    return pval
   
def conduct_mannutest(set1: pd.DataFrame, set2: pd.DataFrame, column: str) -> float:
    """
    Conducts a Mann-Whitney U test on the given column of two DataFrames.

    Parameters:
        set1 (pandas.DataFrame): The first DataFrame to use in the test.
        set2 (pandas.DataFrame): The second DataFrame to use in the test.
        column (str): The name of the column to use in the test.

    Returns:
        float: The p-value resulting from the test.
    """
    u, pval_u = stats.mannwhitneyu(set1[column], set2[column])
    print(f"Mann-Whitney U test: p-value = {pval_u:.5e}, U = {u:.5e} ({set1[column].mean():.3} (N:{len(set1)}) vs. {set2[column].mean():.3} (N:{len(set2)}))")
    return pval_u

def conduct_f_oneway(set1: pd.DataFrame, set2: pd.DataFrame, column: str) -> float:
    u, pval_f = (stats.f_oneway(set1[column], set2[column]))
    print(f"ANOVA: p-value = {pval_f:.5e}, F = {u:.5e}")
    return pval_f

def conduct_kurtosis_test(set1: pd.DataFrame, set2: pd.DataFrame, column: str, genotype1: str, genotype2: str) -> tuple:
    """
    Conducts a kurtosis test on the given column of two DataFrames.

    Parameters:
        set1 (pandas.DataFrame): The first DataFrame to use in the test.
        set2 (pandas.DataFrame): The second DataFrame to use in the test.
        column (str): The name of the column to use in the test.

    Returns:
        float: The p-value resulting from the test.
    """
    u, pval_k = (stats.kurtosistest(set1[column]))
    print(f"Kurtosis test {genotype1}: p-value = {pval_k:.5e}, U = {u:.5e}")
    u2, pval_k2 = (stats.kurtosistest(set2[column]))
    print(f"Kurtosis test {genotype2}: p-value = {pval_k2:.5e}, U = {u2:.5e}")

    return pval_k, pval_k2

def conduct_bartlet_test(set1: pd.DataFrame, set2: pd.DataFrame, column: str) -> float:
    """
    Conducts a Hartley's test on the given column of two DataFrames.

    Parameters:
        set1 (pandas.DataFrame): The first DataFrame to use in the test.
        set2 (pandas.DataFrame): The second DataFrame to use in the test.
        column (str): The name of the column to use in the test.

    Returns:
        float: The p-value resulting from the test.
    """
    u, pval_h = (stats.bartlett(set1[column], set2[column]))
    print(f"Hartley's test: p-value = {pval_h:.5e}, U = {u:.5e}")
    return pval_h

def conduct_levene_test(set1: pd.DataFrame, set2: pd.DataFrame, column: str) -> float:
    """
    Conducts a Levene's test on the given column of two DataFrames.

    Parameters:
        set1 (pandas.DataFrame): The first DataFrame to use in the test.
        set2 (pandas.DataFrame): The second DataFrame to use in the test.
        column (str): The name of the column to use in the test.

    Returns:
        float: The p-value resulting from the test.
    """
    u, pval_l = (stats.levene(set1[column], set2[column]))
    print(f"Levene's test: p-value = {pval_l:.5e}, U = {u:.5e}")
    return pval_l

def conduct_2_sample_kolmogorov_smirnov_test(set1: pd.DataFrame, set2: pd.DataFrame, column: str) -> float:
    """
    Conducts a two-sample Kolmogorov-Smirnov test on the given column of two DataFrames.

    Parameters:
        set1 (pandas.DataFrame): The first DataFrame to use in the test.
        set2 (pandas.DataFrame): The second DataFrame to use in the test.
        column (str): The name of the column to use in the test.

    Returns:
        float: The p-value resulting from the test.
    """
    u, pval_ks = (stats.ks_2samp(set1[column], set2[column]))
    print(f"Kolmogorov-Smirnov test: p-value = {pval_ks:.5e}, U = {u:.5e}")
    return pval_ks

def pval_to_stars(pval: float, max_stars: int = 10) -> str:
    """
    Converts a p-value to significance asterisks:
    p < 0.05 -> *, p < 0.01 -> **, p < 0.001 -> ***, p < 0.0001 -> ****, and one more
    asterisk for every further order of magnitude. Returns "ns" when p >= 0.05 (or p is NaN).
    """
    if pval is None or pd.isna(pval) or pval >= 0.05:
        return "ns"
    stars, threshold = 1, 0.01
    while pval < threshold and stars < max_stars:
        stars += 1
        threshold /= 10
    return "*" * stars

def add_significance_brackets(
        ax,
        comparisons: list,
        rename_dict: dict | None = None,
        show_ns: bool = False,
        stars_only: bool = False,
        ) -> None:
    """
    Draws a bracket above the violins for each tested pair, with significance asterisks on top.
    With stars_only=True no brackets are drawn: the asterisks go in one row above the plot,
    over the second genotype of each pair (the one compared against the reference).

    Parameters:
        ax: matplotlib Axes holding the violin plot (categorical x axis).
        comparisons (list): (genotype_a, genotype_b, pval) tuples, using the original genotype names.
        rename_dict (dict): mapping from original genotype names to the x tick labels used on the plot.
        show_ns (bool): also draw brackets labelled "ns" for non-significant pairs. Defaults to False.
        stars_only (bool): skip the brackets and only write asterisks (for busy plots). Defaults to False.
    """
    rename_dict = rename_dict or {}
    x_of = {t.get_text(): t.get_position()[0] for t in ax.get_xticklabels()}

    brackets = []
    for genotype_a, genotype_b, pval in comparisons:
        label = pval_to_stars(pval)
        if label == "ns" and not show_ns:
            continue
        xa = x_of.get(rename_dict.get(genotype_a, genotype_a))
        xb = x_of.get(rename_dict.get(genotype_b, genotype_b))
        if xa is None or xb is None or xa == xb:
            continue
        brackets.append((min(xa, xb), max(xa, xb), label, xb))
    if not brackets:
        return

    # work in "scaled" y units so brackets are evenly spaced on linear and log axes
    log_scale = ax.get_yscale() == "log"
    import numpy as np
    to_scaled = np.log10 if log_scale else (lambda v: v)
    from_scaled = (lambda v: 10 ** v) if log_scale else (lambda v: v)

    ymin, ymax = ax.get_ylim()
    s_min, s_max = to_scaled(ymin), to_scaled(ymax)
    s_span = s_max - s_min
    step = 0.12 * s_span      # vertical distance between stacked brackets
    tick = 0.02 * s_span      # height of the bracket's end ticks

    if stars_only:
        for _, _, label, x in brackets:
            ax.text(x, from_scaled(s_max), label, ha="center", va="bottom",
                    fontsize=12 if label != "ns" else 10, fontweight="bold", color="black")
        ax.set_ylim(top=from_scaled(s_max + 0.08 * s_span))
        if not log_scale:
            ax.set_yticks([t for t in ax.get_yticks() if ymin <= t <= ymax])
            ax.set_ylim(bottom=ymin)
        return

    # shortest brackets go lowest so longer ones stack over them
    brackets.sort(key=lambda b: (b[1] - b[0], b[0]))

    # each bracket goes on the lowest level that doesn't overlap a bracket already on that level
    levels = []  # list of lists of (x1, x2)
    for x1, x2, label, _ in brackets:
        level = 0
        while level < len(levels) and any(not (x2 < a or x1 > b) for a, b in levels[level]):
            level += 1
        if level == len(levels):
            levels.append([])
        levels[level].append((x1, x2))

        s_y = s_max + level * step
        y_low, y_high = from_scaled(s_y), from_scaled(s_y + tick)
        ax.plot([x1, x1, x2, x2], [y_low, y_high, y_high, y_low], color="black", lw=1.5, clip_on=False)
        ax.text((x1 + x2) / 2, y_high, label, ha="center", va="bottom",
                fontsize=14 if label != "ns" else 11, fontweight="bold", color="black")

    # make room for the brackets and their labels
    ax.set_ylim(top=from_scaled(s_max + len(levels) * step + 0.03 * s_span))

    # keep y tick labels within the data range, not in the bracket area
    if not log_scale:
        ax.set_yticks([t for t in ax.get_yticks() if ymin <= t <= ymax])
        ax.set_ylim(bottom=ymin)

def add_median_labels(ax, df: pd.DataFrame, column: str, by: str = "genotype") -> None:
    """
    Writes each group's median just above the top of its violin.

    Parameters:
        ax: matplotlib Axes holding the violin plot (categorical x axis).
        df (pandas.DataFrame): the plotted data; group names must match the x tick labels.
        column (str): the plotted column.
        by (str): the grouping column. Defaults to 'genotype'.
    """
    import numpy as np
    from matplotlib.collections import PolyCollection

    log_scale = ax.get_yscale() == "log"
    to_scaled = np.log10 if log_scale else (lambda v: v)
    from_scaled = (lambda v: 10 ** v) if log_scale else (lambda v: v)

    # top of each violin body, keyed by its x center
    violin_tops = {}
    for coll in ax.collections:
        if not isinstance(coll, PolyCollection):
            continue
        for path in coll.get_paths():
            v = path.vertices
            if len(v):
                violin_tops[round((v[:, 0].min() + v[:, 0].max()) / 2)] = v[:, 1].max()

    medians = df.groupby(by)[column].median()
    ymin, ymax = ax.get_ylim()
    s_span = to_scaled(ymax) - to_scaled(ymin)
    highest = to_scaled(ymax)
    for t in ax.get_xticklabels():
        name, x = t.get_text(), t.get_position()[0]
        if name not in medians.index or pd.isna(medians[name]):
            continue
        m = medians[name]
        top = violin_tops.get(round(x), df.loc[df[by] == name, column].max())
        s_y = to_scaled(top) + 0.02 * s_span
        # whole numbers (and large values) without decimals, small fractional values with one decimal
        label = f"{m:,.0f}" if float(m).is_integer() or abs(m) >= 100 else f"{m:,.1f}"
        ax.text(x, from_scaled(s_y), label, ha="center", va="bottom", fontsize=11, fontweight="bold", color="dimgray")
        highest = max(highest, s_y + 0.07 * s_span)  # leave room for the text itself

    if highest > to_scaled(ymax):
        ax.set_ylim(top=from_scaled(highest))

def plot_violin_plt(
        df: pd.DataFrame,
        column: str,
        save_name: str | None = "violin_plot.png",
        title: str = "",
        y_label: str = "",
        save: bool = True,
        show: bool = False,
        comparisons: list | None = None,
        show_ns: bool = False,
        show_medians: bool = True,
        max_genotypes_for_brackets: int = 8,
        **kwargs: dict
        ) -> None:
    """
    Creates a violin plot of the given column in the given DataFrame, grouped by genotype.

    Parameters:
        df (pandas.DataFrame): The DataFrame containing the data to plot.
        column (str): The name of the column to plot.
        title (str): The title to use for the plot.
        y_label (str): The label to use for the y axis.
        pval (float): The p-value to display on the plot.
        genotype1 (str): The name of the first genotype to compare.
        genotype2 (str): The name of the second genotype to compare.
        by (str): The column to group the plot by. Defaults to 'Genotype'.
        root_fig_name (str): The name to use when saving the plot. Defaults to '.png'.
        save (bool): Whether or not to save the plot. Defaults to True.
        comparisons (list): (genotype_a, genotype_b, pval) tuples; a bracket with significance
            asterisks is drawn above the violins for each one. Defaults to None (no brackets).
        show_ns (bool): also draw brackets for non-significant (p >= 0.05) pairs, labelled "ns".
        show_medians (bool): write each genotype's median above its violin. Defaults to True.
        max_genotypes_for_brackets (int): with more genotypes than this on the plot, only asterisks are
            written above the violins (no bracket lines), to keep busy plots readable. Defaults to 8.

    Returns:
    None
    """
    from matplotlib.ticker import FuncFormatter
    from betterbeeswarm import Beeswarm
    import seaborn as sns

    #rename spaces in genotype names with new lines (for better visualization)
    rename_dict = {
        'ung1∆'                         : 'ung1∆',
        'ung1∆ EV'                      : 'ung1∆ EV',
        'ung1∆msh2∆'                    : 'ung1∆msh2∆',
        'ung1∆mlh1∆'                    : 'ung1∆mlh1∆',
        'ung1∆NAT'                      : 'ung1∆NAT',
        'UNG1'                          : 'UNG1',
        'exo1-nd'                       : 'ung1∆\nexo1-nd',
        'pol32∆'                        : 'ung1∆\npol32∆',
        'exo1-ndpol32∆'                 : 'ung1∆\nexo1-nd\npol32∆',
        'ung1∆ premeiotic non-selected' : 'ung1∆\npremeiotic\nnon-selected',
        'ung1∆ non-selected'            : 'ung1∆\nnon-selected',
        'inc. Tetrad ung1∆ non-selected': 'inc. Tetrad ung1∆\nnon-selected',
        'spo13∆'                        : 'ung1∆\nspo13∆',
        'spo13∆spo11∆'                  : 'ung1∆\nspo13∆\nspo11∆',
        'exo1-ndsgs1∆C'                 : 'ung1∆\nexo1-nd\nsgs1∆C',
        'sgs1∆C'                        : 'ung1∆\nsgs1∆C',
         }

    #rename the genotype names
    df['genotype'] = df['genotype'].replace(rename_dict)


    if "order" in kwargs:
        order = kwargs["order"]
        order = [rename_dict.get(item, item) for item in order]
    else:
        order = None

    if "palette" in kwargs:
        palette = kwargs["palette"]
        color = None
        hue = 'genotype'
        palette = {rename_dict.get(k, k): v for k, v in palette.items()}
    else:
        palette = None
        color = 'tab:blue'
        hue = None

    if "figsize" in kwargs:
        figsize = kwargs["figsize"]
    else:
        figsize = (8, 6)

    num_genotypes = len(df['genotype'].unique())
    print(f"Number of genotypes: {num_genotypes}")

    plt.figure(figsize=figsize)
    sns.set_style('whitegrid')
    sns.set_context("notebook", font_scale=3)

    sns.violinplot(x='genotype', y=column, data=df, inner='quartile', color='whitesmoke', density_norm='width', order=order)
    sns.swarmplot(x='genotype', y=column, data=df, alpha=0.8, size=4, color=color, order=order, palette=palette, hue=hue, overflow="shrink") #marker="$\circ$"

    plt.title(f"{title}", fontsize=16, fontweight='bold')
    #plt.title(f"p-value = {pval:.4e}", fontsize=11, color='gray')
    plt.xlabel("Genotype", fontsize=12, fontweight='bold')
    plt.xticks(fontweight='bold', fontsize=12, fontstyle='italic', rotation=45, ha='right') #set the x axis font size
    plt.ylabel(y_label, fontsize=12, fontweight='bold')
    #plt.gca().invert_xaxis() #reverse the order of the x axis
    plt.yticks(fontsize=12, fontweight='bold') #set the y axis font size
    plt.subplots_adjust(left=0.15) #add extra padding to the left of the plot
    plt.locator_params(axis='y', nbins=10) #set the number of ticks on the y axis
    plt.gca().get_yaxis().set_major_formatter(FuncFormatter(lambda x, loc: "{:,}".format(int(x)) if x >= 0 else ""))
    
    # OPTIONAL: enable log scale when caller passes yscale="log" in kwargs
    if kwargs.get("yscale", "linear") == "log":
        ax = plt.gca()
        # replace non-positive bottom with a small positive value to avoid errors
        ymin, ymax = ax.get_ylim()
        if ymin <= 0:
            # choose a sensible small positive lower bound (adjust if needed)
            ymin = 1e-3
            ax.set_ylim(bottom=ymin, top=ymax)
        ax.set_yscale("log")
        # adjust y-axis formatter for logs (optional)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{y:g}"))
        
    #get the span of the y axis, and add a 3% padding to the bottom
    ymin, ymax = plt.gca().get_ylim()
    y_span = ymax - ymin
    y_span_percent_pad_bottom = y_span * 0.03
    plt.gca().set_ylim(bottom=0 - y_span_percent_pad_bottom)

    y_span_percent_pad_top = y_span * 0.08
    plt.gca().set_ylim(top=ymax + y_span_percent_pad_top)

    if show_medians:
        add_median_labels(plt.gca(), df, column)

    if comparisons:
        n_plotted = len(plt.gca().get_xticklabels())
        add_significance_brackets(plt.gca(), comparisons, rename_dict=rename_dict, show_ns=show_ns,
                                  stars_only=n_plotted > max_genotypes_for_brackets)

    plt.tight_layout()

    if save:
        plt.savefig(save_name, dpi=500, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.close()

def plot_plotly(
        df: pd.DataFrame, 
        column: str, 
        title: str, 
        genotype1: str, 
        genotype2: str, 
        by: str = "Genotype", 
        root_fig_name: str = "premeiotic_final.html", 
        save:bool = True
        ) -> None:
    """
    Creates a violin plot of the given column in the given DataFrame, grouped by genotype, using Plotly.

    Parameters:
        df (pandas.DataFrame): The DataFrame containing the data to plot.
        column (str): The name of the column to plot.
        title (str): The title to use for the plot.
        genotype1 (str): The name of the first genotype to compare.
        genotype2 (str): The name of the second genotype to compare.
        root_fig_name (str): The name to use when saving the plot. Defaults to 'premeiotic_final.png'.
        by (str): The column to group the plot by. Defaults to 'Genotype'.
        save (bool): Whether or not to save the plot. Defaults to True.

    Returns:
        None
    """
    #change ∆ to &#916; (html char) in genotype column, otherwise plotly will not save the figure
    df = df.copy()
    df['genotype'] = df['genotype'].str.replace('∆', '&#916;')

    fig = px.violin(df, x="genotype", y=column, box=True, points="all", hover_data={'Sample': True}, color='genotype')
    fig.update_layout(title=f"{title} by {by}",
                      xaxis_title=by,
                      yaxis_title=f"{title}",
                      font=dict(
                        #family="Courier New, monospace",
                        size=18,
                        color="#7f7f7f")
                      )
    fig.update_traces(meanline_visible=True,
                      scalemode='count')
    
    if save:
        name=f"{genotype1}_vs_{genotype2}_{column}_{root_fig_name}"
        fig.write_html(f"figures/clusters/plotly/{name}", auto_open=False)

def convert_columns_to_correct_type(df: pd.DataFrame) -> pd.DataFrame:
    """
    Converts the columns of the DataFrame to the correct data type.

    Parameters:
        

    Returns:
        pandas.DataFrame: The DataFrame with the columns converted to the correct data type.
    """
    df['Avg_Cluster_Length'] = df['Avg_Cluster_Length'].astype(float)
    df['Avg_Cluster_Size_Mutations'] = df['Avg_Cluster_Size_Mutations'].astype(float)
    df['mutations'] = df['mutations'].astype(int)
    df['clusters'] = df['clusters'].astype(int)
    return df

def groupby_genotypes_and_print_mean_std(
        df: pd.DataFrame, 
        set1: pd.DataFrame, 
        set2: pd.DataFrame, 
        column: str, 
        genotype1: str, 
        genotype2: str
        ) -> tuple:
    """
    Groups by genotypes and prints mean and standard deviation, and conducts t-test and Mann-Whitney U test.

    Parameters:
        df (pandas.DataFrame): The DataFrame to group by.
        set1 (pandas.DataFrame): The first DataFrame to use in the test.
        set2 (pandas.DataFrame): The second DataFrame to use in the test.
        column (str): The name of the column to use in the test.
    """
    print(f"\n***COMPARISON OF: {column}***")
    print(df.groupby(['genotype'])[column].agg(['mean', 'median', 'std', 'min', 'max']))
    pval = conduct_ttest(set1, set2, column)
    pval_u = conduct_mannutest(set1, set2, column)
    #pval_f = conduct_f_oneway(set1, set2, column)
    try:
        pval_k = conduct_kurtosis_test(set1, set2, column, genotype1, genotype2)
    except Exception as e:
        print(e)
        pval_k = None
    #conduct_bartlet_test(set1, set2, column)
    conduct_levene_test(set1, set2, column)
    conduct_2_sample_kolmogorov_smirnov_test(set1, set2, column)
    print()
    return pval, pval_u

def  conduct_cluster_data_analysis(
        cluster_caller: str,
        snp_num: int,
        sig: float | str,
        association: bool | None,
        association_feature: Literal["GC", "CO", "both"],
        ploidy: Literal["all", "haploid", "aneuploid"],
        genotype1: str, 
        genotype2: str, 
        save_plt: bool = True, 
        save_plotly: bool = True,
        genotypes_include: None | list = None,
        remove_RTG: bool = True,
        **kwargs: dict
        ) -> None:
    
    from mayo.settings import haploid_samples
    from mayo.settings import rtg_list
    from mayo.settings import config as cfg
    from mayo.settings import genotype_dict as genotype_dict_master
    genotype_dict = genotype_dict_master[cfg["genotype_dict"]]
    
    if association == True:
        association_suffix = '_associated'
    elif association == False:
        association_suffix = '_not_associated'
    else:
        association_suffix = ''

    filename = f"outputs/clusters/{cluster_caller}/{snp_num}snp/all_{cluster_caller}_cluster_data_pooled_master_{snp_num}snp_{sig}_{association_feature}_summary{association_suffix}.csv"
    root_fig_name = f'{cluster_caller}_clusters_{sig}_{snp_num}snp_ploidy_{ploidy}_{association_feature}{association_suffix}'

    df = pd.read_csv(f"{filename}")
    print(f"Conducting cluster data analysis for {filename}")
    df = convert_columns_to_correct_type(df)

    #fill in missing samples
    
    df = fill_in_missing_samples(df, genotype_dict)

    if remove_RTG:
        df = df[~df['Sample'].isin(rtg_list)].copy()

    #make sure samples are/ are not haploid if specified
    if ploidy == "all":
        logging.info("Using all samples")
    elif ploidy == "haploid":
        logging.info("Using haploid samples only")
        df = df[df['Sample'].isin(haploid_samples)].copy()
    elif ploidy == "aneuploid":
        logging.info("Using aneuploid samples only")
        df = df[~df['Sample'].isin(haploid_samples)].copy()

    #filter by genotype
    genotype_list = [genotype1, genotype2]
    if genotypes_include:
        genotype_list = genotypes_include
    df = df[df['genotype'].isin(genotype_list)].copy()
    df["clustered_mutations"] = df["clusters"]*df["Avg_Cluster_Size_Mutations"]

    #for the pairplot
    columns_to_plot = [
        'Avg_Cluster_Length',
        'Avg_Cluster_Size_Mutations',
        'mutations',
        'scattered_mutations',
        'clusters',
        'total_ssDNA_kb'
        ]
    
    if association != None:
        columns_to_plot.remove('mutations')

    plt.show() #this is needed to show the plot, otherwise it will close before it is shown, not sure why
    #plot_pairplot_plt(df, columns_to_plot, save_name=f'{genotype1}_vs_{genotype2}_pairplot.png')

    #create a dictionary of the column names and their titles to be used in the plot (violins)
    columns_to_plot_titles_dict = {
        'Avg_Cluster_Length': 'Average Cluster Length',
        'Avg_Cluster_Size_Mutations': 'Average Cluster Mutations',
        'mutations': 'Mutations per Genome',
        'scattered_mutations': 'Scattered Mutations per Genome',
        'clusters': 'Clusters per Genome',
        'total_ssDNA_kb': 'Mutagenic ssDNA per Genome'}
    
    columns_to_plot_y_labels_dict = {
        'Avg_Cluster_Length': 'Average Cluster Length (bp)',
        'Avg_Cluster_Size_Mutations': 'Average Cluster Mutations',
        'mutations': 'Number of Mutations',
        'scattered_mutations': 'Number of Scattered Mutations',
        'clusters': 'Number of Clusters',
        'total_ssDNA_kb': 'Total ssDNA (kb)'}
    
    # genotype1 is the REFERENCE genotype: every other genotype is compared against it.
    if genotype1 not in set(df['genotype']):
        raise ValueError(
            f"Reference genotype '{genotype1}' has no samples in the data being analyzed "
            f"(genotypes present: {sorted(set(df['genotype']))}). Add it to the genotype list or pick another reference.")
    if not genotype2 or genotype2 not in set(df['genotype']):
        # genotype2 is only used for the first, detailed comparison and file names; default to the first other genotype
        others = [g for g in (genotypes_include or sorted(set(df['genotype']))) if g != genotype1 and g in set(df['genotype'])]
        genotype2 = others[0] if others else genotype1

    set1 = df[df['genotype'] == genotype1]
    set2 = df[df['genotype'] == genotype2]

    print(f"Reference genotype: {genotype1}")
    print(f"{genotype1} = {len(set1)}")
    print(f"{genotype2} = {len(set2)}")

    all_genotypes = [g for g in (genotypes_include or sorted(set(df['genotype']))) if g in set(df['genotype']) and g != genotype1]
    stats_rows = []  # every reference-vs-genotype Mann-Whitney p-value, saved to CSV at the end

    for column in columns_to_plot:

        print(df.groupby(['genotype'])[column].count())

        #groupby genotypes and print mean and standard deviation, and conduct t-test, and mann-whitney u test
        pval, pval_u = groupby_genotypes_and_print_mean_std(df, set1, set2, column, genotype1, genotype2)

        for genotype in all_genotypes:

            print(f"\n###{genotype1} vs {genotype}###")
            set1_secondary = df[df['genotype'] == genotype1]
            set2_secondary = df[df['genotype'] == genotype]
            row = {"metric": column, "reference": genotype1, "genotype": genotype,
                   "n_reference": len(set1_secondary), "n_genotype": len(set2_secondary),
                   "median_reference": set1_secondary[column].median(), "median_genotype": set2_secondary[column].median(),
                   "mean_reference": set1_secondary[column].mean(), "mean_genotype": set2_secondary[column].mean()}
            try:

                print(f"ALL SAMPLES:")
                row["p_mannwhitney_all"] = conduct_mannutest(set1_secondary, set2_secondary, column)
            except Exception as e:
                print("Could not conduct Mann-Whitney U test on secondary comparison, skipping...")

            try:
                print(f"NO 20% OUTLIERS:")

                set1_secondary_no_outliers = set1_secondary[set1_secondary[column] > set1_secondary[column].quantile(0.2)]
                set2_secondary_no_outliers = set2_secondary[set2_secondary[column] > set2_secondary[column].quantile(0.2)]

                row["p_mannwhitney_no_lowest_20pct"] = conduct_mannutest(set1_secondary_no_outliers, set2_secondary_no_outliers, column)
            except Exception as e:
                print("Could not conduct Mann-Whitney U test on secondary comparison: no outliers, skipping...")


            try:
                print(f"NO UNINDUCED:")

                set_1_secondary_uninduced = set1_secondary[set1_secondary['mutations'] > 10]
                set_2_secondary_uninduced = set2_secondary[set2_secondary['mutations'] > 10]

                row["p_mannwhitney_induced_only"] = conduct_mannutest(set_1_secondary_uninduced, set_2_secondary_uninduced, column)
            except Exception as e:
                print("Could not conduct Mann-Whitney U test on secondary comparison: uninduced, skipping...")
            stats_rows.append(row)

        # brackets/asterisks on the violin plot use the all-samples Mann-Whitney U p-values
        comparisons = [(r["reference"], r["genotype"], r["p_mannwhitney_all"])
                       for r in stats_rows if r["metric"] == column and "p_mannwhitney_all" in r]

        plot_violin_plt(
            df=df.copy(),
            column=column,
            save_name=f"figures/clusters/{genotype1}_vs_{genotype2}_{column}_{root_fig_name}.png",
            title=columns_to_plot_titles_dict[column],
            y_label=columns_to_plot_y_labels_dict[column],
            save=save_plt,
            comparisons=comparisons,
            **kwargs
            )
        
        plot_plotly(
            df, 
            column, 
            columns_to_plot_titles_dict[column], 
            genotype1, 
            genotype2, 
            by="Genotype", 
            root_fig_name=f'{root_fig_name}.html', 
            save=save_plotly)
        

    # save every reference-vs-genotype comparison in one table
    stats_path = f"outputs/clusters/stats_vs_{genotype1}_{root_fig_name}.csv"
    os.makedirs(os.path.dirname(stats_path), exist_ok=True)
    pd.DataFrame(stats_rows).to_csv(stats_path, index=False, encoding="utf-8")
    print(f"\nSaved all {genotype1}-vs-genotype Mann-Whitney U p-values to {stats_path}")
    logging.info(f"Saved reference ({genotype1}) comparisons to {stats_path}")
