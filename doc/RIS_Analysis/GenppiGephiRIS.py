#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, kendalltau, wilcoxon, ks_2samp
from itertools import combinations
from sklearn.metrics import jaccard_score
from matplotlib.backends.backend_pdf import PdfPages
from colorama import Fore, Style
from scipy.stats import kendalltau, ks_2samp
from collections import defaultdict

MGP = {}

# Function to verify the number of files in a directory
def count_files_in_directory(DIRECTORY):    
    items = os.listdir(DIRECTORY)
    files = [item for item in items if os.path.isfile(os.path.join(DIRECTORY, item))]
    return len(files)

# Function to load CSV files
def load_csv_files(DIRECTORY, base_name, num_files=None):
    files = sorted([f for f in os.listdir(DIRECTORY) if f.startswith(base_name) and f.endswith(".csv")])
    if num_files is not None:
        files = files[:num_files]
    data = {}
    for file in files:
        path = os.path.join(DIRECTORY, file)
        df = pd.read_csv(path)
        expected_columns = {'degree', 'betweenesscentrality', 'bridgingcentrality', 'Label'}
        existing_columns = set(df.columns)
        if not expected_columns.issubset(existing_columns):
            print(f"Error in file {file}: Expected columns not found.")
            print(f"Existing columns: {existing_columns}")
            continue
        data[file] = df
    return data

# Function to select a metric
def select_metric(option):
    metrics = {1: 'degree', 2: 'betweenesscentrality', 3: 'bridgingcentrality'}
    return metrics.get(option, 'degree')

# Function to get top nodes per file
def get_top_nodes_per_file(data, TOP_N, metric, base_name):
    top_nodes_per_file = {}
    for file, df in data.items():
        file_number = file.split(base_name)[1].split(".csv")[0]
        top_n_nodes = df.nlargest(TOP_N, metric)[['Label', metric]].set_index('Label').to_dict()[metric]
        top_nodes_per_file[file_number] = top_n_nodes
    return top_nodes_per_file


# In[2]:


# Function to return a formatted string with top nodes and their values
def print_top_nodes_with_values(top_nodes_per_file, metric, metric_name, TOP_P):
    # Setting larger spacing for columns
    label_column = 30  # Space for "Label" column
    metric_column = 30  # Space for "Metric" column
    file_column = 20    # Space for "File" column
    
    # Header formatting
    result = "{:>{}} {:>{}} {:>{}}\n".format(
        "Label", 30,
        "File", 25,
        metric_name, 30
    )
    result += "             " + "-" * (label_column + metric_column + file_column + 2) + "\n"
    
    for i, (file_number, node_values) in enumerate(top_nodes_per_file.items()):
        if i >= 10:
            result += "\n   ...          continues    ...      (formatted sample for printing)\n" 
            break    
      
        # Sort nodes by metric value in descending order
        sorted_nodes = sorted(node_values.items(), key=lambda x: x[1], reverse=True)           
        # Select only Top_P nodes
        top_p_nodes = sorted_nodes[:TOP_P]        
        if metric_name == "Degree":
            for node, value in top_p_nodes:
                result += "{:>{}} {:>{}} {:>{}}\n".format(
                    node, label_column,                    
                    file_number, file_column,
                    value, metric_column
                )
        else:
            for node, value in top_p_nodes:
                result += "{:>{}} {:>{}}{:>{}.4f}\n".format(
                    node, label_column,
                    file_number, file_column,
                    value, metric_column
                )
    
    return result


# In[3]:


# Function to return a string with top_nodes with presence percentage
def print_top_nodes_with_percentage(top_nodes_per_file, total_files, visual=0):
    result = {}
    for file_number, top_n_nodes in top_nodes_per_file.items():
        for node in top_n_nodes:
            if node not in result:
                result[node] = [file_number]
            else:
                result[node].append(file_number)
    
    result_str = "\nTopN Nodes\n\n{:<15} {:<15} {:>15}\n".format("% Presence", "Label", "Files")
    result_str += "-" * 45 + "\n"
    
    # Dynamic calculation of line limit
    if visual == 0:
        # Accesses global variables (assuming they exist in scope)
        product = TOP_N * num_files_to_print
        adjustment = round((product) ** (1/3))  # Cubic root
        max_lines = max(5, 30 - adjustment)     # Minimum of 5 lines # Indicated of 30 lines
    else:
        max_lines = float('inf')  # Show all lines

    # Generates formatted node list
    for i, (node, files) in enumerate(result.items()):
        if visual == 0 and i >= max_lines:
            result_str += "\n   ...          continues    ...      (formatted sample for printing)\n" 
            break
            
        percentage = (len(files) / total_files) * 100        
        formatted_percentage = f"{percentage:06.2f} %"
        result_str += "{:<15}{:<15}{:<15}\n".format(formatted_percentage, node, ', '.join(files))
        
     # Weighted average calculation
    weighted_sum = 0.0
    total_weights = 0
    
    for node, files in result.items():
        weight = len(files)
        percentage = (weight / total_files) * 100
        weighted_sum += percentage * weight
        total_weights += weight
    
    global_average = (weighted_sum / total_weights) if total_weights != 0 else 0.0
    result_str += f"\nGlobal node Presence Mean (Weighted): {global_average:.2f}%\n"
    
    return result_str


# In[4]:


def calculate_presence_percentage(top_nodes_per_file, total_files):
    result = {}
    for file_number, top_n_nodes in top_nodes_per_file.items():
        for node in top_n_nodes:
            if node not in result:
                result[node] = [file_number]                
            else:
                result[node].append(file_number)   
       
     # Weighted average calculation
    weighted_sum = 0.0
    total_weights = 0
    
    for node, files in result.items():
        weight = len(files)
        percentage = (weight / total_files) * 100
        weighted_sum += percentage * weight
        total_weights += weight
    
    global_average = (weighted_sum / total_weights) if total_weights != 0 else 0.0
      
    return global_average


# In[5]:


from collections import defaultdict

def calculate_global_metrics(top_nodes_per_file, original_files, global_metrics_df):
    results = []
    files = list(top_nodes_per_file.keys())
    ks_stats = []  # To store KS statistics for all pairs

    # Iterate over all combinations of file pairs
    for (file1, file2) in combinations(files, 2):
        set1 = set(top_nodes_per_file[file1].keys())
        set2 = set(top_nodes_per_file[file2].keys())
        
        # Calculate Jaccard coefficient
        jaccard = len(set1 & set2) / len(set1 | set2) if len(set1 | set2) > 0 else 0
        
        # Calculate overlap
        overlap = len(set1 & set2) / min(len(set1), len(set2)) if min(len(set1), len(set2)) > 0 else 0

        # Calculate Kolmogorov-Smirnov test (KS)
        values_file1 = list(top_nodes_per_file[file1].values())
        values_file2 = list(top_nodes_per_file[file2].values())
        
        if len(values_file1) > 1 and len(values_file2) > 1:
            ks_stat, ks_p = ks_2samp(values_file1, values_file2)
            ks_stats.append(ks_stat)  # Store KS statistic
        else:
            ks_stat, ks_p = np.nan, np.nan
        
        # Calculate Kendall correlation
        nodes_file1 = top_nodes_per_file[file1]
        nodes_file2 = top_nodes_per_file[file2]
        common_nodes = set(nodes_file1.keys()) & set(nodes_file2.keys())
        
        if len(common_nodes) >= 2:
            values1 = [nodes_file1[node] for node in common_nodes]
            values2 = [nodes_file2[node] for node in common_nodes]
            kendall_corr, _ = kendalltau(values1, values2)
        else:
            kendall_corr = np.nan
        
        # Map file numbers to full names
        name_file1 = original_files[int(file1) - 1] if int(file1) <= len(original_files) else f"File {file1}"
        name_file2 = original_files[int(file2) - 1] if int(file2) <= len(original_files) else f"File {file2}"
        
        # Store results
        results.append({
            'File1': name_file1, 
            'File2': name_file2,
            'Jaccard': jaccard, 
            'Overlap': overlap,
            'KS_p': ks_p,
            'KS_stat': ks_stat,  # Adds KS statistic
            'Kendall': kendall_corr
        })
    
    # Convert results to DataFrame and store globally
    results_df = pd.DataFrame(results)
    global_metrics_df = results_df  # Updates global variable
    
    # Calculate global metrics
    # 1. Mean Jaccard
    mean_jaccard = results_df['Jaccard'].mean()
    
    # Calculate Fleiss' Kappa
    fleiss_kappa = calculate_fleiss_kappa(top_nodes_per_file)
    
    # 3. Mean KS Statistics (D̄)
    mean_ks_stat = np.nanmean(ks_stats) if ks_stats else np.nan
    
    # 4. KS Test for Multiple Samples (D_global)
    # Calculate global cumulative distribution
    all_values = np.concatenate([list(top_nodes_per_file[file].values()) for file in files])
    global_cdf = np.sort(all_values)

    # 5. Global Metrics for Kendall Tau
    mean_kendall = results_df['Kendall'].mean()
    median_kendall = results_df['Kendall'].median()
    positive_kendall_ratio = (results_df['Kendall'] > 0).mean()
    
    # Calculate D_global for each file
    d_global_values = []
    for file in files:
        file_values = np.sort(list(top_nodes_per_file[file].values()))
        ecdf_file = np.searchsorted(file_values, global_cdf, side='right') / len(file_values)
        ecdf_global = np.searchsorted(global_cdf, global_cdf, side='right') / len(global_cdf)
        d_global = np.max(np.abs(ecdf_file - ecdf_global))
        d_global_values.append(d_global)
    
    d_global_mean = np.nanmean(d_global_values) if d_global_values else np.nan

    p_values = []
    for (file1, file2) in combinations(files, 2):
        values_file1 = list(top_nodes_per_file[file1].values())
        values_file2 = list(top_nodes_per_file[file2].values())
        _, p = ks_2samp(values_file1, values_file2)
        p_values.append(p)

    ks_p_multiple = np.nanmean(p_values) if p_values else np.nan

    ks_multiple_d_mean, ks_multiple_p_mean = calculate_multiple_sample_ks(top_nodes_per_file)
        
    # Create dictionary to store general values
    general_values = {
        'mean_jaccard': mean_jaccard,
        'fleiss_kappa': fleiss_kappa,
        'mean_ks_stat': mean_ks_stat,
        'mean_kendall': mean_kendall,
        'median_kendall': median_kendall,
        'positive_kendall_ratio': positive_kendall_ratio,
        'd_global_mean': d_global_mean,
        'ks_p_multiple': ks_p_multiple,
        'ks_multiple_d_mean': ks_multiple_d_mean,
        'ks_multiple_p_mean': ks_multiple_p_mean
    }
    
    return general_values


# In[6]:


# Global variable to store the DataFrame
global_metrics_df = None

def calculate_metrics(top_nodes_per_file, original_files):
    global global_metrics_df  # Access global variable
    results = []
    files = list(top_nodes_per_file.keys())
    ks_stats = []  # To store KS statistics for all pairs
    
    # Iterate over all combinations of file pairs
    for (file1, file2) in combinations(files, 2):
        set1 = set(top_nodes_per_file[file1].keys())
        set2 = set(top_nodes_per_file[file2].keys())
        
        # Calculate Jaccard coefficient
        jaccard = len(set1 & set2) / len(set1 | set2) if len(set1 | set2) > 0 else 0
        
        # Calculate overlap
        overlap = len(set1 & set2) / min(len(set1), len(set2)) if min(len(set1), len(set2)) > 0 else 0

        # Calculate Kolmogorov-Smirnov test (KS)
        values_file1 = list(top_nodes_per_file[file1].values())
        values_file2 = list(top_nodes_per_file[file2].values())
        
        if len(values_file1) > 1 and len(values_file2) > 1:
            ks_stat, ks_p = ks_2samp(values_file1, values_file2)
            ks_stats.append(ks_stat)  # Store KS statistic
        else:
            ks_stat, ks_p = np.nan, np.nan
        
        # Calculate Kendall correlation
        nodes_file1 = top_nodes_per_file[file1]
        nodes_file2 = top_nodes_per_file[file2]
        common_nodes = set(nodes_file1.keys()) & set(nodes_file2.keys())
        
        if len(common_nodes) >= 2:
            values1 = [nodes_file1[node] for node in common_nodes]
            values2 = [nodes_file2[node] for node in common_nodes]
            kendall_corr, _ = kendalltau(values1, values2)
        else:
            kendall_corr = np.nan
        
        # Map file numbers to full names
        name_file1 = original_files[int(file1) - 1] if int(file1) <= len(original_files) else f"File {file1}"
        name_file2 = original_files[int(file2) - 1] if int(file2) <= len(original_files) else f"File {file2}"
        
        # Store results
        results.append({
            'File1': name_file1, 
            'File2': name_file2,
            'Jaccard': jaccard, 
            'Overlap': overlap,
            'KS_p': ks_p,
            'KS_stat': ks_stat,
            'Kendall': kendall_corr
        })
    
    # Convert results to DataFrame and store globally
    results_df = pd.DataFrame(results)
    global_metrics_df = results_df  # Update global variable
    
    # Calculate global metrics
    mean_jaccard = results_df['Jaccard'].mean()
    fleiss_kappa = calculate_fleiss_kappa(top_nodes_per_file)
    mean_ks_stat = np.nanmean(ks_stats) if ks_stats else np.nan
    
    # Calculate global cumulative distribution
    all_values = np.concatenate([list(top_nodes_per_file[file].values()) for file in files])
    global_cdf = np.sort(all_values)

    # Kendall Tau global metrics
    mean_kendall = results_df['Kendall'].mean()
    median_kendall = results_df['Kendall'].median()
    positive_kendall_ratio = (results_df['Kendall'] > 0).mean()

    # Multiple sample KS test (D_global)
    d_global_values = []
    for file in files:
        file_values = np.sort(list(top_nodes_per_file[file].values()))
        ecdf_file = np.searchsorted(file_values, global_cdf, side='right') / len(file_values)
        ecdf_global = np.searchsorted(global_cdf, global_cdf, side='right') / len(global_cdf)
        d_global = np.max(np.abs(ecdf_file - ecdf_global))
        d_global_values.append(d_global)
    
    d_global_mean = np.nanmean(d_global_values) if d_global_values else np.nan

    p_values = []
    for (file1, file2) in combinations(files, 2):
        values_file1 = list(top_nodes_per_file[file1].values())
        values_file2 = list(top_nodes_per_file[file2].values())
        _, p = ks_2samp(values_file1, values_file2)
        p_values.append(p)

    ks_p_multiple = np.nanmean(p_values) if p_values else np.nan
    ks_multiple_d_mean, ks_multiple_p_mean = calculate_multiple_sample_ks(top_nodes_per_file)
    
    # Format results table
    table_str = "Application of Jaccard and Overlap Coefficients, Kolmogorov-Smirnov Test (KSp)\nand Kendall Tau Correlation\n\n"
    table_str += "{:<25} {:<24} {:<15} {:<15} {:<19} {:<15}\n".format(
        "File1", "File2", "Jaccard", "Overlap", "KS_p", "Kendall"
    )
    table_str += "-" * 120 + "\n"
    for i, result in enumerate(results):
        if i >= 50: # Suggested number: 50
            table_str += "\n   ...          continues    ...      (formatted sample for printing)\n" 
            break
        table_str += "{:<20} {:<20} {:<15.4f} {:<15.4f} {:<15.4f} {:>10.4f}\n".format(
            result['File1'], 
            result['File2'], 
            result['Jaccard'], 
            result['Overlap'],
            result['KS_p'] if not np.isnan(result['KS_p']) else np.nan,
            result['Kendall'] if not np.isnan(result['Kendall']) else np.nan
        )
    
    # Format global metrics
    fleiss_kappa_str = f"{fleiss_kappa:.4f}" if not np.isnan(fleiss_kappa) else "N/A"
    mean_ks_stat_str = f"{mean_ks_stat:.4f}" if not np.isnan(mean_ks_stat) else "N/A"
    d_global_mean_str = f"{d_global_mean:.4f}" if not np.isnan(d_global_mean) else "N/A"
    ks_p_multiple_str = f"{ks_p_multiple:.4f}" if not np.isnan(ks_p_multiple) else "N/A"
    ks_multiple_d_str = f"{ks_multiple_d_mean:.4f}" if not np.isnan(ks_multiple_d_mean) else "N/A"
    ks_multiple_p_str = f"{ks_multiple_p_mean:.4f}" if not np.isnan(ks_multiple_p_mean) else "N/A"
    mean_kendall_str = f"{mean_kendall:.4f}" if not np.isnan(mean_kendall) else "N/A"
    median_kendall_str = f"{median_kendall:.4f}" if not np.isnan(median_kendall) else "N/A"
    positive_kendall_str = f"{positive_kendall_ratio:.2%}" if not np.isnan(positive_kendall_ratio) else "N/A"

    # Add global metrics to table
    table_str += "\n\nGlobal Metrics:\n"
    table_str += f"Mean Jaccard Coefficient (J̄): {mean_jaccard:.4f}\n"
    table_str += f"Fleiss' Kappa Agreement Index (κF): {fleiss_kappa_str}\n"
    table_str += f"Mean KS Distance Between Pairs (D̄): {mean_ks_stat_str}\n"
    table_str += f"Mean p-value for KS Test Pairs: {ks_p_multiple_str}\n"
    table_str += f"Mean KS Distance for Multiple Samples (D̄_mult): {ks_multiple_d_str}\n"
    table_str += f"Mean p-value for Multiple Sample KS Test (p̄_mult): {ks_multiple_p_str}\n"
    table_str += f"Mean Kendall Tau (τ̄): {mean_kendall_str}\n"
    table_str += f"Median Kendall Tau (τ̃): {median_kendall_str}\n"
    table_str += f"Percentage of Pairs with τ > 0: {positive_kendall_str}"
    
    return table_str


# In[7]:


import numpy as np

def calculate_fleiss_kappa(top_nodes_per_file):
    """
    Calculates Fleiss' Kappa for agreement between files about the presence of top nodes.
    
    Args:
        top_nodes_per_file (dict): Dictionary where keys are file identifiers and values
                                 are dictionaries of top nodes and their metrics
    
    Returns:
        float: Fleiss' Kappa value between -1 and 1
    
    Raises:
        ValueError: If input is invalid
    """
    # Input validation
    if not top_nodes_per_file or not isinstance(top_nodes_per_file, dict):
        raise ValueError("Invalid input: 'top_nodes_per_file' must be a non-empty dictionary.")
    
    # 1. List all unique nodes
    all_nodes = set()
    for file_data in top_nodes_per_file.values():
        if not isinstance(file_data, dict):
            raise ValueError("Values in 'top_nodes_per_file' must be dictionaries.")
        all_nodes.update(file_data.keys())
    all_nodes = list(all_nodes)
    
    # 2. Rating matrix: each row is a node, each column is a file
    k = len(top_nodes_per_file)  # Number of raters (files)
    N = len(all_nodes)  # Number of items (nodes)
    rating_matrix = np.zeros((N, k), dtype=int)
    
    for i, node in enumerate(all_nodes):
        for j, file_data in enumerate(top_nodes_per_file.values()):
            rating_matrix[i, j] = 1 if node in file_data else 0
    
    # 3. Count ratings per category (0 or 1) for each item
    n_ik = np.zeros((N, 2), dtype=int)
    for i in range(N):
        counts = np.bincount(rating_matrix[i], minlength=2)
        n_ik[i, :] = counts[:2]
    
    # 4. Observed agreement (P_e)
    P_e = 0.0
    for i in range(N):
        sum_nik2 = np.sum(n_ik[i] * (n_ik[i] - 1))  # Corrected: Subtracts n_ij to avoid duplicates
        P_e += sum_nik2 / (k * (k - 1))  # Normalizes by rater pairs
    P_e /= N
    
    # 5. Expected agreement (P_c)
    p_k = np.sum(n_ik, axis=0) / (N * k)  # Proportion of ratings in each category
    P_c = np.sum(p_k ** 2)
    
    # 6. Fleiss Kappa
    if (1 - P_c) != 0:
        kappa = (P_e - P_c) / (1 - P_c)
    else:
        kappa = 1.0  # Perfect agreement
    
    # Ensure result is in [-1, 1] range
    kappa = max(-1.0, min(1.0, kappa))
    
    return kappa


# In[8]:


from scipy.stats import ks_2samp
import numpy as np

def calculate_multiple_sample_ks(top_nodes_per_file):
    """
    Calculates the Kolmogorov-Smirnov test for multiple samples.
    Returns the mean D statistics and mean p-value.
    
    Args:
        top_nodes_per_file (dict): Dictionary where keys are file identifiers and values
                                 are dictionaries of top nodes and their metrics
    
    Returns:
        tuple: (mean_ks_statistic, mean_p_value)
    """
    files = list(top_nodes_per_file.keys())
    
    # 1. Calculate global distribution
    all_values = np.concatenate([list(top_nodes_per_file[file].values()) for file in files])
    global_cdf = np.sort(all_values)
    
    # 2. Compare each sample with global distribution
    d_values = []
    p_values = []
    for file in files:
        file_values = np.sort(list(top_nodes_per_file[file].values()))
        d_stat, p_val = ks_2samp(file_values, global_cdf)
        d_values.append(d_stat)
        p_values.append(p_val)
    
    # 3. Aggregate results
    mean_ks_statistic = np.nanmean(d_values) if d_values else np.nan
    mean_p_value = np.nanmean(p_values) if p_values else np.nan
    
    return mean_ks_statistic, mean_p_value


# In[9]:


def add_page(header, pdf_pages, content, title, paginate=False):
    """
    Adds a page to the PDF. If paginate=True, splits content across multiple pages.
    """
    lines_per_page = 40  # Number of lines per page (adjust as needed)
    lines = content.split('\n')  # Split content into lines
    
    if paginate and len(lines) > lines_per_page:
        # Split content into smaller chunks
        for i in range(0, len(lines), lines_per_page):
            chunk = '\n'.join(lines[i:i + lines_per_page])
            _add_single_page(header, pdf_pages, chunk, title if i == 0 else f"{title} (Continuation)")
    else:
        # Add content in a single page
        _add_single_page(header, pdf_pages, content, title)

def _add_single_page(header, pdf_pages, content, title):
    """
    Helper function to add a single page to the PDF.
    """
    plt.figure(figsize=(8, 11))
    
    # Add header
    plt.figtext(0.5, 0.97, header.split('\n')[1], fontsize=14, wrap=True, va='top', ha='center', weight='bold')  
    plt.figtext(0.1, 0.93, header.split('\n')[2], fontsize=10, wrap=True, va='top', color='red')
    plt.figtext(0.225, 0.93, header.split('\n')[3], fontsize=10, wrap=True, va='top')  
    plt.figtext(0.1, 0.91, header.split('\n')[4], fontsize=10, wrap=True, va='top')
    
    # Add page title
    plt.figtext(0.5, 0.88, title, fontsize=12, wrap=True, weight='bold', va='top', ha='center', color='green')
    
    # Add content
    plt.figtext(0.1, 0.83, content, fontsize=10, wrap=True, va='top')
    
    # Turn off axes
    plt.axis('off')
    
    # Save page to PDF
    pdf_pages.savefig()
    plt.close()  

def add_single_page_with_table(header, pdf_pages, title, table):
    """
    Helper function to add a single page with a table to the PDF.
    
    :param pdf_pages: PdfPages object to save the page
    :param title: Page title
    :param table: Dictionary containing table data (cellText, colLabels, rowLabels, etc.)
    :param header: Header text to display on the page
    """
    plt.figure(figsize=(8, 11))
    
    # Add header
    plt.figtext(0.5, 0.97, header.split('\n')[1], fontsize=14, wrap=True, va='top', ha='center', weight='bold')
    plt.figtext(0.1, 0.93, header.split('\n')[2], fontsize=10, wrap=True, va='top', color='red')
    plt.figtext(0.225, 0.93, header.split('\n')[3], fontsize=10, wrap=True, va='top')  
    plt.figtext(0.1, 0.91, header.split('\n')[4], fontsize=10, wrap=True, va='top')

    if num_files_to_print < 15:
        digits = 3
        font_size = 9
        label_font_size = 12
    elif num_files_to_print < 25:
        digits = 2
        font_size = 9
        label_font_size = 12
    elif num_files_to_print < 35:
        digits = 2
        font_size = 8
        label_font_size = 12
    else:
        digits = 1
        font_size = 7
        label_font_size = 10
    
    # Add title
    plt.figtext(0.5, 0.82, title, fontsize=12, wrap=True, weight='bold', va='top', ha='center', color='green')
    
    # Format table values with fixed decimal places
    formatted_cell_text = []
    for row in table['cellText']:
        formatted_row = [f"{value:.{digits}f}" if isinstance(value, (int, float)) else value for value in row]
        formatted_cell_text.append(formatted_row)
    
    # Add table
    plt.subplots_adjust(top=0.8)
    table_plot = plt.table(
        cellText=formatted_cell_text,
        colLabels=table['colLabels'],
        rowLabels=table['rowLabels'],
        loc='center',
        cellLoc='center',
        colColours=table['colColours'],
        colWidths=table['colWidths'],
        bbox=[0.0, 0.0, 1.0, 1.0],            
    )
    
    # Set table font sizes
    table_plot.auto_set_font_size(False)
    table_plot.set_fontsize(font_size)

    for key, cell in table_plot.get_celld().items():
        if key[0] == 0:  # Column headers
            cell.set_fontsize(label_font_size)
        if key[1] == -1:  # Row headers
            cell.set_fontsize(label_font_size)
                    
    plt.axis('off')
    pdf_pages.savefig()
    plt.close()

def add_single_page_with_heatmap(header, pdf_pages, title, heatmap_data, min_val, max_val, column):
    """
    Helper function to add a single page with a heatmap to the PDF.
    
    :param pdf_pages: PdfPages object to save the page
    :param title: Page title
    :param heatmap_data: DataFrame containing heatmap data
    :param column: Column name to use for heatmap
    :param header: Header text to display on the page
    """
    plt.figure(figsize=(8, 11))
    
    # Add header
    plt.figtext(0.5, 0.97, header.split('\n')[1], fontsize=14, wrap=True, va='top', ha='center', weight='bold') 
    plt.figtext(0.1, 0.93, header.split('\n')[2], fontsize=10, wrap=True, va='top', color='red')
    plt.figtext(0.225, 0.93, header.split('\n')[3], fontsize=10, wrap=True, va='top')  
    plt.figtext(0.1, 0.91, header.split('\n')[4], fontsize=10, wrap=True, va='top')
    
    # Add title
    plt.figtext(0.5, 0.85, title, fontsize=12, wrap=True, weight='bold', va='top', ha='center', color='green')
    
    # Add heatmap
    plt.subplots_adjust(top=0.82, bottom=0.3)
    pivot_table = heatmap_data.pivot(index='file1', columns='file2', values=column)

    # Set decimal places and font size based on number of files
    if num_files_to_print < 10:
        digits = 2
        font_size = 9
        label_font_size = 10
    elif num_files_to_print < 15:
        digits = 1
        font_size = 9
        label_font_size = 10
    elif num_files_to_print < 20:
        digits = 1
        font_size = 8
        label_font_size = 9
    else:
        digits = 1
        font_size = 7
        label_font_size = 8
    
    ax = plt.axes([0.14, 0.1, 0.85, 0.7])
    show_values = num_files_to_print <= 20
    
    sns.heatmap(
        pivot_table,
        annot=show_values,
        cmap='coolwarm',
        fmt=f'.{digits}f' if show_values else None,
        annot_kws={'size': font_size} if show_values else None,
        vmin=min_val,
        vmax=max_val,
        ax=ax
    )
    
    # Set axis labels
    ax.set_xlabel('File2', color='green', fontsize=label_font_size)
    ax.set_ylabel('File1', color='green', fontsize=label_font_size)
    
    # Set tick colors
    plt.xticks(color='red', fontsize=label_font_size)
    plt.yticks(color='red', fontsize=label_font_size)
    
    pdf_pages.savefig()
    plt.close()


# In[10]:


def append_to_final_pdf(implementation="Temporary.pdf", output_pdf="Combined_Results.pdf"):
    from PyPDF2 import PdfMerger
    import os
    
    # Temporary files
    temp_pdf = "temp_implementation.pdf"   # Optional
    
    # Merge PDF with existing file (if it exists)
    if os.path.exists(output_pdf):
        merger = PdfMerger()
        
        # Add existing file to merger
        with open(output_pdf, "rb") as existing_pdf:
            merger.append(existing_pdf)
        
        # Add new file to merger
        with open(implementation, "rb") as new_pdf:
            merger.append(new_pdf)
        
        # Save result to final file
        with open(output_pdf, "wb") as merged_pdf:
            merger.write(merged_pdf)
        
        merger.close()
        os.remove(implementation)  # Remove temporary file after merging
    else:
        # If file doesn't exist, rename implementation file
        os.rename(implementation, output_pdf)


# In[11]:


def run_calculations(mode, top_nodes, calculation_mode, metric_choice, num_files, implementation_num, output_filename="Test.pdf"):
    global TOP_N, DIRECTORY, base_name, mode_name, actual_metric, metric, header, metric_name, num_files_to_print
    global title, top7, presence_percentage, combined_metrics, combined_heatmaps

    TOP_P = 7
    DIRECTORY = ""
    base_name = ""
    mode_name = ""
    SELECTED_METRIC = metric_choice
    actual_metric = 1
    TOP_N = top_nodes
    metric = 1
    # implementation_num já está definido como parâmetro
    header = ""
    metric_name = ""
    num_files_to_print = num_files
    
    CANCEL = 0
    if calculation_mode == 1:
        DIRECTORY = "/home/william/Documentos/DisserTemp/TestesBA/NCSV"
        base_name = "F"
        mode_name = "Features"
    elif calculation_mode == 2:
        DIRECTORY = "/home/william/Documentos/DisserTemp/TestesBA/MLCSV"
        base_name = "ML"
        mode_name = "Machine Learning"
    else:
        CANCEL = 1
    
    if CANCEL == 1:
        print("Incorrect Parameters: Cannot continue.\n")
        sys.exit("Program terminated due to error.")
    
    if SELECTED_METRIC not in [2, 3]:
        actual_metric = 1
        metric = 1
    else:
        actual_metric = SELECTED_METRIC
        metric = SELECTED_METRIC
    
    actual_num_files = count_files_in_directory(DIRECTORY)
    
    if actual_num_files < num_files_to_print:
        num_files_to_print = actual_num_files
    
    data = load_csv_files(DIRECTORY, base_name, num_files_to_print)
    metric = select_metric(actual_metric)
    metric_name = metric[0].capitalize() + "".join(met.lower() for met in metric[1:])
    
    header = f"""
    Implementation Number {implementation_num}
    Parameters: 
    Top_N = {TOP_N:<{25}} Mode: {mode_name}
    Number of files = {num_files_to_print:<{33}} Selected metric: {metric_name}    
    """
    
    top_nodes_per_file = get_top_nodes_per_file(data, TOP_N, metric, base_name)

    if mode == 3:
        global_avg = calculate_presence_percentage(top_nodes_per_file, len(data))
        general_metrics = calculate_global_metrics(top_nodes_per_file, list(data.keys()), global_metrics_df)
        save_MGP(implementation_num, TOP_N, mode_name, num_files_to_print, global_avg, actual_metric, general_metrics)

    else:
        pdf_pages = PdfPages(output_filename)
        if title:
            print_title_page(pdf_pages, implementation_num)
        if top7:
            print_top7(header, pdf_pages, top_nodes_per_file, metric, metric_name, TOP_P) 
        if presence_percentage:
            print_presence_percentage(header, pdf_pages, data, top_nodes_per_file, 0)
        if combined_metrics:
            print_combined_metrics(header, pdf_pages, top_nodes_per_file, data)
        if combined_heatmaps:      
            print_combined_heatmaps(header, pdf_pages, top_nodes_per_file, data, 0, 1, 1)
            print_combined_heatmaps(header, pdf_pages, top_nodes_per_file, data, 0, 1, 3)
            print_combined_heatmaps(header, pdf_pages, top_nodes_per_file, data, 0, 1, 2)
            print_combined_heatmaps(header, pdf_pages, top_nodes_per_file, data, -1, 1, 4)        
        pdf_pages.close()

    if mode == 2:
        append_to_final_pdf(implementation=output_filename, output_pdf=final_filename)


# In[12]:


def print_top7(header, pdf_pages, top_nodes_per_file, metric, metric_name, TOP_P):
    # Top Nodes with Values
    top_nodes_values_text = print_top_nodes_with_values(top_nodes_per_file, metric, metric_name, TOP_P)
    add_page(header, pdf_pages, top_nodes_values_text, f"Top {TOP_P} Nodes per file", paginate=True)

def print_presence_percentage(header, pdf_pages, data, top_nodes_per_file, view_mode):
    # Top Nodes with Percentage
    if view_mode == 0:
        # Used for sample printing (prints everything on 1 page) (Default configuration)
        top_nodes_percentage_text = print_top_nodes_with_percentage(top_nodes_per_file, len(data))
        add_page(header, pdf_pages, top_nodes_percentage_text, "Top Nodes with Percentage", paginate=True)
    if view_mode == 1:
        # Used for complete printing (need to comment previous commands)
        top_nodes_percentage_text = print_top_nodes_with_percentage(top_nodes_per_file, len(data), 1)
        add_page(header, pdf_pages, top_nodes_percentage_text, "Top Nodes with Percentage", paginate=True)

def print_combined_metrics(header, pdf_pages, top_nodes_per_file, data):
    # Various metrics
    metrics_data = calculate_metrics(top_nodes_per_file, list(data.keys()))
    metrics_text = str(metrics_data)
    add_page(header, pdf_pages, metrics_text, 
             "Similarity, Correlation and Distribution Tests \nbetween files", 
             paginate=True)


# In[13]:


def print_combined_heatmaps(header, pdf_pages, top_nodes_per_file, data, min_val, max_val, desired_calculation):
    """
    Function to generate a heatmap based on the desired calculation.

    :param header: Header text to display on the page
    :param pdf_pages: PdfPages object to save the page
    :param top_nodes_per_file: Dictionary with top nodes per file
    :param data: Dictionary with complete file data
    :param desired_calculation: String or integer defining the calculation to display.
                               Can be "Jaccard", "Overlap", "KS_p", "Kendall" or 1, 2, 3, 4.
    """
    global global_metrics_df  # Accesses the global DataFrame

    # Mapping of calculations to DataFrame columns
    calculation_mapping = {
        "Jaccard": ("Jaccard", r"Heatmap of $Jaccard$ Coefficient"),
        "Overlap": ("Overlap", r"Heatmap of $Overlap$ Coefficient"),
        "KS_p": ("KS_p", r"Heatmap of $KS_p$ Value"),
        "Kendall": ("Kendall", r"Heatmap of $Kendall$ Tau Correlation"),
        1: ("Jaccard", r"Heatmap of Jaccard Coefficient"),
        2: ("Overlap", r"Heatmap of Overlap Coefficient"),
        3: ("KS_p", r"Heatmap of Kolmogorov-Smirnov Test ($KS_p$)"),
        4: ("Kendall", r"Heatmap of Kendall Tau Correlation")
    }

    # Validate the desired calculation
    if desired_calculation not in calculation_mapping:
        raise ValueError(f"Invalid calculation. Options: {list(calculation_mapping.keys())}")

    # Get corresponding column and title
    column, title = calculation_mapping[desired_calculation]

    # Check if global DataFrame exists
    if global_metrics_df is not None:
        # Rename columns for heatmap function compatibility
        metrics_df = global_metrics_df.rename(columns={'File1': 'file1', 'File2': 'file2'})

        # Generate heatmap
        add_single_page_with_heatmap(
            header,
            pdf_pages,
            title,
            metrics_df[['file1', 'file2', column]],
            min_val,
            max_val,
            column,
        )
    else:
        print("Warning: No metric calculations performed. Run calculate_metrics first.")


# In[14]:


def save_MGP(implementation_num, TOP_N, mode_name, num_files_to_print, global_avg, actual_metric, general_metrics):
    """
    Stores the weighted global presence average (MGP) in the MGP dictionary.
    The key is implementation_num, and the value is a dictionary with parameters and separated metrics.
    """
    global MGP  # Indicates we're using the global MGP variable
    
    MGP[implementation_num] = {
        "TOP_N": TOP_N,
        "Mode": mode_name,
        "Metric": actual_metric,
        "num_files": num_files_to_print,
        "MGP": global_avg,
        "mean_jaccard": general_metrics['mean_jaccard'],
        "fleiss_kappa": general_metrics['fleiss_kappa'],
        "mean_ks_stat": general_metrics['mean_ks_stat'],
        "mean_kendall": general_metrics['mean_kendall'],
        "median_kendall": general_metrics['median_kendall'],
        "positive_kendall_ratio": general_metrics['positive_kendall_ratio'],
        "d_global_mean": general_metrics['d_global_mean'],
        "ks_p_multiple": general_metrics['ks_p_multiple'],
        "ks_multiple_d_mean": general_metrics['ks_multiple_d_mean'],
        "ks_multiple_p_mean": general_metrics['ks_multiple_p_mean']
    }


# In[15]:


def get_MGP_summary2(MGP):
    """
    Returns a string with general results of all averages stored in the MGP dictionary.
    - Formats implementation_num, TOP_N, num_files and MGP to show 3 digits on left (except MGP which has 3 before and 2 after decimal point).
    - Reorganizes columns so "Num Files" appears before "Mode".
    """
    result = "\n\nLegend:\n\n"
    # Header with reorganized columns
    result += "MGP = Weighted Global Presence Mean.\n"
    result += "Jaccard = Global Jaccard Coefficient Average.\n"
    result += "Fleiss = Fleiss Kappa Agreement Index.\n"
    result += "R-Kendall = Positive Kendall Tau Ratio.\n"
    result += "KS_p = Average p-value of Kolmogorov-Smirnov Test (KSp) for Multiple Samples\n"
    result += "Impl = Implementation Number.\n"
    result += "Num Files = Number of Files in Implementation.\n"
    result += "Met = Implementation Metric.\n          1. Degree\n          2. Betweenness Centrality\n          3. Bridging Centrality\n\n"
    result += "Mode: Features.\n\n"
    result += "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<10} {:<10} \n".format(
        "MGP", "Jaccard", "Fleiss", "R-Kendall", "KS_p", "TOP_N", "Num Files", "Met", "Impl"
    )
    result += "-" * 130 + "\n"
    
    for implementation_num, data in MGP.items():
        # Format implementation_num, TOP_N and num_files with 3 left digits
        formatted_impl = f"{implementation_num:03}"
        formatted_top_n = f"{data['TOP_N']:03}"
        formatted_num_files = f"{data['num_files']:03}"
        
        # Format MGP with 3 digits before and 2 after decimal point
        formatted_mgp = f"{data['MGP']:06.2f}%"  # 6 digits total (3 before + 1 point + 2 after)
        
        # Add formatted data to result
        result += "{:<10} {:<10.4f} {:<12.4f} {:<10.4f} {:<12.4f} {:<15} {:<13} {:<12} {:<12} \n".format(
            formatted_mgp, data["mean_jaccard"], data["fleiss_kappa"], data["positive_kendall_ratio"],
            data["d_global_mean"], formatted_top_n, formatted_num_files, data["Metric"], formatted_impl
        )
    
    return result


# In[16]:


def get_MGP_summary(MGP):
    """
    Returns a string with general results of all averages stored in the MGP dictionary.
    - Formats implementation_num, TOP_N, num_files and MGP to show 3 digits on left (except MGP which has 3 before and 2 after decimal point).
    - Reorganizes columns so "Num Files" appears before "Mode".
    """
    result = "\n\nLegend:\n\n"
    # Header with reorganized columns
    result += "MGP = Weighted Global Presence Mean.\n"
    result += "Jaccard = Global Jaccard Coefficient Average.\n"
    result += "Fleiss = Fleiss Kappa Agreement Index.\n"
    result += "R-Kendall = Positive Kendall Tau Ratio.\n"
    result += "KS_Multi = Average p-value of Kolmogorov-Smirnov Test (KSp) for Multiple Samples\n"
    result += "KS_p = Kolmogorov-Smirnov Test (KSp) for Global Distribution\n"
    result += "Impl = Implementation Number.\n"
    result += "Num Files = Number of Files in Implementation.\n"
    result += "Met = Implementation Metric.\n          1. Degree\n          2. Betweenness Centrality\n          3. Bridging Centrality\n\n"
    result += "Mode: Features.\n\n"
    result += "{:<12} {:<12} {:<12} {:<12} {:<10} {:<10} {:<10} {:<8} {:<8} {:<8} \n".format(
        "MGP", "Jaccard", "Fleiss", "R-Kendall", "KS_p", "KS_Multi", "TOP_N", "Num Files", "Met", "Impl"
    )
    result += "-" * 130 + "\n"
    
    for implementation_num, data in MGP.items():
        # Format implementation_num, TOP_N and num_files with 3 left digits
        formatted_impl = f"{implementation_num:03}"
        formatted_top_n = f"{data['TOP_N']:03}"
        formatted_num_files = f"{data['num_files']:03}"
        
        # Format MGP with 3 digits before and 2 after decimal point
        formatted_mgp = f"{data['MGP']:06.2f}%"  # 6 digits total (3 before + 1 point + 2 after)
        
        # Add formatted data to result
        result += "{:<10} {:<10.4f} {:<12.4f} {:<10.4f} {:<12.4f} {:<12.4f} {:<10} {:<10} {:<10} {:<10} \n".format(
            formatted_mgp, data["mean_jaccard"], data["fleiss_kappa"], data["positive_kendall_ratio"],
            data["d_global_mean"], data["ks_multiple_p_mean"], formatted_top_n, formatted_num_files, 
            data["Metric"], formatted_impl
        )
    
    return result


# In[17]:


def print_title_page(pdf_pages, implementation_num):
    # Page 0: Title Page
    plt.figure(figsize=(8, 11))  # Page size (8x11 inches)
    
    # Add centered title
    plt.text(0.5, 0.7, f"Implementation {implementation_num}\n", 
             fontsize=48, 
             color='blue',              
             ha='center',  # Horizontal alignment
             va='center',  # Vertical alignment
             weight='bold')
             
    plt.text(0.5, 0.6, "Similarity, Correlation \n and Distribution Tests", 
             fontsize=37, 
             color='blue',              
             ha='center',
             va='center',
             weight='bold')
             
    # Add centered subtitle
    plt.text(0.5, 0.3, f"Mode: {mode_name} \nMetric: {metric_name}", 
             fontsize=32, 
             color='red',             
             ha='center',
             va='center',
             style='italic')
             
    plt.text(0.5, 0.0, f"Top Nodes: {TOP_N} \nNumber of Files: {num_files_to_print}", 
             fontsize=16, 
             color='black',             
             ha='center',
             va='center',
             weight='bold')
              
    # Remove axes and borders
    plt.axis('off')    
    
    # Save page to PDF
    pdf_pages.savefig()
    plt.close()


# In[18]:


from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

def create_pdf(text, title, output_filename="Global_Results.pdf", paginate=True):
    """
    Creates a PDF with the provided text, splitting it across multiple pages if needed.
    
    :param text: Text to include in the PDF
    :param title: PDF title
    :param paginate: If True, splits text across multiple pages
    """
    # Page settings
    lines_per_page = 50  # Number of lines that fit on one page
    fontsize = 10  # Text font size
    top_margin = 0.9  # Top margin for title
    bottom_margin = 0.1  # Bottom margin

    # Split text into lines
    lines = text.split('\n')

    # Initialize PDF
    pdf_pages = PdfPages(output_filename)

    # Page creation loop
    current_page = 1
    while lines:
        # Create new page
        plt.figure(figsize=(8, 11))
        plt.figtext(0.5, top_margin, title, fontsize=12, wrap=True, 
                   weight='bold', va='top', ha='center', color='green')

        # Select lines for current page
        page_lines = lines[:lines_per_page]
        lines = lines[lines_per_page:]

        # Add text to page
        page_text = "\n".join(page_lines)
        plt.figtext(0.1, top_margin - 0.05, page_text, 
                   fontsize=fontsize, wrap=True, va='top')

        # Turn off axes
        plt.axis('off')

        # Save page to PDF
        pdf_pages.savefig()
        plt.close()

        # Update page number
        current_page += 1

    # Close PDF
    pdf_pages.close()


# In[19]:


def generate_commands(mode, output_filename, num_files_list, mode_list, metric_list, top_n_list, start_implementation_num, combinations):
    """
    Generates and executes commands for running calculations with different parameter combinations.
    
    Args:
        mode: Execution mode
        output_filename: Output PDF filename
        num_files_list: List of number of files to process
        mode_list: List of processing modes
        metric_list: List of metrics to calculate
        top_n_list: List of top N values
        start_implementation_num: Starting implementation number
        combinations: Total number of combinations
    """
    implementation_num = start_implementation_num
    execution_count = 0
    for processing_mode in mode_list:
        for num_files in num_files_list:
            for metric in metric_list:
                for top_n in top_n_list:
                    run_calculations(mode, top_n, processing_mode, metric, num_files, implementation_num, output_filename)
                    implementation_num += 1
                    execution_count += 1
                    print(f"{((execution_count) / combinations * 100):.2f} % completed." + " " * 5)                    
    print("Execution finished. \n Final results files completed.")
    print("Procedure completed. File saved or Error occurred.")


# In[20]:


import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from sklearn.preprocessing import MinMaxScaler

def generate_and_save_plot(MGP, output_pdf="final_plot.pdf"):
    """
    Generates a scatter plot with multiple series from MGP dictionary data.
    Normalizes metrics to the same scale (0 to 1).
    Saves the plot to a PDF file.

    Parameters:
    - MGP: Dictionary containing metric data
    - output_pdf: Output PDF filename (default: "final_plot.pdf")
    """
    # Create DataFrame from MGP dictionary
    data = []
    for impl_num, values in MGP.items():
        data.append({
            "Top_N": values["TOP_N"],
            "MGP": values["MGP"],
            "Jaccard": values["mean_jaccard"],
            "Fleiss": values["fleiss_kappa"],
            "R_Kendall": values["positive_kendall_ratio"],
            "KS_Multi": values["ks_multiple_d_mean"],
            "Metric": values["Metric"]
        })
    df = pd.DataFrame(data)

    # Normalize metrics to [0, 1] scale
    scaler = MinMaxScaler()
    df[["MGP", "Jaccard", "Fleiss", "R_Kendall", "KS_Multi"]] = scaler.fit_transform(
        df[["MGP", "Jaccard", "Fleiss", "R_Kendall", "KS_Multi"]]
    )

    # Plot configuration
    plt.figure(figsize=(12, 8))
    sns.set_style("whitegrid")
    sns.set_palette("tab10")

    # Plot series
    sns.scatterplot(data=df, x="Top_N", y="MGP", label="MGP", s=100, marker="o")
    sns.scatterplot(data=df, x="Top_N", y="Jaccard", label="Jaccard", s=100, marker="s")
    sns.scatterplot(data=df, x="Top_N", y="Fleiss", label="Fleiss", s=100, marker="D")
    sns.scatterplot(data=df, x="Top_N", y="R_Kendall", label="R-Kendall", s=100, marker="X")
    sns.scatterplot(data=df, x="Top_N", y="KS_Multi", label="KS_Multi", s=100, marker="*")

    # Trend lines (no confidence interval)
    sns.lineplot(data=df, x="Top_N", y="MGP", errorbar=None, color="blue", alpha=0.5)
    sns.lineplot(data=df, x="Top_N", y="Jaccard", errorbar=None, color="orange", alpha=0.5)
    sns.lineplot(data=df, x="Top_N", y="Fleiss", errorbar=None, color="green", alpha=0.5)
    sns.lineplot(data=df, x="Top_N", y="R_Kendall", errorbar=None, color="red", alpha=0.5)
    sns.lineplot(data=df, x="Top_N", y="KS_Multi", errorbar=None, color="purple", alpha=0.5)

    # Aesthetic adjustments
    plt.title("Trend of Normalized Metrics vs Top_N: \nMultimetric Analysis with Topological Implementations in Machine Learning Mode", fontsize=14)
    plt.xlabel("Top_N (Number of Top Nodes in File Rankings)", fontsize=12)
    plt.ylabel("Normalized Metric Value", fontsize=12)

    # Legend configuration
    legend = plt.legend(
        title="Metrics", 
        loc='upper left', 
        bbox_to_anchor=(-0.005, 1.007),
        frameon=True, 
        framealpha=0.9,
        facecolor="#F5F5F5",
        borderpad=1,
        labelspacing=0.6,
        handletextpad=0.5,
        handlelength=2
    )
    
    legend.get_frame().set_width(2.2)
    legend.get_frame().set_height(4.8)
    
    # Explanatory notes
    notes = [
    "                              - Reference Values -\n",
    "• MGP (Weighted Global Presence Mean):",
    "  - Desirable: Close to 1.0",
    "  - Maximum: 1.0 (Total ranking consistency)",
    "  - Minimum: 0.0 (Total inconsistency)\n",
    "• Jaccard (Similarity Index):",
    "  - Desirable: Close to 1.0",
    "  - Maximum: 1.0 (Identical sets)",
    "  - Minimum: 0.0 (No intersection between sets)\n",
    "• R-Kendall (Positive Correlation Ratio):",
    "  - Desirable: > 0.8 (High proportion of concordant pairs)",
    "  - Maximum: 1.0 (All pairs have positive correlation)",
    "  - Minimum: 0.0 (No pairs with positive correlation)\n",
    "• Fleiss (Fleiss' Kappa):",
    "  - Desirable > 0.6",
    "  - Maximum: 1.0 (Perfect agreement between raters)",
    "  - Minimum: 0.0 (Agreement equivalent to chance)\n",
    "• KS_Multi (KS Distance (Kolmogorov-Smirnov)",
    "                    for Multiple Samples):",    
    "  - Desirable: Close to 0.0",
    "  - Maximum: 1.0 (Completely different distributions)",
    "  - Minimum: 0.0 (Identical distributions)\n"
    ]
    
    plt.figtext(
        0.615, 0.35, '\n'.join(notes),
        ha='left', va='center', fontsize=8, wrap=True,
        bbox=dict(boxstyle="round,pad=0.3", edgecolor="gray", facecolor="white"),
    )

    plt.tight_layout()

    # Save plot to PDF
    plt.savefig(output_pdf, format="pdf", bbox_inches="tight")
    plt.close()

    print(f"Plot successfully saved to {output_pdf}!")


# In[21]:


def main():
    global title, top7, presence_percentage, combined_metrics, combined_heatmaps, final_filename, MGP
    
    # Configuration of pages to include
    title = True
    top7 = True
    presence_percentage = True
    combined_metrics = True
    combined_heatmaps = True
    final_filename = "Combined_Results"
    text = ""
    
    # Predefined parameters for batch execution
    # num_files_list = [50]
    # mode_list = [2]
    # metric_list = [1]
    # top_n_list = [100, 200]
    num_files_list = [10, 20, 30, 40, 50]
    mode_list = [1]  # Machine Learning mode
    metric_list = [1, 2, 3]  # Degree, Betweenness, Bridging
    top_n_list = [10, 20, 30, 50, 100, 200]

    # Calculate total number of combinations
    total_combinations = len(num_files_list) * len(mode_list) * len(metric_list) * len(top_n_list)

    choice = input("Select execution mode:\n1 - Single\n2 - Batch\n3 - Global Results\nOption: ").strip()

    if choice == '1':
        filename = input("Output filename (e.g.: Result.pdf): ")
        # Predefined parameters for single execution
        run_calculations(mode=0, top_nodes=20, calculation_mode=1, metric_choice=3, 
                        num_files=40, implementation_num=68, output_filename=filename)        
    elif choice == '2':
        final_filename = input("Final output filename (e.g.: Combined_Results.pdf): ")        
        generate_commands(mode=2, output_filename="Temporary.pdf", 
                         num_files_list=num_files_list, 
                         mode_list=mode_list, 
                         metric_list=metric_list, 
                         top_n_list=top_n_list, 
                         start_implementation_num=1,
                         combinations=total_combinations)
    elif choice == '3':
        final_filename = input("Final output filename (e.g.: Global_Results.pdf): ")
        plot_filename = input("Final plot filename (e.g.: Plot.pdf): ")
        generate_commands(mode=3, output_filename="Temporary", 
                         num_files_list=num_files_list, 
                         mode_list=mode_list, 
                         metric_list=metric_list, 
                         top_n_list=top_n_list, 
                         start_implementation_num=1,
                         combinations=total_combinations)                                    
        text = get_MGP_summary(MGP)        
        create_pdf(text, "Global Results of Similarity, Correlation\nand Distribution Tests for Machine Learning Mode", final_filename)
        generate_and_save_plot(MGP, output_pdf=plot_filename)
    else:
        print("Invalid option.")
    print("Procedure completed successfully.\n")

if __name__ == "__main__":
    main()    


# In[ ]:




