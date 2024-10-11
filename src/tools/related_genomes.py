#!/usr/bin/python3
import sys
from collections import Counter

def count_similar_genomes(input_file, target_genome):
    """
    Count the occurrences of similar genomes to the target genome.
    
    Args:
    input_file (str): Path to the input file (pan-genome/pan-genome-format-1.txt)
    target_genome (str): The reference genome to compare against
    
    Returns:
    Counter: A Counter object with genome names as keys and their occurrence counts as values
    """
    genome_counts = Counter()
    target_scope = False

    with open(input_file, 'r') as file:
        for line in file:
            if line.startswith("Protein:") and line.endswith(target_genome + "\n"):
                target_scope = True
            elif target_scope and line.startswith("Protein:") and "Genome:" in line:
                genome_name = line.split("Genome:")[1].strip()
                if genome_name != target_genome + ";":
                    genome_counts[genome_name] += 1
            elif line.startswith("----"):
                target_scope = False

    return genome_counts

def print_help():
    """Print the help message explaining the program's functionality and usage."""
    help_message = """
    Related Genomes Counter

    Description:
    This program analyzes a pan-genome file and generates a sorted list of genomes that are 
    most similar to a specified reference genome. This information can help decide which 
    genomes should be used to generate an interaction network. Highly similar genomes should 
    be removed from the genppi execution directory to avoid excessively high degree nodes 
    or proteins.

    Usage:
    python related_genomes.py <input_file> <target_genome>

    Arguments:
    <input_file>    : Path to the input file (expected to be pan-genome/pan-genome-format-1.txt)
    <target_genome> : The reference genome to compare against

    Output:
    A sorted list of genomes most similar to the target genome, along with their occurrence counts.

    Example:
    python related_genomes.py pan-genome/pan-genome-format-1.txt Escherichia_coli_K12
    """
    print(help_message)

def main():
    if len(sys.argv) != 3:
        print_help()
        sys.exit(1)

    input_file = sys.argv[1]
    target_genome = sys.argv[2]

    genome_counts = count_similar_genomes(input_file, target_genome)

    print(f"Genomes most similar to {target_genome}:")
    for genome, count in genome_counts.most_common():
        if genome != target_genome:
            print(f"{genome[:-1]}\t{count}")

if __name__ == "__main__":
    main()