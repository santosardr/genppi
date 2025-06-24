#!/usr/bin/env python3
"""
Download sample files for GenPPI testing
"""

import os
import urllib.request
from pathlib import Path

def download_samples(output_dir='samples'):
    """
    Download sample files from GitHub and save them to the specified directory
    
    Args:
        output_dir (str): Directory to save the sample files
    """
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Base URL for the sample files
    base_url = "https://raw.githubusercontent.com/santosardr/genppi/master/test/Buchnera_aphidicola/assemblies/"
    
    # List of sample files to download
    sample_files = [
        "Ba_Ak.faa",
        "Ba_Bp.faa",
        "Ba_G002.faa",
        "Ba_Sg.faa",
        "Ba_Ua.faa"
    ]
    
    print(f"Downloading sample files to {output_dir}/...")
    
    # Download each sample file
    for file in sample_files:
        url = base_url + file
        output_path = os.path.join(output_dir, file)
        
        print(f"  Downloading {file}...")
        try:
            urllib.request.urlretrieve(url, output_path)
            print(f"  Successfully downloaded {file}")
        except Exception as e:
            print(f"  Error downloading {file}: {e}")
    
    print("Sample files download complete.")

def main():
    """
    Command-line entry point
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="Download sample files for GenPPI testing")
    parser.add_argument("--output-dir", "-o", default="samples", help="Directory to save the sample files")
    
    args = parser.parse_args()
    
    download_samples(args.output_dir)
    
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())