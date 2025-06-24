#!/usr/bin/env python3
"""
Setup script for GenPPI Python interface
"""

import os
import sys
import platform
import urllib.request
import tempfile
import zipfile
from setuptools import setup, find_packages
from setuptools.command.install import install
import shutil
from pathlib import Path


class CustomInstall(install):
    """
    Custom installation class to download and install the appropriate executables
    """
    # Base URL for downloading executables directly from GitHub repository
    GITHUB_RAW_URL = "https://raw.githubusercontent.com/santosardr/genppi/master/"
    GITHUB_API_URL = "https://api.github.com/repos/santosardr/genppi/contents/"
    
    # Mapping of platform to executable names
    PLATFORM_EXECUTABLES = {
        "Linux": ["genppi32g-Linux", "genppidb32g-Linux"],
        "Darwin": ["genppi32g-Mac.x", "genppidb32g-Mac.x"],  # macOS
        "Windows": ["genppi32g-Win.exe", "genppidb32g-Win.exe"]
    }
    
    # Model.dat 7zip parts
    MODEL_7Z_PARTS = ["model.7z.001", "model.7z.002", "model.7z.003"]
    
    def run(self):
        # Run the standard installation
        install.run(self)
        
        # Get the package installation directory
        install_dir = Path(self.install_lib) / 'genppi_py'
        bin_dir = install_dir / 'bin'
        
        # Create the bin directory if it doesn't exist
        os.makedirs(bin_dir, exist_ok=True)
        
        # Determine the current platform
        current_platform = platform.system()
        if current_platform not in self.PLATFORM_EXECUTABLES:
            print(f"Warning: Unsupported platform {current_platform}. GenPPI executables may not work.")
            return
        
        # Try to download and install executables
        self.download_and_install_executables(bin_dir, current_platform)
        
        # Try to download model.dat file
        self.download_model_file(bin_dir)
    
    def download_and_install_executables(self, bin_dir, platform):
        """Download and install executables for the current platform"""
        executables = self.PLATFORM_EXECUTABLES.get(platform, [])
        if not executables:
            print(f"Warning: No executables defined for platform {platform}")
            return
        
        print(f"Downloading GenPPI executables for {platform}...")
        
        # Create the bin directory if it doesn't exist
        os.makedirs(bin_dir, exist_ok=True)
        
        # Download each executable
        for executable in executables:
            download_url = f"{self.GITHUB_RAW_URL}binaries/{executable}"
            output_path = bin_dir / executable
            
            try:
                print(f"  Downloading {executable}...")
                urllib.request.urlretrieve(download_url, output_path)
                
                # Make executable on Unix-like systems
                if platform != "Windows" and not executable.endswith('.exe'):
                    output_path.chmod(0o755)
                
                print(f"  Successfully downloaded {executable}")
            except Exception as e:
                print(f"  Error downloading {executable}: {e}")
                print(f"  You may need to manually download {executable} from:")
                print(f"  {download_url}")
                print(f"  and place it in: {bin_dir}")
        
        print(f"Finished downloading executables for {platform}")
    
    def download_model_file(self, bin_dir):
        """Download and extract the model.dat file from 7zip parts"""
        model_path = bin_dir / "model.dat"
        
        # Skip if model.dat already exists
        if model_path.exists():
            print("model.dat file already exists, skipping download")
            return
        
        print("Downloading and extracting model.dat file...")
        
        # Create a temporary directory for the 7z parts
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir_path = Path(temp_dir)
            
            # Download all 7z parts
            for part in self.MODEL_7Z_PARTS:
                part_url = f"{self.GITHUB_RAW_URL}src/{part}"
                part_path = temp_dir_path / part
                
                try:
                    print(f"  Downloading {part}...")
                    urllib.request.urlretrieve(part_url, part_path)
                    print(f"  Successfully downloaded {part}")
                except Exception as e:
                    print(f"  Error downloading {part}: {e}")
                    print(f"  You may need to manually download the model.dat file.")
                    return
            
            # Extract the model.dat file using py7zr
            try:
                print("  Extracting model.dat using py7zr...")
                import py7zr
                from multivolumefile import MultiVolume
                
                # Get the base name for multivolume (remove .001 extension)
                first_part_path = temp_dir_path / self.MODEL_7Z_PARTS[0]
                archive_base_name = str(first_part_path).rsplit('.', 1)[0]
                
                # Extract using multivolume support
                with MultiVolume(archive_base_name, mode='rb') as vol:
                    with py7zr.SevenZipFile(vol, mode='r') as archive:
                        # Get list of files in archive
                        filenames_in_archive = archive.getnames()
                        
                        if len(filenames_in_archive) != 1:
                            raise Exception(f"Expected exactly 1 file in archive, found {len(filenames_in_archive)}: {filenames_in_archive}")
                        
                        original_filename = filenames_in_archive[0]
                        print(f"  Found file in archive: {original_filename}")
                        
                        # Extract to bin directory
                        archive.extractall(path=str(bin_dir))
                        
                        # Check if we need to rename the extracted file
                        extracted_file = bin_dir / original_filename
                        if extracted_file.exists() and original_filename != "model.dat":
                            print(f"  Renaming {original_filename} to model.dat")
                            extracted_file.rename(model_path)
                
                print("  Successfully extracted model.dat file")
                
                # Verify that model.dat was extracted
                if not model_path.exists():
                    raise Exception("model.dat was not extracted properly")
                
            except Exception as e:
                print(f"  Error extracting model.dat: {e}")
                print("  Automatic extraction failed. Please extract manually:")
                print("  1. Download the following files to the same directory:")
                for part in self.MODEL_7Z_PARTS:
                    print(f"     - {self.GITHUB_RAW_URL}src/{part}")
                print("  2. Extract using py7zr with the following Python code:")
                print("     import py7zr")
                print("     from multivolumefile import MultiVolume")
                print("     archive_base = 'model.7z'  # without .001 extension")
                print("     with MultiVolume(archive_base, mode='rb') as vol:")
                print("         with py7zr.SevenZipFile(vol, mode='r') as archive:")
                print("             archive.extractall(path='.')")
                print(f"  3. Place the extracted model.dat file in: {bin_dir}")
                print("  4. Ensure py7zr is installed: pip install --upgrade --force-reinstall \"py7zr>=0.20.0,<1.0.0\"")


setup(
    cmdclass={
        'install': CustomInstall,
    },
)