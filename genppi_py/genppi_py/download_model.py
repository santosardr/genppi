#!/usr/bin/env python3
"""
Download and extract model.dat file for GenPPI
"""

import os
import sys
import tempfile
import urllib.request
import shutil
from pathlib import Path
# NOVO: Import para localizar arquivos de forma robusta
from importlib import resources

# Função auxiliar para encontrar o diretório 'bin' de forma robusta
def _get_package_bin_dir():
    """Retorna o Path para o diretório 'bin' dentro do pacote."""
    try:
        # Método moderno e preferencial para pacotes instalados
        return resources.files('genppi_py') / 'bin'
    except (AttributeError, ModuleNotFoundError):
        # Fallback para ambientes de desenvolvimento
        return Path(__file__).parent.absolute() / 'bin'

def download_model_file(bin_dir=None):
    """
    Download and extract the model.dat file from 7zip parts
    
    Args:
        bin_dir (str or Path, optional): Directory to save the model.dat file.
            If None, use the bin directory in the package.
    
    Returns:
        bool: True if successful, False otherwise
    """
    if bin_dir is None:
        bin_dir = _get_package_bin_dir()
    else:
        bin_dir = Path(bin_dir)
    
    model_path = bin_dir / "model.dat"
    
    if model_path.exists():
        print("model.dat file already exists.")
        return True
    
    print("Model file not found. Downloading and extracting model.dat...")
    
    github_raw_url = "https://github.com/santosardr/genppi/raw/master/src/"
    model_7z_parts = ["model.7z.001", "model.7z.002", "model.7z.003"]
    
    os.makedirs(bin_dir, exist_ok=True)
    
    temp_dir_path = bin_dir / "model_parts"
    os.makedirs(temp_dir_path, exist_ok=True)
    print(f"Saving 7z parts to {temp_dir_path}")
    
    try:
        for part in model_7z_parts:
            part_url = f"{github_raw_url}{part}"
            part_path = temp_dir_path / part
            
            try:
                print(f"Downloading {part}...")
                urllib.request.urlretrieve(part_url, part_path)
                print(f"Successfully downloaded {part}")
            except Exception as e:
                print(f"Error downloading {part}: {e}", file=sys.stderr)
                shutil.rmtree(temp_dir_path, ignore_errors=True)
                return False
        
        print("All parts downloaded. Extracting model.dat...")
        
        # Extract using py7zr with multivolume support
        extracted = False
        
        try:
            print("Extracting model.dat using py7zr...")
            import py7zr
            from multivolumefile import MultiVolume
            
            # Get the base name for multivolume (remove .001 extension)
            first_part_path = temp_dir_path / model_7z_parts[0]
            archive_base_name = str(first_part_path).rsplit('.', 1)[0]
            
            # Extract using multivolume support
            with MultiVolume(archive_base_name, mode='rb') as vol:
                with py7zr.SevenZipFile(vol, mode='r') as archive:
                    # Get list of files in archive
                    filenames_in_archive = archive.getnames()
                    
                    if len(filenames_in_archive) != 1:
                        raise Exception(f"Expected exactly 1 file in archive, found {len(filenames_in_archive)}: {filenames_in_archive}")
                    
                    original_filename = filenames_in_archive[0]
                    print(f"Found file in archive: {original_filename}")
                    
                    # Extract to bin directory
                    archive.extractall(path=str(bin_dir))
                    
                    # Check if we need to rename the extracted file
                    extracted_file = bin_dir / original_filename
                    if extracted_file.exists() and original_filename != "model.dat":
                        print(f"Renaming {original_filename} to model.dat")
                        extracted_file.rename(model_path)
            
            if model_path.exists():
                print(f"model.dat successfully extracted to {model_path}")
                extracted = True
            else:
                raise Exception("model.dat was not found after extraction")
                
        except Exception as e:
            print(f"Error with py7zr extraction: {e}", file=sys.stderr)

        # Limpa os arquivos .7z baixados após a extração
        shutil.rmtree(temp_dir_path, ignore_errors=True)
        
        if extracted:
            print("model.dat extraction successful.")
            return True
        else:
            print("\nExtraction failed.", file=sys.stderr)
            print("Please extract manually using the following steps:", file=sys.stderr)
            print("1. Ensure py7zr is installed: pip install --upgrade --force-reinstall \"py7zr>=0.20.0,<1.0.0\"", file=sys.stderr)
            print("2. Download the 7z parts to the same directory:", file=sys.stderr)
            for part in model_7z_parts:
                print(f"   - {github_raw_url}{part}", file=sys.stderr)
            print("3. Extract using py7zr with the following Python code:", file=sys.stderr)
            print("   import py7zr", file=sys.stderr)
            print("   from multivolumefile import MultiVolume", file=sys.stderr)
            print("   archive_base = 'model.7z'  # without .001 extension", file=sys.stderr)
            print("   with MultiVolume(archive_base, mode='rb') as vol:", file=sys.stderr)
            print("       with py7zr.SevenZipFile(vol, mode='r') as archive:", file=sys.stderr)
            print("           archive.extractall(path='.')", file=sys.stderr)
            print(f"4. Place the extracted model.dat file in: {bin_dir}", file=sys.stderr)
            return False
            
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        shutil.rmtree(temp_dir_path, ignore_errors=True)
        return False

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Download and extract model.dat file for GenPPI")
    parser.add_argument("--output-dir", "-o", help="Directory to save the model.dat file")
    
    args = parser.parse_args()
    
    success = download_model_file(args.output_dir)
    
    if success:
        print("\nModel file download and extraction completed successfully.")
        return 0
    else:
        print("\nModel file download or extraction failed.", file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())