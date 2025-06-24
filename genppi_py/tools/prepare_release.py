#!/usr/bin/env python3
"""
Script to prepare release files for GenPPI Python interface
"""

import os
import zipfile
import shutil
from pathlib import Path

# Path to the root of the GenPPI repository
GENPPI_ROOT = Path(__file__).parent.parent.absolute()

# Path to the binaries directory
BINARIES_DIR = GENPPI_ROOT / 'binaries'

# Path to the src directory (for model.dat)
SRC_DIR = GENPPI_ROOT / 'src'

# Path to the output directory for release files
OUTPUT_DIR = Path(__file__).parent / 'release'

# Create the output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Function to create a zip file for a specific platform
def create_platform_zip(platform, executables):
    """Create a zip file containing executables for a specific platform"""
    zip_filename = f"genppi_{platform.lower()}_executables.zip"
    zip_path = OUTPUT_DIR / zip_filename
    
    print(f"Creating {zip_filename}...")
    
    with zipfile.ZipFile(zip_path, 'w') as zipf:
        for executable in executables:
            executable_path = BINARIES_DIR / executable
            if executable_path.exists():
                print(f"  Adding {executable}")
                zipf.write(executable_path, executable)
            else:
                print(f"  Warning: {executable} not found")
    
    print(f"Created {zip_filename}")

# Create zip files for each platform
create_platform_zip('linux', [
    'genppi32g-Linux',
    'genppidb32g-Linux'
])

create_platform_zip('mac', [
    'genppi32g-Mac.x',
    'genppidb32g-Mac.x'
])

create_platform_zip('windows', [
    'genppi32g-Win.exe',
    'genppidb32g-Win.exe'
])

# Copy model.dat file
model_src = SRC_DIR / 'model.dat'
model_dst = OUTPUT_DIR / 'model.dat'
if model_src.exists():
    print(f"Copying model.dat...")
    shutil.copy2(model_src, model_dst)
    print(f"Copied model.dat to {model_dst}")
else:
    print(f"Warning: model.dat not found at {model_src}")

print("\nRelease files prepared successfully!")
print(f"Files are available in: {OUTPUT_DIR}")
print("\nPlease upload these files to the GitHub releases page:")
for file in OUTPUT_DIR.glob('*'):
    print(f"- {file.name}")