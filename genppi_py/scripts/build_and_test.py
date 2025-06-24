#!/usr/bin/env python3
"""
Build and test script for GenPPI Python interface
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path

def run_command(cmd, description):
    """Run a command and handle errors"""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Command: {' '.join(cmd)}")
    print(f"{'='*60}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("SUCCESS")
        if result.stdout:
            print("STDOUT:", result.stdout)
        return True
    except subprocess.CalledProcessError as e:
        print("FAILED")
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
        return False
    except FileNotFoundError:
        print(f"FAILED: Command not found: {cmd[0]}")
        return False

def clean_build_artifacts():
    """Clean build artifacts"""
    print("\nCleaning build artifacts...")
    artifacts = [
        "build",
        "dist", 
        "*.egg-info",
        "__pycache__",
        "genppi_py/__pycache__",
        "genppi_py/bin/__pycache__"
    ]
    
    for artifact in artifacts:
        if "*" in artifact:
            # Handle glob patterns
            import glob
            for path in glob.glob(artifact):
                if os.path.exists(path):
                    if os.path.isdir(path):
                        shutil.rmtree(path)
                    else:
                        os.remove(path)
                    print(f"  Removed: {path}")
        else:
            if os.path.exists(artifact):
                if os.path.isdir(artifact):
                    shutil.rmtree(artifact)
                else:
                    os.remove(artifact)
                print(f"  Removed: {artifact}")

def main():
    """Main build and test function"""
    print("GenPPI Python Interface - Build and Test Script")
    print("=" * 60)
    
    # Change to project root (parent of scripts directory)
    project_root = Path(__file__).parent.parent
    os.chdir(project_root)
    
    success = True
    
    # Clean previous builds
    clean_build_artifacts()
    
    # 1. Check Python syntax
    success &= run_command([
        sys.executable, "-m", "py_compile", "setup.py"
    ], "Check setup.py syntax")
    
    success &= run_command([
        sys.executable, "-m", "py_compile", "genppi_py/genppi.py"
    ], "Check main module syntax")
    
    success &= run_command([
        sys.executable, "-m", "py_compile", "genppi_py/download_model.py"
    ], "Check download_model syntax")
    
    # 2. Check dependencies
    success &= run_command([
        sys.executable, "-c", "import py7zr; from multivolumefile import MultiVolume; print('Dependencies OK')"
    ], "Check required dependencies")
    
    # 3. Build source distribution
    success &= run_command([
        sys.executable, "-m", "build", "--sdist"
    ], "Build source distribution")
    
    # 4. Build wheel
    success &= run_command([
        sys.executable, "-m", "build", "--wheel"
    ], "Build wheel distribution")
    
    # 5. Check distributions
    success &= run_command([
        sys.executable, "-m", "twine", "check", "dist/*"
    ], "Check distribution files")
    
    # 6. Test installation in virtual environment (if available)
    if shutil.which("python3") and shutil.which("pip"):
        print("\n" + "="*60)
        print("Testing installation in clean environment...")
        print("="*60)
        
        # Create temporary virtual environment
        venv_path = project_root / "test_venv"
        if venv_path.exists():
            shutil.rmtree(venv_path)
        
        success &= run_command([
            sys.executable, "-m", "venv", str(venv_path)
        ], "Create test virtual environment")
        
        if success:
            # Determine venv python path
            if os.name == 'nt':  # Windows
                venv_python = venv_path / "Scripts" / "python.exe"
                venv_pip = venv_path / "Scripts" / "pip.exe"
            else:  # Unix-like
                venv_python = venv_path / "bin" / "python"
                venv_pip = venv_path / "bin" / "pip"
            
            # Install package in venv
            success &= run_command([
                str(venv_pip), "install", "--upgrade", "pip"
            ], "Upgrade pip in test environment")
            
            success &= run_command([
                str(venv_pip), "install", "."
            ], "Install package in test environment")
            
            # Test basic functionality
            if success:
                success &= run_command([
                    str(venv_python), "-c", "import genppi_py; print('Package import OK')"
                ], "Test package import")
                
                success &= run_command([
                    str(venv_python), "-c", "from genppi_py.genppi import main; print('Main function import OK')"
                ], "Test main function import")
        
        # Clean up test environment
        if venv_path.exists():
            shutil.rmtree(venv_path)
            print("Cleaned up test virtual environment")
    
    # 7. Run unit tests
    if os.path.exists("tests/test_genppi.py"):
        success &= run_command([
            sys.executable, "tests/test_genppi.py"
        ], "Run unit tests")
    
    # Summary
    print("\n" + "="*60)
    if success:
        print("✅ BUILD AND TEST SUCCESSFUL!")
        print("\nGenerated files:")
        if os.path.exists("dist"):
            for file in os.listdir("dist"):
                print(f"  - dist/{file}")
        print("\nReady for deployment!")
    else:
        print("❌ BUILD AND TEST FAILED!")
        print("Please fix the issues above before deployment.")
    print("="*60)
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())