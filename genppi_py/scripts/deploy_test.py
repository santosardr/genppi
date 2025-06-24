#!/usr/bin/env python3
"""
Deployment script for GenPPI Python interface - Test Environment
"""

import os
import sys
import subprocess
import shutil
from pathlib import Path
import tempfile

def run_command(cmd, description, check=True):
    """Run a command and handle errors"""
    print(f"\n{'='*50}")
    print(f"Running: {description}")
    print(f"Command: {' '.join(cmd)}")
    print(f"{'='*50}")
    
    try:
        result = subprocess.run(cmd, check=check, capture_output=True, text=True)
        if result.returncode == 0:
            print("✅ SUCCESS")
        else:
            print("⚠️  WARNING: Non-zero exit code")
        
        if result.stdout:
            print("STDOUT:", result.stdout)
        if result.stderr and result.returncode != 0:
            print("STDERR:", result.stderr)
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print("❌ FAILED")
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
        return False
    except FileNotFoundError:
        print(f"❌ FAILED: Command not found: {cmd[0]}")
        return False

def check_prerequisites():
    """Check if all prerequisites are installed"""
    print("Checking prerequisites...")
    
    required_tools = [
        ("python", "Python interpreter"),
        ("pip", "Python package installer"),
        ("build", "Python build tool (pip install build)"),
        ("twine", "Python package uploader (pip install twine)")
    ]
    
    missing = []
    for tool, description in required_tools:
        if tool in ["build", "twine"]:
            # Check if these are available as Python modules
            try:
                subprocess.run([sys.executable, "-m", tool, "--help"], 
                             capture_output=True, check=True)
                print(f"✅ {description}")
            except (subprocess.CalledProcessError, FileNotFoundError):
                missing.append((tool, description))
                print(f"❌ {description}")
        else:
            if shutil.which(tool):
                print(f"✅ {description}")
            else:
                missing.append((tool, description))
                print(f"❌ {description}")
    
    if missing:
        print("\n❌ Missing prerequisites:")
        for tool, description in missing:
            if tool in ["build", "twine"]:
                print(f"  Install with: pip install {tool}")
            else:
                print(f"  {description}")
        return False
    
    print("✅ All prerequisites available")
    return True

def clean_environment():
    """Clean the build environment"""
    print("\nCleaning build environment...")
    
    artifacts = ["build", "dist", "*.egg-info"]
    for artifact in artifacts:
        if "*" in artifact:
            import glob
            for path in glob.glob(artifact):
                if os.path.exists(path):
                    shutil.rmtree(path) if os.path.isdir(path) else os.remove(path)
                    print(f"  Removed: {path}")
        else:
            if os.path.exists(artifact):
                shutil.rmtree(artifact) if os.path.isdir(artifact) else os.remove(artifact)
                print(f"  Removed: {artifact}")

def build_package():
    """Build the package"""
    print("\nBuilding package...")
    
    success = True
    
    # Build source distribution
    success &= run_command([
        sys.executable, "-m", "build", "--sdist"
    ], "Build source distribution")
    
    # Build wheel
    success &= run_command([
        sys.executable, "-m", "build", "--wheel"
    ], "Build wheel distribution")
    
    return success

def validate_package():
    """Validate the built package"""
    print("\nValidating package...")
    
    success = True
    
    # Check with twine
    success &= run_command([
        sys.executable, "-m", "twine", "check", "dist/*"
    ], "Validate package with twine")
    
    return success

def test_installation():
    """Test installation in a clean environment"""
    print("\nTesting installation...")
    
    with tempfile.TemporaryDirectory() as temp_dir:
        venv_path = Path(temp_dir) / "test_venv"
        
        # Create virtual environment
        success = run_command([
            sys.executable, "-m", "venv", str(venv_path)
        ], "Create test virtual environment")
        
        if not success:
            return False
        
        # Determine venv paths
        if os.name == 'nt':  # Windows
            venv_python = venv_path / "Scripts" / "python.exe"
            venv_pip = venv_path / "Scripts" / "pip.exe"
        else:  # Unix-like
            venv_python = venv_path / "bin" / "python"
            venv_pip = venv_path / "bin" / "pip"
        
        # Upgrade pip
        success &= run_command([
            str(venv_pip), "install", "--upgrade", "pip"
        ], "Upgrade pip in test environment")
        
        # Install our package
        wheel_files = list(Path("dist").glob("*.whl"))
        if wheel_files:
            wheel_file = wheel_files[0]
            success &= run_command([
                str(venv_pip), "install", str(wheel_file)
            ], f"Install package from wheel: {wheel_file.name}")
        else:
            print("❌ No wheel file found")
            return False
        
        # Test basic functionality
        if success:
            success &= run_command([
                str(venv_python), "-c", 
                "import genppi_py; print(f'GenPPI version: {genppi_py.__version__}')"
            ], "Test package import")
            
            success &= run_command([
                str(venv_python), "-c", 
                "from genppi_py.genppi import main; print('Main function available')"
            ], "Test main function import")
            
            # Test console scripts
            success &= run_command([
                str(venv_python), "-m", "genppi_py.genppi", "--help"
            ], "Test genppi module execution", check=False)
            
            # Test dependencies
            success &= run_command([
                str(venv_python), "-c", 
                "import py7zr; from multivolumefile import MultiVolume; print('Dependencies OK')"
            ], "Test required dependencies")
    
    return success

def deploy_to_test_pypi():
    """Deploy to Test PyPI (optional)"""
    print("\nDeployment to Test PyPI...")
    print("Note: This requires Test PyPI credentials")
    
    response = input("Do you want to upload to Test PyPI? (y/N): ").strip().lower()
    if response != 'y':
        print("Skipping Test PyPI upload")
        return True
    
    # Upload to Test PyPI
    success = run_command([
        sys.executable, "-m", "twine", "upload", 
        "--repository", "testpypi", 
        "dist/*"
    ], "Upload to Test PyPI", check=False)
    
    if success:
        print("\n✅ Package uploaded to Test PyPI!")
        print("Install with: pip install --index-url https://test.pypi.org/simple/ genppi-py")
    else:
        print("\n⚠️  Upload failed. Check credentials and try again.")
    
    return success

def main():
    """Main deployment function"""
    print("GenPPI Python Interface - Test Deployment Script")
    print("=" * 60)
    
    # Change to project root (parent of scripts directory)
    project_root = Path(__file__).parent.parent
    os.chdir(project_root)
    
    # Check prerequisites
    if not check_prerequisites():
        print("\n❌ Prerequisites not met. Please install missing tools.")
        return 1
    
    # Clean environment
    clean_environment()
    
    # Build package
    if not build_package():
        print("\n❌ Package build failed")
        return 1
    
    # Validate package
    if not validate_package():
        print("\n❌ Package validation failed")
        return 1
    
    # Test installation
    if not test_installation():
        print("\n❌ Installation test failed")
        return 1
    
    # Optional: Deploy to Test PyPI
    deploy_to_test_pypi()
    
    # Summary
    print("\n" + "="*60)
    print("✅ TEST DEPLOYMENT COMPLETED SUCCESSFULLY!")
    print("\nGenerated files:")
    if os.path.exists("dist"):
        for file in sorted(os.listdir("dist")):
            file_path = Path("dist") / file
            size = file_path.stat().st_size
            print(f"  - {file} ({size:,} bytes)")
    
    print("\nNext steps:")
    print("1. Test the package thoroughly")
    print("2. If everything works, deploy to production PyPI")
    print("3. Create GitHub release with the built artifacts")
    print("="*60)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())