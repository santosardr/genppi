#!/usr/bin/env python3
"""
Test upload to TestPyPI
"""

import subprocess
import sys
import os
from pathlib import Path

def check_dist_files():
    """Check if distribution files exist"""
    dist_dir = Path("dist")
    if not dist_dir.exists():
        print("❌ dist/ directory not found. Run build first.")
        return False
    
    files = list(dist_dir.glob("*.whl")) + list(dist_dir.glob("*.tar.gz"))
    if not files:
        print("❌ No distribution files found in dist/")
        return False
    
    print("✅ Distribution files found:")
    for file in files:
        size = file.stat().st_size
        print(f"  - {file.name} ({size:,} bytes)")
    
    return True

def test_twine_config():
    """Test twine configuration"""
    try:
        result = subprocess.run([
            sys.executable, '-m', 'twine', 'check', 'dist/*'
        ], capture_output=True, text=True, check=True)
        print("✅ Twine check passed")
        return True
    except subprocess.CalledProcessError as e:
        print("❌ Twine check failed")
        print("STDERR:", e.stderr)
        return False

def upload_to_testpypi():
    """Upload to TestPyPI"""
    print("\nUploading to TestPyPI...")
    print("This will upload the package to https://test.pypi.org")
    
    try:
        result = subprocess.run([
            sys.executable, '-m', 'twine', 'upload',
            '--repository', 'testpypi',
            'dist/*'
        ], text=True, check=False)
        
        if result.returncode == 0:
            print("✅ Upload successful!")
            print("\nTo test installation:")
            print("pip install --index-url https://test.pypi.org/simple/ genppi-py")
            return True
        else:
            print("❌ Upload failed")
            return False
            
    except Exception as e:
        print(f"❌ Upload error: {e}")
        return False

def main():
    print("GenPPI - Test Upload to TestPyPI")
    print("=" * 50)
    
    # Change to project root (parent of scripts directory)
    project_root = Path(__file__).parent.parent
    os.chdir(project_root)
    
    # Check prerequisites
    if not check_dist_files():
        print("\nRun this first: python -m build")
        return 1
    
    if not test_twine_config():
        return 1
    
    # Confirm upload
    response = input("\nDo you want to upload to TestPyPI? (y/N): ").strip().lower()
    if response != 'y':
        print("Upload cancelled.")
        return 0
    
    if upload_to_testpypi():
        print("\n" + "="*50)
        print("✅ SUCCESS! Package uploaded to TestPyPI")
        print("\nNext steps:")
        print("1. Test installation:")
        print("   pip install --index-url https://test.pypi.org/simple/ genppi-py")
        print("2. Test functionality:")
        print("   genppi --help")
        print("3. If everything works, deploy to production PyPI")
        print("="*50)
        return 0
    else:
        print("\n" + "="*50)
        print("❌ Upload failed. Check the error messages above.")
        print("="*50)
        return 1

if __name__ == "__main__":
    sys.exit(main())