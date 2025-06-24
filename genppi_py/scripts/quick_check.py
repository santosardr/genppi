#!/usr/bin/env python3
"""
Quick check script for GenPPI Python interface
Run from project root: python scripts/quick_check.py
"""

import sys
import importlib.util

def check_import(module_name, description):
    """Check if a module can be imported"""
    try:
        spec = importlib.util.find_spec(module_name)
        if spec is not None:
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            print(f"✅ {description}")
            return True
        else:
            print(f"❌ {description} - Module not found")
            return False
    except Exception as e:
        print(f"❌ {description} - Error: {e}")
        return False

def main():
    """Main check function"""
    print("GenPPI Python Interface - Quick Check")
    print("=" * 50)
    
    success = True
    
    # Check core dependencies
    success &= check_import("py7zr", "py7zr library")
    success &= check_import("multivolumefile", "multivolumefile library")
    
    # Check package modules
    success &= check_import("genppi_py", "GenPPI package")
    success &= check_import("genppi_py.genppi", "GenPPI main module")
    success &= check_import("genppi_py.download_model", "Download model module")
    success &= check_import("genppi_py.download_samples", "Download samples module")
    
    # Check specific functions
    try:
        from genppi_py.genppi import main as genppi_main
        print("✅ GenPPI main function")
    except Exception as e:
        print(f"❌ GenPPI main function - Error: {e}")
        success = False
    
    try:
        from genppi_py.download_model import download_model_file
        print("✅ Download model function")
    except Exception as e:
        print(f"❌ Download model function - Error: {e}")
        success = False
    
    # Check version
    try:
        import genppi_py
        print(f"✅ Package version: {genppi_py.__version__}")
    except Exception as e:
        print(f"❌ Package version - Error: {e}")
        success = False
    
    print("\n" + "=" * 50)
    if success:
        print("✅ ALL CHECKS PASSED!")
        print("Package is ready for deployment.")
    else:
        print("❌ SOME CHECKS FAILED!")
        print("Please fix the issues before deployment.")
    print("=" * 50)
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())