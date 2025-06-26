#!/usr/bin/env python3
"""
Simplified build script for GenPPI Python interface
This script provides easy access to common development tasks.
"""

import sys
import subprocess
from pathlib import Path

def run_script(script_name, description):
    """Run a script from the scripts directory"""
    script_path = Path("scripts") / script_name
    if not script_path.exists():
        print(f"❌ Script not found: {script_path}")
        return False
    
    print(f"\n🚀 {description}")
    print("=" * 50)
    
    try:
        result = subprocess.run([sys.executable, str(script_path)], check=True)
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print(f"❌ Script failed with exit code {e.returncode}")
        return False

def main():
    """Main function"""
    if len(sys.argv) < 2:
        print("GenPPI Python Interface - Build Script")
        print("=" * 50)
        print("Usage: python dev.py <command>")
        print("")
        print("Available commands:")
        print("  check       - Quick dependency and functionality check")
        print("  test        - Full build and test")
        print("  deploy-test - Deploy to TestPyPI")
        print("  deploy-prod - Deploy to PyPI (production)")
        print("  install-test- Test installation")
        print("")
        print("Examples:")
        print("  python dev.py check")
        print("  python dev.py test")
        print("  python dev.py deploy-test")
        return 1
    
    command = sys.argv[1].lower()
    
    if command == "check":
        return 0 if run_script("quick_check.py", "Quick Check") else 1
    elif command == "test":
        return 0 if run_script("build_and_test.py", "Build and Test") else 1
    elif command == "deploy-test":
        return 0 if run_script("deploy_test.py", "Deploy to TestPyPI") else 1
    elif command == "deploy-prod":
        return 0 if run_script("deploy_production.py", "Deploy to Production") else 1
    elif command == "install-test":
        return 0 if run_script("test_installation.py", "Test Installation") else 1
    else:
        print(f"❌ Unknown command: {command}")
        print("Run 'python dev.py' for help")
        return 1

if __name__ == "__main__":
    sys.exit(main())