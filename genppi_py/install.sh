#!/bin/bash
# Installation script for GenPPI Python interface

# Check Python version
python_version=$(python3 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
required_version="3.8"

if [ "$(printf '%s\n' "$required_version" "$python_version" | sort -V | head -n1)" != "$required_version" ]; then
    echo "Error: Python $required_version or higher is required."
    echo "Current Python version: $python_version"
    echo "Please install Python $required_version or higher from https://www.python.org/downloads/"
    exit 1
fi

echo "Python version $python_version is compatible."

# Check if pip is available
if ! command -v pip &> /dev/null; then
    echo "Error: pip is not installed. Please install pip first."
    exit 1
fi

# Install the package in development mode
pip install -e .

echo "GenPPI Python interface installed successfully!"
echo "You can now use the 'genppi' command from the command line."
echo "Note: Use 'genppi' command, not 'genppi.py'."