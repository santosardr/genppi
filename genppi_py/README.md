# GenPPI Python Interface

A Python interface for GenPPI (Genomic-based Protein-Protein Interaction prediction).

## Installation

### From PyPI (Recommended)
```bash
pip install genppi-py
```

### From TestPyPI (Development)
```bash
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ genppi-py
```

### From Source
```bash
git clone https://github.com/santosardr/genppi.git
cd genppi
pip install .
```

## Project Structure

```
genppi/
├── genppi_py/          # Main package code
├── docs/               # Documentation (LaTeX/PDF)
├── scripts/            # Development and deployment scripts
├── tests/              # Automated tests
├── tools/              # Development tools and configurations
├── README.md           # This file
├── setup.py            # Installation configuration
└── pyproject.toml      # Modern project configuration
```

### For Developers

- `docs/` - Technical documentation in LaTeX/PDF
- `scripts/` - Scripts for build, test and deployment
- `tests/` - Automated tests
- `tools/` - Development tools and configurations

## Usage

### Command Line Interface

After installation, you can use GenPPI from the command line:

```bash
# Basic usage
genppi.py -dir /path/to/protein/files

# With gene fusion
genppi.py -dir /path/to/protein/files -genefusion

# Using machine learning
genppi.py -dir /path/to/protein/files -ml
```

### Download Model and Samples

```bash
# Download machine learning model
genppi-download-model

# Download sample data
genppi-download-samples
```

### Python API

```python
from genppi_py.genppi import run_genppi

# Run GenPPI with arguments
run_genppi(['-dir', '/path/to/protein/files', '-genefusion'])
```

## How it Works

When you install the package, it will:

1. Download the appropriate pre-compiled GenPPI executables for your operating system
2. Download and extract the machine learning model file (model.dat) using py7zr
3. Install the Python wrapper that provides a command-line interface

The package uses py7zr library exclusively for extracting multivolume 7z archives (model.7z.001, model.7z.002, model.7z.003).

## Dependencies

- Python 3.6+
- py7zr >= 0.20.0 (for 7z extraction)
- multivolumefile >= 0.2.3 (for multivolume support)

All dependencies are installed automatically.

## Development

### Quick Setup

```bash
# Clone repository
git clone https://github.com/santosardr/genppi.git
cd genppi

# Install development dependencies
pip install -r tools/requirements-dev.txt

# Install in development mode
pip install -e .
```

### Development Scripts

```bash
# Quick check
python scripts/quick_check.py

# Build and test
python scripts/build_and_test.py

# Deploy to TestPyPI
python scripts/deploy_test.py

# Deploy to PyPI
python scripts/deploy_production.py
```

### Running Tests

```bash
# Run all tests
python -m pytest tests/

# Run specific test
python tests/test_genppi.py
```

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.