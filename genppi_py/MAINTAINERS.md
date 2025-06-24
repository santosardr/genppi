# Maintainer's Guide for GenPPI Python Interface

This document provides instructions for maintainers of the GenPPI Python interface package.

## Preparing a New Release

1. Update the version number in `genppi_py/__init__.py` and `setup.py`.

2. Prepare the release files:

   ```bash
   python prepare_release.py
   ```

   This will create the following files in the `release` directory:
   - `genppi_linux_executables.zip`
   - `genppi_mac_executables.zip`
   - `genppi_windows_executables.zip`
   - `model.dat`

3. Create a new release on GitHub:
   - Go to https://github.com/santosardr/genppi/releases/new
   - Tag version: `v1.0.0` (replace with the actual version)
   - Release title: `GenPPI v1.0.0` (replace with the actual version)
   - Upload the files from the `release` directory as assets

4. Update the download URLs in the code if necessary:
   - In `setup.py`, update the `DOWNLOAD_BASE_URL` if the release tag changes
   - In `genppi_py/genppi.py`, update the `download_base_url` in both functions if the release tag changes

5. Build and publish the package to PyPI:

   ```bash
   # Build the package
   python -m build

   # Upload to PyPI
   python -m twine upload dist/*
   ```

## Testing the Package

Before releasing, make sure to test the package:

1. Install the package in development mode:

   ```bash
   pip install -e .
   ```

2. Run the test script:

   ```bash
   python test_genppi.py
   ```

3. Test with actual protein files:

   ```bash
   genppi -dir /path/to/protein/files -genefusion -pp -cn
   ```

## Troubleshooting

If users report issues with downloading executables:

1. Check that the release files are correctly uploaded to GitHub
2. Verify that the download URLs in the code are correct
3. Check that the file names in the code match the actual file names in the release

If users report issues with the executables not working:

1. Make sure the executables have the correct permissions (executable bit set on Unix-like systems)
2. Check that the executables are compatible with the target platforms
3. Verify that the model.dat file is correctly downloaded when using the -ml option