#!/usr/bin/env python3
"""
GenPPI - A Python interface for GenPPI (Genomic-based Protein-Protein Interaction prediction)

This module provides a Python interface to the GenPPI executables compiled from SBCL.
"""

import os
import sys
import platform
import subprocess
import shutil
import tempfile
import urllib.request
from pathlib import Path
import argparse
from importlib import resources

def _get_bin_dir():
    """Retorna o Path para o diretório 'bin' dentro do pacote."""
    try:
        return resources.files('genppi_py') / 'bin'
    except (AttributeError, ModuleNotFoundError):
        return Path(__file__).parent.absolute() / 'bin'

def get_executable_path(use_db=False):
    """
    Get the path to the appropriate GenPPI executable.
    """
    system = platform.system()
    bin_dir = _get_bin_dir()
    
    if system == 'Linux':
        exe = 'genppidb32g-Linux' if use_db else 'genppi32g-Linux'
    elif system == 'Darwin':
        exe = 'genppidb32g-Mac.x' if use_db else 'genppi32g-Mac.x'
    elif system == 'Windows':
        exe = 'genppidb32g-Win.exe' if use_db else 'genppi32g-Win.exe'
    else:
        raise RuntimeError(f"Unsupported operating system: {system}")
    
    executable = bin_dir / exe
    
    if not executable.exists():
        try:
            download_missing_executable(executable, system, use_db)
        except Exception as e:
            raise FileNotFoundError(
                f"GenPPI executable not found and could not be downloaded automatically.\n"
                f"Please ensure you have internet access or manually download the executable."
            )
    
    return executable

def download_missing_executable(executable_path, system, use_db):
    """
    Download a missing executable file. (Esta função está correta e não precisa de mudanças)
    """
    import urllib.request
    github_raw_url = "https://raw.githubusercontent.com/santosardr/genppi/master/"
    if system == "Linux":
        exe_name = "genppidb32g-Linux" if use_db else "genppi32g-Linux"
    elif system == "Darwin":
        exe_name = "genppidb32g-Mac.x" if use_db else "genppi32g-Mac.x"
    elif system == "Windows":
        exe_name = "genppidb32g-Win.exe" if use_db else "genppi32g-Win.exe"
    else:
        raise RuntimeError(f"Unsupported platform: {system}")
    download_url = f"{github_raw_url}binaries/{exe_name}"
    os.makedirs(executable_path.parent, exist_ok=True)
    try:
        urllib.request.urlretrieve(download_url, executable_path)
        if system != "Windows":
            executable_path.chmod(0o755)
    except Exception as e:
        raise RuntimeError(f"Failed to download executable from {download_url}: {e}")

def check_and_download_model_file():
    """
    Check if model.dat exists and download if needed. (Esta função está correta)
    """
    from .download_model import download_model_file
    bin_dir = _get_bin_dir()
    model_path = bin_dir / 'model.dat'
    if model_path.exists():
        return
    success = download_model_file()
    if not success:
        print("Failed to download or extract model.dat.", file=sys.stderr)
        print("Please run 'genppi-download-model' to manually download the model.", file=sys.stderr)
        sys.exit(1)

def run_genppi(args, use_db=False):
    """
    Run the GenPPI executable with the given arguments. (Esta função está correta)
    """
    bin_dir = _get_bin_dir()
    executable = get_executable_path(use_db)
    
    if platform.system() != 'Windows':
        executable.chmod(0o755)
    
    if '-ml' in args:
        check_and_download_model_file()
    
    cmd = [str(executable)] + args
    
    try:
        process = subprocess.run(
            cmd,
            cwd=str(bin_dir),
            stderr=subprocess.PIPE,
            text=True,
            check=False
        )
        if process.returncode != 0 and process.stderr:
            print(process.stderr, end='', file=sys.stderr)
        return process.returncode
    except Exception as e:
        print(f"Error executing GenPPI: {e}", file=sys.stderr)
        return 1

def main():
    """
    Main entry point for the GenPPI command line interface.
    """
    parser = argparse.ArgumentParser(
        description='GenPPI - Genomic-based Protein-Protein Interaction prediction',
        add_help=False
    )
    
    parser.add_argument('-h', '--help', '-help', action='store_true', dest='show_help', help='Show this help message and exit')
    parser.add_argument('--use-db', action='store_true', help='Use the database version of GenPPI (genppidb)')
    parser.add_argument('-dir', help='Working directory containing the protein files')
    # ... (todos os outros argumentos estão aqui)
    parser.add_argument('-genefusion', action='store_true', help='Make PPI predictions by gene fusion')
    parser.add_argument('-gfp', type=float, help='Gene fusion score percentage (default: 5)')
    parser.add_argument('-pp', action='store_true', help='Enable phylogenetic profiles (Note: Method 1 runs by default with multiple genomes)')
    parser.add_argument('-ppp', type=float, help='Phylogenetic profiles score percentage (default: 30)')
    parser.add_argument('-ppdifftolerated', type=int, help='Phylogenetic profiles difference tolerated (default: 0)')
    parser.add_argument('-pphistofilter', action='store_true', help='Apply histogram filter to phylogenetic profiles')
    parser.add_argument('-cnp', type=float, help='Conserved neighborhood percentage (default: 65)')
    parser.add_argument('-expt', choices=['fixed', 'dynamic'], help='Expansion type (default: fixed)')
    parser.add_argument('-ws', type=int, help='Window size (default: 1)')
    parser.add_argument('-w1', type=int, help='Window size 1 (default: 10)')
    parser.add_argument('-cw1', type=int, help='Conservation window 1 (default: 4)')
    parser.add_argument('-w2', type=int, help='Window size 2 (default: 7)')
    parser.add_argument('-cw2', type=int, help='Conservation window 2 (default: 3)')
    parser.add_argument('-w3', type=int, help='Window size 3 (default: 5)')
    parser.add_argument('-cw3', type=int, help='Conservation window 3 (default: 2)')
    parser.add_argument('-aadifflimit', type=int, choices=[0, 1], help='Amino acid difference limit (default: 1)')
    parser.add_argument('-aacheckminlimit', type=int, help='Minimum amount of amino acids to check (default: 25)')
    parser.add_argument('-ml', action='store_true', help='Use machine learning model to predict protein similarity')
    parser.add_argument('-ppmethod', type=int, choices=range(1, 8), help='Phylogenetic profiles method (1-7): 1=No filters, 2=Only conserved neighborhood interactions, 3=Limit interactions, 4=Limit by weight, 5=Threshold filter, 6=Delete groups by weight, 7=Exclude unwanted profiles')
    parser.add_argument('-ppcomplete', action='store_true', help='Complete phylogenetic profiles')
    parser.add_argument('-pplimit', action='store_true', help='Limit phylogenetic profiles')
    parser.add_argument('-ppiterlimit', type=int, help='Phylogenetic profiles iteration limit')
    parser.add_argument('-pptrim', action='store_true', help='Trim phylogenetic profiles')
    parser.add_argument('-trim', type=int, help='Trim value')
    parser.add_argument('-ppthreshold', action='store_true', help='Apply threshold to phylogenetic profiles')
    parser.add_argument('-threshold', type=float, help='Threshold value')
    parser.add_argument('-plusminus', help='Plus/minus value')
    parser.add_argument('-ppdeletegroup', action='store_true', help='Delete group from phylogenetic profiles')
    parser.add_argument('-grouplimit', type=int, help='Group limit')
    parser.add_argument('-ppdeleteprofile', action='store_true', help='Delete profile from phylogenetic profiles')
    parser.add_argument('-profiles', help='Profiles to delete')
    
    args, unknown = parser.parse_known_args()
    
    if args.show_help:
        parser.print_help()
        return 0

    # <<< CORREÇÃO CRÍTICA AQUI >>>
    # Resolve o caminho do argumento -dir para um caminho absoluto.
    # Isso é essencial porque mudamos o CWD para a pasta 'bin'.
    # Sem isso, o executável Lisp não encontraria o diretório de dados.
    if args.dir:
        args.dir = os.path.abspath(args.dir)
    # <<< FIM DA CORREÇÃO >>>

    cmd_args = []
    
    # Handle ppmethod parameter mapping to individual Lisp parameters
    if args.ppmethod is not None:
        # Check for conflicts with individual method parameters
        method_params = ['ppcomplete', 'ppcn', 'pplimit', 'pptrim', 'ppthreshold', 'ppdeletegroup', 'ppdeleteprofile']
        if any(getattr(args, param, False) for param in method_params):
            print("Warning: -ppmethod parameter overrides individual method parameters", file=sys.stderr)
        
        if args.ppmethod == 1:
            cmd_args.extend(["-ppcomplete"])
        elif args.ppmethod == 2:
            cmd_args.extend(["-ppcn"])
        elif args.ppmethod == 3:
            cmd_args.extend(["-pplimit"])
            if args.ppiterlimit is not None:
                cmd_args.extend(["-ppiterlimit", str(args.ppiterlimit)])
            else:
                cmd_args.extend(["-ppiterlimit", "500000"])
        elif args.ppmethod == 4:
            cmd_args.extend(["-pptrim"])
            if args.trim is not None:
                cmd_args.extend(["-trim", str(args.trim)])
            else:
                cmd_args.extend(["-trim", "45000"])
        elif args.ppmethod == 5:
            cmd_args.extend(["-ppthreshold"])
            if args.threshold is not None and args.plusminus is not None:
                cmd_args.extend(["-threshold", str(args.threshold), "-plusminus", args.plusminus])
            else:
                print("Error: Method 5 requires -threshold and -plusminus parameters", file=sys.stderr)
                return 1
        elif args.ppmethod == 6:
            cmd_args.extend(["-ppdeletegroup"])
            if args.grouplimit is not None:
                cmd_args.extend(["-grouplimit", str(args.grouplimit)])
            else:
                cmd_args.extend(["-grouplimit", "45000"])
        elif args.ppmethod == 7:
            cmd_args.extend(["-ppdeleteprofile"])
            if args.profiles is not None:
                cmd_args.extend(["-profiles", args.profiles])
            else:
                print("Error: Method 7 requires -profiles parameter", file=sys.stderr)
                return 1
    
    for arg_name, arg_value in vars(args).items():
        if arg_name in ['use_db', 'show_help', 'ppmethod']:
            continue
        if arg_value is not None:
            if isinstance(arg_value, bool):
                if arg_value:
                    cmd_args.append(f"-{arg_name}")
            else:
                cmd_args.append(f"-{arg_name}")
                cmd_args.append(str(arg_value))
    
    cmd_args.extend(unknown)
    
    if not cmd_args and not unknown:
        parser.print_help()
        return 0
    
    return run_genppi(cmd_args, args.use_db)

if __name__ == '__main__':
    sys.exit(main())