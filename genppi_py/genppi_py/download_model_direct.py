#!/usr/bin/env python3
"""
Download model.dat file directly (without 7zip extraction)
"""

import os
import sys
import urllib.request
from pathlib import Path
# NOVO: Import para localizar arquivos de forma robusta
from importlib import resources

# URL para download direto do model.dat
MODEL_DAT_URL = "https://github.com/santosardr/genppi/releases/download/model/model.dat"

# Função auxiliar para encontrar o diretório 'bin' de forma robusta
def _get_package_bin_dir():
    """Retorna o Path para o diretório 'bin' dentro do pacote."""
    try:
        # Método moderno e preferencial para pacotes instalados
        return resources.files('genppi_py') / 'bin'
    except (AttributeError, ModuleNotFoundError):
        # Fallback para ambientes de desenvolvimento
        return Path(__file__).parent.absolute() / 'bin'

def download_model_direct(bin_dir=None):
    """
    Download model.dat file directly from a pre-built URL
    
    Args:
        bin_dir (str or Path, optional): Directory to save the model.dat file.
            If None, use the bin directory in the package.
    
    Returns:
        bool: True if successful, False otherwise
    """
    if bin_dir is None:
        # MUDANÇA: Usar a nova função para encontrar o diretório bin do pacote
        bin_dir = _get_package_bin_dir()
    else:
        bin_dir = Path(bin_dir)
    
    model_path = bin_dir / "model.dat"
    
    # Pula se o model.dat já existe
    if model_path.exists():
        print("model.dat file already exists.")
        return True
    
    print(f"Downloading model.dat directly from {MODEL_DAT_URL}...")
    
    # Cria o diretório bin se não existir
    os.makedirs(bin_dir, exist_ok=True)
    
    try:
        # Baixa o arquivo model.dat
        urllib.request.urlretrieve(MODEL_DAT_URL, model_path)
        
        # Verifica se o model.dat foi baixado
        if model_path.exists() and os.path.getsize(model_path) > 0:
            print(f"model.dat successfully downloaded to {model_path}")
            return True
        else:
            print("model.dat was not downloaded properly.", file=sys.stderr)
            return False
    except Exception as e:
        print(f"Error downloading model.dat: {e}", file=sys.stderr)
        return False

def main():
    """
    Command-line entry point
    """
    import argparse
    
    parser = argparse.ArgumentParser(description="Download model.dat file directly (without 7zip extraction)")
    parser.add_argument("--output-dir", "-o", help="Directory to save the model.dat file")
    
    args = parser.parse_args()
    
    success = download_model_direct(args.output_dir)
    
    if success:
        print("Model file download completed successfully.")
        return 0
    else:
        print("Model file download failed.", file=sys.stderr)
        return 1

if __name__ == "__main__":
    sys.exit(main())