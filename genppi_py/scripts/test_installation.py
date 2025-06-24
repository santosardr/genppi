#!/usr/bin/env python3
"""
Script para testar a instalação do genppi-py do TestPyPI
"""

import subprocess
import sys
import tempfile
import os
from pathlib import Path

def run_command(cmd, description, check=True):
    """Run a command and handle errors"""
    print(f"\n{'='*50}")
    print(f"Executando: {description}")
    print(f"Comando: {' '.join(cmd)}")
    print(f"{'='*50}")
    
    try:
        result = subprocess.run(cmd, check=check, capture_output=True, text=True)
        if result.returncode == 0:
            print("✅ SUCESSO")
        else:
            print("⚠️  WARNING: Non-zero exit code")
        
        if result.stdout:
            print("STDOUT:", result.stdout)
        if result.stderr and result.returncode != 0:
            print("STDERR:", result.stderr)
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print("❌ FALHOU")
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
        return False
    except FileNotFoundError:
        print(f"❌ FALHOU: Comando não encontrado: {cmd[0]}")
        return False

def test_installation():
    """Test installation in a clean environment"""
    print("GenPPI - Teste de Instalação do TestPyPI")
    print("=" * 60)
    
    with tempfile.TemporaryDirectory() as temp_dir:
        venv_path = Path(temp_dir) / "test_venv"
        
        # Create virtual environment
        success = run_command([
            sys.executable, "-m", "venv", str(venv_path)
        ], "Criar ambiente virtual de teste")
        
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
        ], "Atualizar pip no ambiente de teste")
        
        # Install genppi-py from TestPyPI with dependencies from PyPI
        success &= run_command([
            str(venv_pip), "install", 
            "--index-url", "https://test.pypi.org/simple/",
            "--extra-index-url", "https://pypi.org/simple/",
            "genppi-py==0.1.5"
        ], "Instalar genppi-py==0.1.5 do TestPyPI")
        
        # Test basic functionality
        if success:
            success &= run_command([
                str(venv_python), "-c", 
                "import genppi_py; print(f'GenPPI versão: {genppi_py.__version__}')"
            ], "Testar importação do pacote")
            
            success &= run_command([
                str(venv_python), "-c", 
                "import py7zr; print(f'py7zr versão: {py7zr.__version__}')"
            ], "Verificar py7zr")
            
            success &= run_command([
                str(venv_python), "-c", 
                "import multivolumefile; print('multivolumefile OK')"
            ], "Verificar multivolumefile")
            
            # Test console scripts
            success &= run_command([
                str(venv_python), "-m", "genppi_py.genppi", "--help"
            ], "Testar execução do módulo genppi", check=False)
    
    return success

def main():
    """Main test function"""
    success = test_installation()
    
    print("\n" + "="*60)
    if success:
        print("✅ TESTE DE INSTALAÇÃO BEM-SUCEDIDO!")
        print("\nComando para instalação:")
        print("pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ genppi-py==0.1.5")
    else:
        print("❌ TESTE DE INSTALAÇÃO FALHOU!")
        print("\nVerifique os erros acima e tente novamente.")
    print("="*60)
    
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())