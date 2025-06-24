#!/usr/bin/env python3
"""
Script para deployment em produção (PyPI oficial)
"""

import subprocess
import sys
import os
from pathlib import Path

def run_command(cmd, description, check=True):
    """Run a command and handle errors"""
    print(f"\n{'='*50}")
    print(f"Executando: {description}")
    print(f"Comando: {' '.join(cmd)}")
    print(f"{'='*50}")
    
    try:
        result = subprocess.run(cmd, check=check, text=True)
        if result.returncode == 0:
            print("✅ SUCESSO")
        else:
            print("⚠️  WARNING: Non-zero exit code")
        return result.returncode == 0
    except subprocess.CalledProcessError as e:
        print("❌ FALHOU")
        return False
    except FileNotFoundError:
        print(f"❌ FALHOU: Comando não encontrado: {cmd[0]}")
        return False

def check_pypirc():
    """Check if .pypirc is configured for production PyPI"""
    pypirc_path = Path.home() / '.pypirc'
    
    if not pypirc_path.exists():
        print("❌ Arquivo .pypirc não encontrado")
        print("Configure suas credenciais do PyPI oficial primeiro")
        return False
    
    with open(pypirc_path, 'r') as f:
        content = f.read()
        if '[pypi]' in content:
            print("✅ Configuração PyPI oficial encontrada")
            return True
        else:
            print("⚠️  Configuração PyPI oficial não encontrada em .pypirc")
            print("Adicione a seção [pypi] com suas credenciais")
            return False

def deploy_to_pypi():
    """Deploy to production PyPI"""
    print("GenPPI - Deployment para PyPI Oficial (Produção)")
    print("=" * 60)
    
    # Change to project root (parent of scripts directory)
    project_root = Path(__file__).parent.parent
    os.chdir(project_root)
    
    # Check prerequisites
    if not check_pypirc():
        print("\nConfigure o .pypirc primeiro:")
        print("""
[distutils]
index-servers =
    pypi
    testpypi

[pypi]
repository = https://upload.pypi.org/legacy/
username = __token__
password = pypi-SEU_TOKEN_AQUI

[testpypi]
repository = https://test.pypi.org/legacy/
username = __token__
password = pypi-SEU_TOKEN_TESTPYPI_AQUI
""")
        return False
    
    # Check if dist files exist
    if not Path("dist").exists() or not list(Path("dist").glob("*.whl")):
        print("❌ Arquivos de distribuição não encontrados")
        print("Execute primeiro: python -m build")
        return False
    
    print("✅ Arquivos de distribuição encontrados:")
    for file in Path("dist").glob("*"):
        size = file.stat().st_size
        print(f"  - {file.name} ({size:,} bytes)")
    
    # Confirm deployment
    print("\n" + "="*60)
    print("⚠️  ATENÇÃO: Você está prestes a fazer deployment para o PyPI OFICIAL!")
    print("Isso tornará o pacote disponível para todos os usuários Python do mundo.")
    print("Comando que os usuários usarão: pip install genppi-py")
    print("="*60)
    
    response = input("\nTem certeza que quer continuar? Digite 'SIM' para confirmar: ").strip()
    if response != 'SIM':
        print("Deployment cancelado.")
        return False
    
    # Upload to PyPI
    success = run_command([
        sys.executable, "-m", "twine", "upload", "dist/*"
    ], "Upload para PyPI oficial")
    
    if success:
        print("\n" + "="*60)
        print("🎉 SUCESSO! Pacote publicado no PyPI oficial!")
        print("\nUsuários finais podem agora instalar com:")
        print("    pip install genppi-py")
        print("\nVerifique em: https://pypi.org/project/genppi-py/")
        print("="*60)
        return True
    else:
        print("\n" + "="*60)
        print("❌ Falha no upload. Verifique os erros acima.")
        print("="*60)
        return False

if __name__ == "__main__":
    success = deploy_to_pypi()
    sys.exit(0 if success else 1)