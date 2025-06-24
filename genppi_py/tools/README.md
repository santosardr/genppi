# Ferramentas de Desenvolvimento

Esta pasta contém ferramentas e configurações para desenvolvimento.

## Arquivos

- `requirements-dev.txt` - Dependências de desenvolvimento
- `tox.ini` - Configuração para testes em múltiplas versões Python
- `prepare_release.py` - Script para preparar releases

## Configuração do Ambiente de Desenvolvimento

```bash
# Instalar dependências de desenvolvimento
pip install -r tools/requirements-dev.txt

# Executar testes com tox
tox

# Preparar release
python tools/prepare_release.py
```

## Dependências de Desenvolvimento

- build - Para construir pacotes
- twine - Para upload ao PyPI
- pytest - Para testes
- tox - Para testes em múltiplas versões
- flake8 - Para análise de código
- black - Para formatação de código