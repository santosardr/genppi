# Scripts de Desenvolvimento

Esta pasta contém scripts para automatizar tarefas de desenvolvimento, build e deployment.

## Scripts Disponíveis

### Build e Teste
- `quick_check.py` - Verificação rápida de dependências e funcionalidades
- `build_and_test.py` - Build completo com testes abrangentes
- `test_installation.py` - Teste de instalação em ambiente limpo

### Deployment
- `deploy_test.py` - Deployment para TestPyPI (ambiente de testes)
- `deploy_production.py` - Deployment para PyPI oficial (produção)
- `test_upload.py` - Teste de upload para TestPyPI

## Como Usar

```bash
# Verificação rápida
python scripts/quick_check.py

# Build e teste completo
python scripts/build_and_test.py

# Deploy para testes
python scripts/deploy_test.py

# Deploy para produção
python scripts/deploy_production.py
```

## Pré-requisitos

- Python 3.6+
- build, twine, pytest instalados
- Credenciais configuradas em ~/.pypirc