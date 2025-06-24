# Testes do GenPPI

Esta pasta contém os testes automatizados do projeto.

## Arquivos de Teste

- `test_genppi.py` - Testes unitários principais

## Executar Testes

```bash
# Executar todos os testes
python -m pytest tests/

# Executar teste específico
python tests/test_genppi.py

# Executar com cobertura
python -m pytest tests/ --cov=genppi_py
```

## Estrutura dos Testes

Os testes cobrem:
- Importação de módulos
- Funcionalidades básicas
- Validação de dependências
- Testes de integração