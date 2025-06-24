# Documentação do GenPPI

Esta pasta contém a documentação técnica do projeto GenPPI Python Interface.

## Arquivos

- `testing.tex` / `testing.pdf` - Guia de testes e validação
- `deployment.tex` / `deployment.pdf` - Guia de deployment e publicação

## Compilação dos PDFs

Para recompilar os PDFs a partir dos arquivos LaTeX:

```bash
cd docs/
pdflatex testing.tex
pdflatex deployment.tex
```

## Requisitos

- pdflatex
- Pacotes LaTeX: inputenc, babel, listings, xcolor, geometry