# nextSSR BMC Bioinformatics Manuscript Directory

This directory contains the complete scientific manuscript describing **nextSSR** formatted for submission to *BMC Bioinformatics* (Software Article / Short Communication).

---

## 📂 Directory Contents

- `bmc_bioinformatics_nextssr.md`: Full manuscript in Markdown format.
- `bmc_nextssr.tex`: Complete LaTeX manuscript adhering to BMC formatting standards.
- `references.bib`: BibTeX file containing all scientific citations.

---

## 📄 How to Compile to PDF

### Option 1: Using `pdflatex` + `bibtex`
```bash
pdflatex bmc_nextssr.tex
bibtex bmc_nextssr
pdflatex bmc_nextssr.tex
pdflatex bmc_nextssr.tex
```

### Option 2: Using `pandoc`
```bash
pandoc bmc_bioinformatics_nextssr.md -o bmc_bioinformatics_nextssr.pdf --citeproc --bibliography=references.bib
```
