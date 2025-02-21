## OSPEAD Milestone II README

This directory contains three items:

- `OSPEAD-milestone-II.pdf`, a PDF report.
- `OSPEAD-docs`, a website.
- `decoyanalysis_0.1.0.tar.gz`, an R package that implements the BJR and Patra-Sen estimator.

`OSPEAD-milestone-II.pdf` is a followup to the "Fully Specified Estimation Plan for OSPEAD" document I submitted for Milestone 1.

The `OSPEAD-docs` website can be accessed by opening `OSPEAD-docs/_book/index.html` in a browser. It is designed to be read by you, the reviewer, but also for members of the Monero community that prefer to read procedures and code more than the mathematical formulas and statistical methodologies that the PDFs have contained. Most of the statistical results of Milestone II are laid out in OSPEAD-docs.

`decoyanalysis_0.1.0.tar.gz` is an R package that is needed to reproduce the results in OSPEAD-docs. See the "Requirements" chapter in OSPEAD-docs for installation instructions.

Some large data files would be necessary for full reproducibility. The files contain the arrival times of transactions in nodes' txpools and isthmus' list of nonstandard transactions. I decided to place them in a separate archive file because their size exceeds 1GB. They can be downloaded from https://rucknium.me/html/OSPEAD-milestone-II-supplementary-data.tar.xz

## Rebuilding OSPEAD-docs

The OSPEAD-docs website is not fully self-contained. It fetches fonts and javascript files from external sources when loaded. A self-contained version can be created, but the HTML files are much larger.

To build a self-contained OSPEAD-docs, open the `_quarto.yml` file and find this section:
```yml
format:
  html:
```

Under it, add:
```yml
    embed-resources: true
    self-contained-math: true
```

Save the file.
Install Quarto: https://quarto.org/docs/get-started/
In the OSPEAD-docs directory, run this shell command:
```bash
quarto render --to html
```
