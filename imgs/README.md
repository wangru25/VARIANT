# VARIANT Images Directory

This directory contains all generated visualization files from VARIANT analysis.

## Structure

```
imgs/
├── visualizations/          # Generated genome organization plots
│   ├── SARS-CoV-2/         # SARS-CoV-2 specific visualizations
│   │   ├── EPI_ISL_*.html  # Interactive HTML visualizations
│   │   └── EPI_ISL_*.pdf   # Static PDF visualizations
│   ├── HIV-1/              # HIV-1 specific visualizations
│   │   ├── MW881698.1_*.html # Interactive HTML visualizations
│   │   └── MW881698.1_*.pdf  # Static PDF visualizations
│   ├── Chikungunya/        # Chikungunya specific visualizations
│   ├── ZaireEbola/         # ZaireEbola specific visualizations
│   └── ...                 # Other virus folders
└── README.md               # This file
```

## File Types

- **HTML files**: Interactive visualizations with hover details and zoom capabilities
- **PDF files**: Static, publication-ready figures for scientific papers

## Naming Convention

Files are named using the pattern: `{Virus}_{GenomeID}_genome_organization.{extension}`

Examples:
- `SARS-CoV-2_EPI_ISL_16127650_genome_organization.html`
- `HIV-1_MW881698.1_genome_organization.pdf`

## Usage

Generated files can be:
- Opened in web browsers (HTML files)
- Included in scientific papers (PDF files)
- Used for presentations and reports
- Shared with collaborators

## Note

This directory is gitignored to avoid committing large generated files to the repository.
