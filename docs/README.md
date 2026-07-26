# VARIANT Documentation

This folder contains all documentation files for the VARIANT (Viral mutAtion trackeR aImed At GeNome and proTein-level) package.

## 📚 Documentation Index

### 🚀 Getting Started
- **[README.md](../README.md)** - Main project README (in root directory)
- **[README_WEB.md](README_WEB.md)** - Web application user guide and interface documentation

### 🛠️ Deployment & Setup
- **[RAILWAY_DEPLOYMENT.md](RAILWAY_DEPLOYMENT.md)** - Railway.app deployment guide
- **[plot_visualization_guide.md](plot_visualization_guide.md)** - Plotting and visualization guide

### 📊 Analysis & Visualization
- **[COMPREHENSIVE_VISUALIZATION_SUMMARY.md](COMPREHENSIVE_VISUALIZATION_SUMMARY.md)** - Comprehensive visualization analysis summary
- **[PROTEIN_VISUALIZATION_SUMMARY.md](PROTEIN_VISUALIZATION_SUMMARY.md)** - Protein-level visualization analysis summary

### 🧬 PRF Detection
- **[src/core/prf/README.md](../src/core/prf/README.md)** - Programmed ribosomal frameshifting (+1/−1) detection module: algorithm, CLI flags, folding tools (RNAfold/PKNOTS/ProbKnot/NUPACK), and output columns
- **[data/tRNA_abundance_README.md](../data/tRNA_abundance_README.md)** - Format and provenance of the optional tRNA-abundance proxy table

## 📁 Project Structure

```
VARIANT/
├── README.md                    # Main project documentation
├── docs/                        # 📚 All documentation files
│   ├── README.md               # This documentation index
│   ├── README_WEB.md           # Web application guide
│   ├── RAILWAY_DEPLOYMENT.md   # Deployment guide
│   ├── plot_visualization_guide.md # Plotting guide
│   ├── COMPREHENSIVE_VISUALIZATION_SUMMARY.md
│   └── PROTEIN_VISUALIZATION_SUMMARY.md
├── src/                        # Source code
│   └── core/prf/               # PRF (+1/−1) detection module (see its README.md)
├── data/                       # Virus data
├── result/                     # Analysis results
├── imgs/                       # Generated visualizations (e.g. imgs/visualizations/)
└── templates/                  # Web templates
```

## 🔗 Quick Links

- **Main Project**: [README.md](../README.md)
- **Web Application**: [README_WEB.md](README_WEB.md)
- **Deployment**: [RAILWAY_DEPLOYMENT.md](RAILWAY_DEPLOYMENT.md)
- **Visualization**: [plot_visualization_guide.md](plot_visualization_guide.md)

---

*Last updated: 2026-07-25*
