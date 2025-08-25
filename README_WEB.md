# VARIANT Web Application

## Overview

VARIANT (Viral mutAtion trackeR aImed At GeNome and proTein-level) is now available as a web-based application, providing an intuitive interface for virus mutation analysis without requiring command-line expertise. The application features a modern NYU-themed interface with streamlined file management and bulk download capabilities.

## Features

- 🌐 **Web-based Interface**: Modern, responsive web interface with NYU branding
- 🔬 **Multi-Virus Support**: Analyze SARS-CoV-2, HIV-1, H3N2, Chikungunya, and Zaire Ebola
- 📊 **Real-time Analysis**: Background job processing with real-time status updates
- 📁 **File Upload**: Upload custom MSA files for analysis
- 🗜️ **Bulk ZIP Downloads**: Download complete datasets as organized ZIP files
- 📈 **Results Visualization**: Interactive plots and downloadable results
- 🔄 **Job History**: Track and manage previous analysis jobs
- 🎨 **NYU Branding**: Professional interface using NYU's official color palette
- 🚀 **Railway Deployment**: Easy cloud deployment with free tier

## Quick Start

### Option 1: Railway Deployment (Recommended)

1. **Clone the repository**:
   ```bash
   git clone https://github.com/wangru25/VARIANT.git
   cd VARIANT
   ```

2. **Deploy to Railway**:
   - Go to [railway.app](https://railway.app)
   - Sign up with GitHub
   - Create new project → Deploy from GitHub repo
   - Select your VARIANT repository
   - Railway will auto-deploy your application

3. **Access the application**:
   - Web Interface: https://your-app-name.up.railway.app
   - API Documentation: https://your-app-name.up.railway.app/docs

📖 **[Complete Railway Deployment Guide](RAILWAY_DEPLOYMENT.md)**

### Option 2: Manual Installation

1. **Install dependencies**:
   ```bash
   pip install -r requirements_web.txt
   ```

2. **Create necessary directories**:
   ```bash
   mkdir -p static templates uploads results data result plot
   ```

3. **Run the application**:
   ```bash
   python web_app.py
   ```

## Interface Features

### **🎨 NYU-Themed Design**
- **Professional Interface**: NYU Violet and Ultra Violet color scheme
- **Clean Layout**: Simplified file listings without redundant information
- **Responsive Design**: Works on desktop, tablet, and mobile devices

### **📁 Enhanced File Management**
- **Bulk ZIP Downloads**: Download complete datasets as organized ZIP files
- **Upload Progress Tracking**: Visual indicators for file upload status
- **Clean File Display**: Simplified file listings showing only essential information
- **File Type Badges**: Clear visual indicators for different file types

## Usage Guide

### Web Interface

The VARIANT web application features a clean, NYU-themed interface with four main tabs:

#### 1. **Data Upload Tab**

##### **Option 1: Using Pre-configured Viruses**
1. **Select a pre-configured virus** from the dropdown (SARS-CoV-2, HIV-1, etc.)
2. **Upload your files** following the steps below

##### **Option 2: Creating Your Own Custom Virus (Recommended for New Users)**
1. **Create a Custom Virus:**
   - Enter a unique virus name (e.g., "MyCustomVirus")
   - Choose virus type:
     - **Single-segment virus**: Most viruses (SARS-CoV-2, HIV-1, etc.)
     - **Multi-segment virus**: Viruses like Influenza A (H3N2) with multiple segments
   - For multi-segment viruses, enter segment names (e.g., "segment_1,segment_2,segment_3")
   - Click "Create Custom Virus"
2. **Your custom virus will appear** in the virus selection dropdown
3. **Select your custom virus** and proceed with file uploads

##### **Required Files for Analysis**
For each virus, you need to upload **3 types of files**:

**1. Reference Genome**
- **Format**: FASTA (.fasta, .fa, .fna)
- **Content**: Complete genome sequence of the reference strain
- **Example**: `NC_045512.fasta` (SARS-CoV-2 reference genome)

**2. Reference Proteome**
- **Format**: FASTA (.fasta, .fa, .faa)
- **Content**: Protein sequences translated from the genome
- **Example**: `SARS-CoV-2_proteome.fasta`

**3. MSA File (Multiple Sequence Alignment)**
- **Format**: Text or FASTA (.txt, .fasta, .fa, .aln, .clustal)
- **Content**: Multiple sequence alignment of virus strains
- **Example**: `gisaid_hcov-19_20221201_20221230_China_0_msa.txt`

##### **File Upload Process**

**For Single-Segment Viruses:**
1. **Select your virus** from the dropdown
2. **Choose file type** (Reference Genome, Reference Proteome, or MSA)
3. **Select your file** from your computer
4. **Click "Upload File"**
5. **Repeat for all 3 file types**

**For Multi-Segment Viruses:**
1. **Select your virus** from the dropdown
2. **Choose the segment** you want to upload files for
3. **Upload files for each segment** (Reference Genome, Reference Proteome, MSA)
4. **Repeat for all segments**

##### **Upload Checklist**
The interface includes an interactive checklist that shows:
- ✅ **Reference Genome**: Uploaded
- ✅ **Reference Proteome**: Uploaded  
- ✅ **MSA File**: Uploaded

**Remember**: You must upload all 3 file types (reference genome, reference proteome, MSA) before you can start analysis!

##### **File Organization**
Your uploaded files are automatically organized with user-friendly names:
- **Single-segment viruses**: `data/VirusName/refs/` and `data/VirusName/clustalW/`
- **Multi-segment viruses**: `data/VirusName/SegmentName/refs/` and `data/VirusName/SegmentName/clustalW/`

**File Naming Convention:**
- **Reference Genome**: `VirusName_reference_genome.fasta`
- **Reference Proteome**: `VirusName_proteome.fasta` 
- **MSA File**: `VirusName_msa.txt`

If you upload multiple files of the same type, they'll be numbered (e.g., `VirusName_proteome_1.fasta`, `VirusName_proteome_2.fasta`).

#### 2. **Analysis Tab**
- **Virus Selection**: Choose from available viruses or custom viruses
- **Analysis Configuration**: Set genome ID, processing options, and frameshift detection
- **Real-time Monitoring**: Live status updates during analysis
- **Background Processing**: Non-blocking analysis with job queue

#### 3. **Files Tab**
- **Data Files Management**: View and manage uploaded data files
- **Result Files Management**: Access analysis results and reports
- **Bulk ZIP Downloads**: Download complete datasets as organized ZIP files
  - **Data Files ZIP**: Reference Genome, Reference Proteome, MSA files
  - **Result Files ZIP**: Analysis results, reports, and visualizations
- **Clean File Listings**: Simplified file display without redundant information

#### 4. **History Tab**
- **Job Tracking**: View all analysis jobs with status and timestamps
- **Job Management**: Monitor completed, running, and failed jobs
- **Status Indicators**: Visual status badges for easy identification

### File Download Options

#### **Individual ZIP Downloads**
- **"Download All Data Files (ZIP)"**: Downloads all data files for selected virus
- **"Download All Result Files (ZIP)"**: Downloads all analysis results

#### **File Organization in ZIP**
```
virus_name_data_files_20250115_143022.zip
├── data/
│   ├── virus_name/
│   │   ├── refs/
│   │   │   ├── reference_genome.fasta
│   │   │   └── proteome.fasta
│   │   └── clustalW/
│   │       └── msa_file.txt
└── results/
    └── virus_name/
        ├── mutation_summary.csv
        ├── hot_mutations.csv
        └── analysis_report.txt
```

### After Upload and Analysis

Once all files are uploaded:
1. **Go to the "Analysis" tab**
2. **Select your virus** (or custom virus) - the same one you uploaded files for
3. **Configure analysis parameters**:
   - **Process all genomes**: Recommended (processes all sequences in your MSA file)
   - **Detect frameshifts**: Optional (identifies potential frameshifting sites)
4. **Start the analysis** - the system will automatically use your uploaded files
5. **Monitor progress** in the "History" tab
6. **Download results** from the "Files" tab using the ZIP download options

**Note**: You don't need to upload MSA files again in the Analysis tab - the system uses the files you already uploaded in the Data Upload tab.

### Tips for Success

- **Use descriptive virus names** for custom viruses
- **Ensure file formats are correct** (FASTA for genomes/proteomes, text/FASTA for MSA)
- **Check the upload checklist** to ensure all required files are uploaded
- **Use the ZIP download options** for convenient bulk downloads
- **Monitor job progress** in the History tab for long-running analyses

## API Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/` | Web interface |
| GET | `/api/viruses` | List available viruses |
| POST | `/api/analyze` | Start analysis job |
| GET | `/api/job/{job_id}` | Get job status |
| POST | `/api/upload-data` | Upload data files |
| POST | `/api/create-custom-virus` | Create custom virus |
| GET | `/api/jobs` | List all jobs |
| GET | `/api/results/{job_id}/download` | Download job results |
| GET | `/api/download-virus-files/{virus_name}` | Download virus files as ZIP |
| GET | `/api/data-files/{virus_name}` | List data files for virus |
| GET | `/api/result-files/{virus_name}` | List result files for virus |

### Example API Usage

```python
import requests

# List available viruses
response = requests.get('http://localhost:8000/api/viruses')
viruses = response.json()

# Create custom virus
custom_virus = {
    "virus_name": "MyCustomVirus",
    "is_multi_segment": False,
    "segments": None
}
response = requests.post('http://localhost:8000/api/create-custom-virus', json=custom_virus)

# Start analysis
analysis_request = {
    "virus_name": "SARS-CoV-2",
    "genome_id": "EPI_ISL_123456",
    "process_all": False,
    "detect_frameshifts": True
}

response = requests.post('http://localhost:8000/api/analyze', json=analysis_request)
job = response.json()

# Download virus files as ZIP
response = requests.get('http://localhost:8000/api/download-virus-files/SARS-CoV-2?type=data')
# Saves as: SARS-CoV-2_data_files_20250115_143022.zip
```

## NYU Branding

The web application features NYU's official color palette:

### **Primary Colors**
- **NYU Violet** (`#57068c`): Primary branding and headers
- **Ultra Violet** (`#8900e1`): Secondary elements and gradients
- **NYU Black** (`#000000`): Text and dark elements

### **Secondary Colors**
- **Deep Violet** (`#330662`): Hover states and darker elements
- **Light Violet 2** (`#eee6f3`): Background gradients
- **NYU Teal** (`#009b8a`): Success states and positive actions
- **NYU Magenta** (`#fb0f78`): Error/danger states
- **NYU Yellow** (`#f4ec51`): Warning states and file indicators

### **Design Features**
- **Gradient Headers**: NYU Violet to Ultra Violet gradients
- **Professional Cards**: Clean card layouts with NYU color schemes
- **Consistent Buttons**: NYU-themed button styles with hover effects
- **Accessible Design**: High contrast and readable typography

## Configuration

### Virus Configuration

The application uses `virus_config.yaml` to configure available viruses. Example:

```yaml
viruses:
  SARS-CoV-2:
    codon_table_id: 1
    description: SARS-CoV-2 coronavirus
    proteome_file: SARS-CoV-2_proteome.fasta
    reference_genome: NC_045512.fasta
    default_msa_file: SARS-CoV-2_msa.txt
```

### Environment Variables

| Variable | Description | Default |
|----------|-------------|---------|
| `PYTHONPATH` | Python path | `/app` |
| `PYTHONUNBUFFERED` | Python output buffering | `1` |

## Deployment Options

### Local Development

```bash
python web_app.py
```

### Railway Deployment

```bash
# Deploy to Railway (Recommended)
# 1. Go to railway.app
# 2. Connect your GitHub repository
# 3. Railway will auto-deploy your application

# Access your application
# Web Interface: https://your-app-name.up.railway.app
# API Documentation: https://your-app-name.up.railway.app/docs
```

### Production Features

1. **Automatic SSL**: Railway provides HTTPS certificates
2. **Global CDN**: Fast loading worldwide
3. **Auto-scaling**: Handles traffic spikes automatically
4. **Monitoring**: Built-in logs and metrics dashboard

📖 **[Complete Railway Deployment Guide](RAILWAY_DEPLOYMENT.md)**

## Architecture

### Components

- **FastAPI Backend**: RESTful API with automatic documentation
- **Jinja2 Templates**: Server-side templating with NYU branding
- **Background Tasks**: Asynchronous job processing
- **File Management**: Secure file upload and result storage
- **ZIP Generation**: Dynamic ZIP file creation for bulk downloads
- **Job Queue**: In-memory job management (can be extended to Redis/Celery)

### Directory Structure

```
VARIANT/
├── web_app.py              # FastAPI application
├── templates/
│   └── index.html          # NYU-themed web interface template
├── static/                 # Static files (CSS, JS, images)
├── uploads/                # Uploaded MSA files
├── results/                # Analysis results and ZIP files
├── data/                   # Virus data files
├── result/                 # Original analysis outputs
├── virus_config.yaml       # Virus configuration
├── requirements_web.txt    # Web dependencies
├── Procfile               # Railway deployment configuration
├── runtime.txt            # Python version specification
├── railway.json           # Railway configuration
└── RAILWAY_DEPLOYMENT.md  # Railway deployment guide
```

## Troubleshooting

### Common Issues

1. **Port already in use**:
   ```bash
   # Railway automatically handles port configuration
   # No manual port configuration needed
   ```

2. **Permission errors**:
   ```bash
   # Fix file permissions
   chmod +x deploy.sh
   chmod -R 755 uploads results
   ```

3. **Import errors**:
   ```bash
   # Ensure PYTHONPATH is set
   export PYTHONPATH=/path/to/VARIANT
   ```

4. **ZIP download issues**:
   ```bash
   # Check results directory permissions
   chmod -R 755 results/
   # Ensure sufficient disk space for ZIP generation
   ```

### Logs and Debugging

```bash
# View application logs
docker-compose logs -f variant-web

# Check application health
curl http://localhost:8000/docs

# Debug mode
python web_app.py --debug
```

## Development

### Adding New Features

1. **New API endpoints**: Add to `web_app.py`
2. **UI improvements**: Modify `templates/index.html`
3. **Virus support**: Update `virus_config.yaml`
4. **NYU branding**: Follow NYU color palette guidelines

### Testing

```bash
# Run tests
pytest tests/

# API testing
curl -X POST http://localhost:8000/api/analyze \
  -H "Content-Type: application/json" \
  -d '{"virus_name": "SARS-CoV-2"}'

# ZIP download testing
curl -O http://localhost:8000/api/download-virus-files/SARS-CoV-2?type=data
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make changes following the coding guidelines
4. Test thoroughly
5. Submit a pull request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Support

For issues and questions:
- Create an issue on GitHub
- Check the API documentation at `/docs`
- Review the troubleshooting section above
