# VARIANT Web Application

## Overview

VARIANT (Viral mutAtion trackeR aImed At GeNome and proTein-level) is now available as a web-based application, providing an intuitive interface for virus mutation analysis without requiring command-line expertise.

## Features

- 🌐 **Web-based Interface**: Modern, responsive web interface accessible from any browser
- 🔬 **Multi-Virus Support**: Analyze SARS-CoV-2, HIV-1, H3N2, Chikungunya, and Zaire Ebola
- 📊 **Real-time Analysis**: Background job processing with real-time status updates
- 📁 **File Upload**: Upload custom MSA files for analysis
- 📈 **Results Visualization**: Interactive plots and downloadable results
- 🔄 **Job History**: Track and manage previous analysis jobs
- 🐳 **Docker Deployment**: Easy containerized deployment

## Quick Start

### Option 1: Docker Deployment (Recommended)

1. **Clone the repository**:
   ```bash
   git clone <repository-url>
   cd VARIANT
   ```

2. **Run the deployment script**:
   ```bash
   chmod +x deploy.sh
   ./deploy.sh
   ```

3. **Access the application**:
   - Web Interface: http://localhost:8000
   - API Documentation: http://localhost:8000/docs

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

## Usage

### Web Interface

1. **Select Virus Type**: Choose from available viruses (SARS-CoV-2, HIV-1, H3N2, etc.)
2. **Configure Analysis**:
   - Enter Genome ID (optional)
   - Upload MSA file (optional)
   - Enable/disable options (process all genomes, detect frameshifts)
3. **Start Analysis**: Click "Start Analysis" and monitor progress
4. **Download Results**: Download analysis results as ZIP file

### API Endpoints

| Method | Endpoint | Description |
|--------|----------|-------------|
| GET | `/` | Web interface |
| GET | `/api/viruses` | List available viruses |
| POST | `/api/analyze` | Start analysis job |
| GET | `/api/job/{job_id}` | Get job status |
| POST | `/api/upload-msa` | Upload MSA file |
| GET | `/api/jobs` | List all jobs |
| GET | `/api/results/{job_id}/download` | Download results |

### Example API Usage

```python
import requests

# List available viruses
response = requests.get('http://localhost:8000/api/viruses')
viruses = response.json()

# Start analysis
analysis_request = {
    "virus_name": "SARS-CoV-2",
    "genome_id": "EPI_ISL_123456",
    "process_all": False,
    "detect_frameshifts": True
}

response = requests.post('http://localhost:8000/api/analyze', json=analysis_request)
job = response.json()

# Check job status
job_id = job['job_id']
response = requests.get(f'http://localhost:8000/api/job/{job_id}')
status = response.json()
```

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

### Docker Deployment

```bash
# Build and run with Docker Compose
docker-compose up -d

# Or build manually
docker build -t variant-web .
docker run -p 8000:8000 variant-web
```

### Production Deployment

1. **Set up reverse proxy** (nginx):
   ```bash
   docker-compose --profile production up -d
   ```

2. **Configure SSL certificates**:
   - Place certificates in `./ssl/`
   - Update nginx configuration

3. **Set up monitoring**:
   - Application includes health checks
   - Monitor logs: `docker-compose logs -f`

## Architecture

### Components

- **FastAPI Backend**: RESTful API with automatic documentation
- **Jinja2 Templates**: Server-side templating
- **Background Tasks**: Asynchronous job processing
- **File Management**: Secure file upload and result storage
- **Job Queue**: In-memory job management (can be extended to Redis/Celery)

### Directory Structure

```
VARIANT/
├── web_app.py              # FastAPI application
├── templates/
│   └── index.html          # Web interface template
├── static/                 # Static files (CSS, JS, images)
├── uploads/                # Uploaded MSA files
├── results/                # Analysis results
├── data/                   # Virus data files
├── result/                 # Original analysis outputs
├── virus_config.yaml       # Virus configuration
├── requirements_web.txt    # Web dependencies
├── Dockerfile             # Docker configuration
├── docker-compose.yml     # Docker Compose setup
└── deploy.sh              # Deployment script
```

## Troubleshooting

### Common Issues

1. **Port already in use**:
   ```bash
   # Change port in docker-compose.yml or web_app.py
   ports:
     - "8001:8000"  # Use port 8001 instead
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

### Testing

```bash
# Run tests
pytest tests/

# API testing
curl -X POST http://localhost:8000/api/analyze \
  -H "Content-Type: application/json" \
  -d '{"virus_name": "SARS-CoV-2"}'
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
