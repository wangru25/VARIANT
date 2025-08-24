#!/bin/bash
# -*- coding: utf-8 -*-
# Author: Rui Wang
# Date: 2025-08-23 11:00:00
# LastModifiedBy: Rui Wang
# LastEditTime: 2025-08-23 11:00:00
# Email: rw3594@nyu.edu
# FilePath: /VARIANT/deploy.sh
# Description: Deployment script for VARIANT web application.

set -e

echo "🚀 Starting VARIANT Web Application Deployment..."

# Check if Docker is installed
if ! command -v docker &> /dev/null; then
    echo "❌ Docker is not installed. Please install Docker first."
    exit 1
fi

# Check if Docker Compose is installed
if ! command -v docker-compose &> /dev/null; then
    echo "❌ Docker Compose is not installed. Please install Docker Compose first."
    exit 1
fi

# Create necessary directories
echo "📁 Creating necessary directories..."
mkdir -p static templates uploads results data result plot

# Check if virus_config.yaml exists
if [ ! -f "virus_config.yaml" ]; then
    echo "⚠️  Warning: virus_config.yaml not found. Creating a basic template..."
    cat > virus_config.yaml << EOF
global:
  default_codon_table: 1
  supported_codon_tables:
    1: Standard genetic code

viruses:
  SARS-CoV-2:
    codon_table_id: 1
    description: SARS-CoV-2 coronavirus
    proteome_file: SARS-CoV-2_proteome.fasta
    reference_genome: NC_045512.fasta
    default_msa_file: SARS-CoV-2_msa.txt
EOF
fi

# Build and start the application
echo "🔨 Building Docker image..."
docker-compose build

echo "🚀 Starting VARIANT web application..."
docker-compose up -d

# Wait for the application to start
echo "⏳ Waiting for application to start..."
sleep 10

# Check if the application is running
if curl -f http://localhost:8000/docs > /dev/null 2>&1; then
    echo "✅ VARIANT web application is running successfully!"
    echo ""
    echo "🌐 Access the application at:"
    echo "   - Web Interface: http://localhost:8000"
    echo "   - API Documentation: http://localhost:8000/docs"
    echo "   - ReDoc Documentation: http://localhost:8000/redoc"
    echo ""
    echo "📋 Available endpoints:"
    echo "   - GET  /                    - Web interface"
    echo "   - GET  /api/viruses         - List available viruses"
    echo "   - POST /api/analyze         - Start analysis"
    echo "   - GET  /api/job/{job_id}    - Get job status"
    echo "   - POST /api/upload-msa      - Upload MSA file"
    echo "   - GET  /api/jobs            - List all jobs"
    echo ""
    echo "🛑 To stop the application, run: docker-compose down"
    echo "📊 To view logs, run: docker-compose logs -f"
else
    echo "❌ Application failed to start. Check logs with: docker-compose logs"
    exit 1
fi
