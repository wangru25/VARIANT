#!/bin/bash
# -*- coding: utf-8 -*-
# Author: Rui Wang
# Date: 2025-01-15 17:30:00
# LastModifiedBy: Rui Wang
# LastEditTime: 2025-01-15 17:30:00
# Email: rw3594@nyu.edu
# FilePath: /VARIANT/deploy_singularity.sh
# Description: Singularity deployment script for VARIANT web application on university servers.

set -e

echo "🚀 Starting VARIANT Web Application Deployment with Singularity..."

# Check if Singularity is installed
if ! command -v singularity &> /dev/null; then
    echo "❌ Singularity is not installed. Please install Singularity first."
    exit 1
fi

echo "✅ Singularity found: $(singularity --version)"

# Create necessary directories
echo "📁 Creating necessary directories..."
mkdir -p static templates uploads results data result plot logs

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

# Build Singularity container
echo "🔨 Building Singularity container..."
if [ -f "variant.sif" ]; then
    echo "⚠️  Existing variant.sif found. Removing..."
    rm -f variant.sif
fi

singularity build variant.sif Singularity.def

echo "✅ Singularity container built successfully!"

# Create startup script
echo "📝 Creating startup script..."
cat > start_variant.sh << 'EOF'
#!/bin/bash
# Start VARIANT web application with Singularity

# Set environment variables
export PORT=${PORT:-8000}
export HOST=${HOST:-0.0.0.0}

# Create necessary directories if they don't exist
mkdir -p data result results uploads logs

# Start the application
echo "🚀 Starting VARIANT web application..."
echo "📊 Application will be available at: http://$HOST:$PORT"
echo "📚 API Documentation: http://$HOST:$PORT/docs"

# Run the Singularity container
singularity run --bind ./data:/app/data,./result:/app/result,./results:/app/results,./uploads:/app/uploads,./logs:/app/logs variant.sif
EOF

chmod +x start_variant.sh

# Create systemd service file (if running as a service)
echo "🔧 Creating systemd service file..."
cat > variant.service << EOF
[Unit]
Description=VARIANT Web Application
After=network.target

[Service]
Type=simple
User=$USER
WorkingDirectory=$(pwd)
Environment=PORT=8000
Environment=HOST=0.0.0.0
ExecStart=$(pwd)/start_variant.sh
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
EOF

echo "✅ VARIANT Singularity deployment completed!"
echo ""
echo "🌐 To start the application:"
echo "   ./start_variant.sh"
echo ""
echo "🔧 To run as a service (requires sudo):"
echo "   sudo cp variant.service /etc/systemd/system/"
echo "   sudo systemctl enable variant"
echo "   sudo systemctl start variant"
echo ""
echo "📊 To check status:"
echo "   sudo systemctl status variant"
echo ""
echo "📝 To view logs:"
echo "   sudo journalctl -u variant -f"
echo ""
echo "🛑 To stop the service:"
echo "   sudo systemctl stop variant"
