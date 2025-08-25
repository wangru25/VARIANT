#!/bin/bash
# -*- coding: utf-8 -*-
# Author: Rui Wang
# Date: 2025-01-15 17:30:00
# LastModifiedBy: Rui Wang
# LastEditTime: 2025-01-15 17:30:00
# Email: rw3594@nyu.edu
# FilePath: /VARIANT/deploy_direct.sh
# Description: Direct deployment script for VARIANT web application without containerization.

set -e

echo "🚀 Starting VARIANT Web Application Direct Deployment..."

# Check if Python is available
if ! command -v python3 &> /dev/null; then
    echo "❌ Python 3 is not installed. Please install Python 3.8+ first."
    exit 1
fi

echo "✅ Python found: $(python3 --version)"

# Check if pip is available
if ! command -v pip3 &> /dev/null; then
    echo "❌ pip3 is not installed. Please install pip3 first."
    exit 1
fi

echo "✅ pip3 found: $(pip3 --version)"

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

# Create virtual environment
echo "🐍 Creating Python virtual environment..."
if [ ! -d "venv" ]; then
    python3 -m venv venv
fi

# Activate virtual environment
echo "🔧 Activating virtual environment..."
source venv/bin/activate

# Upgrade pip
echo "⬆️  Upgrading pip..."
pip install --upgrade pip

# Install Python dependencies
echo "📦 Installing Python dependencies..."
if [ -f "requirements_web.txt" ]; then
    pip install -r requirements_web.txt
else
    echo "⚠️  requirements_web.txt not found. Installing basic dependencies..."
    pip install fastapi uvicorn jinja2 python-multipart aiofiles
fi

# Create startup script
echo "📝 Creating startup script..."
cat > start_variant_direct.sh << 'EOF'
#!/bin/bash
# Start VARIANT web application directly

# Set environment variables
export PORT=${PORT:-8000}
export HOST=${HOST:-0.0.0.0}

# Create necessary directories if they don't exist
mkdir -p data result results uploads logs

# Activate virtual environment
source venv/bin/activate

# Start the application
echo "🚀 Starting VARIANT web application..."
echo "📊 Application will be available at: http://$HOST:$PORT"
echo "📚 API Documentation: http://$HOST:$PORT/docs"

# Run the application
python web_app.py
EOF

chmod +x start_variant_direct.sh

# Create systemd service file (if running as a service)
echo "🔧 Creating systemd service file..."
cat > variant-direct.service << EOF
[Unit]
Description=VARIANT Web Application (Direct)
After=network.target

[Service]
Type=simple
User=$USER
WorkingDirectory=$(pwd)
Environment=PORT=8000
Environment=HOST=0.0.0.0
ExecStart=$(pwd)/start_variant_direct.sh
Restart=always
RestartSec=10

[Install]
WantedBy=multi-user.target
EOF

echo "✅ VARIANT direct deployment completed!"
echo ""
echo "🌐 To start the application:"
echo "   ./start_variant_direct.sh"
echo ""
echo "🔧 To run as a service (requires sudo):"
echo "   sudo cp variant-direct.service /etc/systemd/system/"
echo "   sudo systemctl enable variant-direct"
echo "   sudo systemctl start variant-direct"
echo ""
echo "📊 To check status:"
echo "   sudo systemctl status variant-direct"
echo ""
echo "📝 To view logs:"
echo "   sudo journalctl -u variant-direct -f"
echo ""
echo "🛑 To stop the service:"
echo "   sudo systemctl stop variant-direct"
echo ""
echo "🧹 To clean up virtual environment:"
echo "   rm -rf venv"
