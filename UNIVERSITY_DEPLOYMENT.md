# VARIANT - University Server Deployment Guide

This guide provides deployment options for hosting VARIANT on university servers and infrastructure.

## 🏫 University Server Options

### **Option 1: University Server with Singularity (Recommended)**

Most university HPC environments support Singularity instead of Docker.

#### **Prerequisites**
- SSH access to university server
- Singularity support (usually available on HPC clusters)
- Domain/subdomain allocation (e.g., `variant.your-university.edu`)

#### **Deployment Steps**
```bash
# 1. Connect to university server
ssh your-username@university-server.edu

# 2. Check Singularity availability
singularity --version

# 3. Clone VARIANT repository
git clone https://github.com/wangru25/VARIANT.git
cd VARIANT

# 4. Deploy using Singularity script
./deploy_singularity.sh

# 5. Start the application
./start_variant.sh
```

#### **Singularity Container Management**
```bash
# Build the container
singularity build variant.sif Singularity.def

# Run the container
singularity run --bind ./data:/app/data,./result:/app/result,./results:/app/results,./uploads:/app/uploads,./logs:/app/logs variant.sif

# Shell into container (for debugging)
singularity shell variant.sif

# Execute commands in container
singularity exec variant.sif python web_app.py
```

#### **Running as a Service**
```bash
# Install as systemd service (if you have sudo access)
sudo cp variant.service /etc/systemd/system/
sudo systemctl enable variant
sudo systemctl start variant

# Check status
sudo systemctl status variant

# View logs
sudo journalctl -u variant -f
```

### **Option 2: University VPS/Dedicated Server (Docker Support)**

Most universities provide VPS or dedicated server access for research projects.

#### **Prerequisites**
- SSH access to university server
- Docker support (usually available)
- Domain/subdomain allocation (e.g., `variant.your-university.edu`)

#### **Deployment Steps**
```bash
# 1. Connect to university server
ssh your-username@university-server.edu

# 2. Install Docker (if not already installed)
sudo apt update
sudo apt install docker.io docker-compose

# 3. Clone VARIANT repository
git clone https://github.com/wangru25/VARIANT.git
cd VARIANT

# 4. Deploy using existing script
./deploy.sh

# 5. Configure nginx reverse proxy (if needed)
sudo apt install nginx
sudo nano /etc/nginx/sites-available/variant
```

#### **Nginx Configuration for University Domain**
```nginx
server {
    listen 80;
    server_name variant.your-university.edu;  # Your allocated domain

    location / {
        proxy_pass http://localhost:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        
        # Increase timeouts for long-running analyses
        proxy_connect_timeout 300s;
        proxy_send_timeout 300s;
        proxy_read_timeout 300s;
    }
}
```

#### **SSL Certificate (Let's Encrypt)**
```bash
# Install certbot
sudo apt install certbot python3-certbot-nginx

# Get SSL certificate
sudo certbot --nginx -d variant.your-university.edu

# Auto-renewal
sudo crontab -e
# Add: 0 12 * * * /usr/bin/certbot renew --quiet
```

### **Option 3: University HPC Cluster**

For high-performance computing clusters with job scheduling systems.

#### **Prerequisites**
- HPC cluster access (SLURM, PBS, etc.)
- Module system for Python/Singularity
- Web server access

#### **SLURM Job Script for Singularity**
```bash
#!/bin/bash
#SBATCH --job-name=variant-web
#SBATCH --output=variant-web-%j.out
#SBATCH --error=variant-web-%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

# Load required modules
module load singularity
module load python/3.8

# Start VARIANT web application with Singularity
cd /path/to/VARIANT
./start_variant.sh
```

#### **Interactive Session**
```bash
# Request interactive session
srun --pty --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=8G bash

# Load modules and start application
module load singularity
cd /path/to/VARIANT
./start_variant.sh
```

### **Option 4: University Cloud Infrastructure**

Many universities have cloud platforms (AWS, Azure, Google Cloud).

#### **AWS (via University Account)**
```bash
# Deploy to AWS using university credentials
aws configure  # Use university AWS credentials

# Create EC2 instance
aws ec2 run-instances \
  --image-id ami-0c02fb55956c7d316 \
  --count 1 \
  --instance-type t3.medium \
  --key-name your-key-pair \
  --security-group-ids sg-xxxxxxxxx

# Deploy application
ssh -i your-key.pem ubuntu@your-instance-ip
git clone https://github.com/wangru25/VARIANT.git
cd VARIANT
./deploy_singularity.sh
```

### **Option 5: University Shared Hosting (cPanel)**

For universities with shared hosting environments.

#### **Prerequisites**
- cPanel access
- Python support enabled
- SSH access (optional)

#### **Deployment Steps**
1. **Upload files via cPanel File Manager**
2. **Create Python application**
3. **Configure WSGI file**

#### **WSGI Configuration**
```python
# passenger_wsgi.py
import sys
import os

# Add application directory to path
sys.path.insert(0, '/home/username/public_html/variant')

# Set environment variables
os.environ['PYTHONPATH'] = '/home/username/public_html/variant'

# Import and run application
from web_app import app

application = app
```

## 🔧 University-Specific Configuration

### **Security Considerations**

#### **University Firewall Configuration**
```bash
# Request firewall rules from IT department
# Port 80 (HTTP)
# Port 443 (HTTPS)
# Port 8000 (Application, if direct access needed)
```

#### **Authentication Integration**
```python
# Add university SSO (if available)
# Example: Shibboleth, CAS, or LDAP integration
from flask_oidc import OpenIDConnect

oidc = OpenIDConnect(app)
```

### **Resource Management**

#### **Memory and CPU Limits for Singularity**
```bash
# Set resource limits when running Singularity
singularity run --bind ./data:/app/data \
  --bind ./result:/app/result \
  --bind ./results:/app/results \
  --bind ./uploads:/app/uploads \
  --bind ./logs:/app/logs \
  variant.sif
```

#### **Storage Configuration**
```bash
# Mount university storage
singularity run --bind /university/storage/variant:/app/data variant.sif
```

### **Monitoring and Logging**

#### **University Logging Standards**
```python
# Add university-compliant logging
import logging
from logging.handlers import RotatingFileHandler

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]',
    handlers=[
        RotatingFileHandler('logs/variant.log', maxBytes=10000000, backupCount=10),
        logging.StreamHandler()
    ]
)
```

## 📋 University IT Requirements Checklist

### **Before Deployment**
- [ ] **Server Access**: SSH or web-based access
- [ ] **Domain/Subdomain**: Allocated by IT department
- [ ] **SSL Certificate**: HTTPS support
- [ ] **Firewall Rules**: Port access approved
- [ ] **Resource Allocation**: CPU, memory, storage
- [ ] **Backup Strategy**: Data retention policies
- [ ] **Monitoring**: System health checks
- [ ] **Documentation**: User guides and technical docs

### **Security Requirements**
- [ ] **Authentication**: University SSO integration
- [ ] **Authorization**: Role-based access control
- [ ] **Data Protection**: FERPA/HIPAA compliance (if applicable)
- [ ] **Audit Logging**: User activity tracking
- [ ] **Backup**: Regular data backups
- [ ] **Updates**: Security patch management

## 🚀 Quick Start for University Deployment

### **Step 1: Contact IT Department**
```bash
# Prepare request with:
# - Application description
# - Resource requirements
# - Security considerations
# - Expected user base
# - Timeline for deployment
```

### **Step 2: Prepare Application**
```bash
# Ensure application is ready for production
git clone https://github.com/wangru25/VARIANT.git
cd VARIANT

# Test locally first (if you have Singularity locally)
./deploy_singularity.sh
./start_variant.sh
curl http://localhost:8000
```

### **Step 3: Deploy to University Server**
```bash
# Follow the appropriate deployment method above
# Based on your university's infrastructure
```

### **Step 4: Configure Domain and SSL**
```bash
# Set up domain routing
# Configure SSL certificate
# Test HTTPS access
```

## 📞 University IT Support

### **Common Contact Points**
- **Research Computing**: For HPC and high-performance resources
- **IT Services**: For general server access and domain management
- **Information Security**: For security reviews and compliance
- **Academic Computing**: For research project support

### **Required Information for IT Request**
1. **Project Description**: Virus mutation analysis tool
2. **Technical Requirements**: Python 3.8+, Singularity, 4GB RAM, 2 CPU cores
3. **User Base**: Researchers, students, collaborators
4. **Data Handling**: No sensitive data, public virus sequences only
5. **Security**: Standard web application security measures
6. **Timeline**: Immediate deployment needed

## 💰 Cost Considerations

### **University Resources (Usually Free)**
- ✅ **Server Access**: Often free for academic projects
- ✅ **Domain/SSL**: Usually provided by IT department
- ✅ **Storage**: University storage systems
- ✅ **Backup**: Institutional backup services
- ✅ **Support**: IT department assistance

### **Potential Costs**
- ❌ **Additional Resources**: If exceeding standard allocations
- ❌ **Custom Domain**: If not using university subdomain
- ❌ **Premium Support**: For urgent issues outside business hours

## 🔄 Maintenance and Updates

### **Regular Maintenance**
```bash
# Weekly: Check logs and performance
# Monthly: Update dependencies
# Quarterly: Security review
# Annually: Full system audit
```

### **Update Process**
```bash
# Pull latest changes
git pull origin main

# Rebuild Singularity container
singularity build --force variant.sif Singularity.def

# Restart application
./start_variant.sh

# Test functionality
curl https://variant.your-university.edu/docs
```

This setup provides a robust, secure, and cost-effective deployment for your university environment!
