# VARIANT - Railway Deployment Guide

This guide will help you deploy your VARIANT web application to Railway, making it publicly accessible with NYU branding and all features.

## 🚀 Quick Deployment Steps

### **Step 1: Prepare Your Repository**

Your repository is already prepared with the necessary files:
- ✅ `Procfile` - Tells Railway how to run your app
- ✅ `runtime.txt` - Specifies Python version
- ✅ `railway.json` - Railway configuration
- ✅ `requirements_web.txt` - Python dependencies
- ✅ `web_app.py` - Updated for Railway environment

### **Step 2: Deploy to Railway**

1. **Go to Railway**: Visit [railway.app](https://railway.app)
2. **Sign up/Login**: Use your GitHub account
3. **Create New Project**: Click "New Project"
4. **Deploy from GitHub**: Select "Deploy from GitHub repo"
5. **Select Repository**: Choose your VARIANT repository
6. **Wait for Deployment**: Railway will automatically build and deploy

### **Step 3: Get Your URL**

After deployment, Railway will provide:
- **Production URL**: `https://variant.up.railway.app`

## 📋 Detailed Deployment Instructions

### **1. Railway Account Setup**

```bash
# No command line setup needed - everything is web-based
```

**Steps:**
1. Go to [railway.app](https://railway.app)
2. Click "Start a New Project"
3. Sign in with GitHub
4. Authorize Railway to access your repositories

### **2. Project Creation**

**In Railway Dashboard:**
1. Click "New Project"
2. Select "Deploy from GitHub repo"
3. Find and select your VARIANT repository
4. Click "Deploy Now"

### **3. Environment Configuration**

Railway will automatically detect your Python app and configure:
- **Build Command**: `pip install -r requirements_web.txt`
- **Start Command**: `python web_app.py`
- **Port**: Automatically set from `PORT` environment variable

### **4. Domain Configuration**

**Add Custom Domain (Optional):**
1. Go to your project settings
2. Click "Domains"
3. Add your custom domain (e.g., `variant.yourdomain.com`)
4. Configure DNS records as instructed

## 🔧 Configuration Files Explained

### **Procfile**
```
web: python web_app.py
```
- Tells Railway this is a web application
- Specifies the command to start your app

### **runtime.txt**
```
python-3.8.18
```
- Specifies Python version for Railway
- Ensures compatibility with your dependencies

### **railway.json**
```json
{
  "$schema": "https://railway.app/railway.schema.json",
  "build": {
    "builder": "NIXPACKS"
  },
  "deploy": {
    "startCommand": "python web_app.py",
    "healthcheckPath": "/",
    "healthcheckTimeout": 300,
    "restartPolicyType": "ON_FAILURE",
    "restartPolicyMaxRetries": 10
  }
}
```

**Configuration Details:**
- **builder**: Uses Nixpacks for automatic dependency detection
- **startCommand**: How to start your application
- **healthcheckPath**: Endpoint to check if app is healthy
- **healthcheckTimeout**: Time to wait for health check
- **restartPolicy**: Automatically restart on failures

## 🌐 Accessing Your Application

### **Default Railway URL**
```
https://variant-production-xxxx.up.railway.app
```

### **Available Endpoints**
- **Main Application**: `/` - NYU-themed VARIANT interface
- **API Documentation**: `/docs` - FastAPI auto-generated docs
- **Health Check**: `/` - Application status

### **Features Available**
- ✅ **Professional Design**: Clean blue theme throughout
- ✅ **File Upload**: Virus sequence uploads
- ✅ **Analysis**: Mutation analysis tools
- ✅ **ZIP Downloads**: Bulk file downloads
- ✅ **Results Management**: File organization
- ✅ **Job History**: Analysis tracking

## 📊 Monitoring and Management

### **Railway Dashboard Features**

**Overview Tab:**
- Real-time deployment status
- Resource usage (CPU, memory)
- Request logs
- Error tracking

**Deployments Tab:**
- Deployment history
- Rollback options
- Build logs
- Environment variables

**Settings Tab:**
- Domain management
- Environment variables
- Team collaboration
- Billing information

### **Logs and Debugging**

**View Logs:**
1. Go to your project in Railway dashboard
2. Click "Deployments" tab
3. Select latest deployment
4. View build and runtime logs

**Common Log Locations:**
- **Build Logs**: During container building
- **Runtime Logs**: Application execution
- **Error Logs**: Application errors and exceptions

## 🔒 Security and Environment Variables

### **Environment Variables**

Railway automatically provides:
- `PORT` - Port number (set by Railway)
- `RAILWAY_ENVIRONMENT` - Environment name
- `RAILWAY_PROJECT_ID` - Project identifier

### **Custom Environment Variables**

**Add in Railway Dashboard:**
1. Go to project settings
2. Click "Variables" tab
3. Add key-value pairs as needed

**Example Variables:**
```bash
DEBUG=false
LOG_LEVEL=INFO
MAX_FILE_SIZE=100MB
```

## 💰 Cost Management

### **Free Tier Limits**
- **500 hours/month** of runtime
- **Auto-sleep** when not in use
- **Unlimited** deployments
- **Custom domains** included

### **Usage Optimization**
- **Auto-sleep**: App sleeps after 15 minutes of inactivity
- **Wake-up**: Automatically wakes when accessed
- **Monitoring**: Track usage in dashboard

### **Upgrade Options**
- **$5/month**: Unlimited hours
- **$20/month**: Pro plan with more resources
- **Custom**: Enterprise plans available

## 🚨 Troubleshooting

### **Common Issues**

**1. Build Failures**
```bash
# Check requirements_web.txt exists
# Verify Python version compatibility
# Check for missing dependencies
```

**2. Runtime Errors**
```bash
# Check application logs
# Verify environment variables
# Test locally first
```

**3. Port Issues**
```bash
# Railway sets PORT automatically
# Don't hardcode port 8000
# Use os.environ.get("PORT", 8000)
```

### **Debug Commands**

**Local Testing:**
```bash
# Test locally before deploying
python web_app.py

# Check if all dependencies are installed
pip install -r requirements_web.txt
```

**Railway CLI (Optional):**
```bash
# Install Railway CLI
npm install -g @railway/cli

# Login to Railway
railway login

# Deploy from command line
railway up
```

## 🔄 Updates and Maintenance

### **Automatic Deployments**

**GitHub Integration:**
- Railway automatically deploys on git push
- No manual deployment needed
- Rollback available if issues occur

### **Manual Updates**

**Update Process:**
1. Make changes to your code
2. Push to GitHub
3. Railway automatically redeploys
4. Monitor deployment status

### **Rollback Process**

**If Issues Occur:**
1. Go to Railway dashboard
2. Click "Deployments" tab
3. Select previous working deployment
4. Click "Rollback"

## 📈 Performance Optimization

### **Railway Optimizations**

**Automatic Features:**
- **CDN**: Global content delivery
- **SSL**: Automatic HTTPS certificates
- **Load Balancing**: Automatic traffic distribution
- **Auto-scaling**: Handles traffic spikes

### **Application Optimizations**

**Code Level:**
- Use async/await for I/O operations
- Implement proper caching
- Optimize database queries
- Minimize memory usage

## 🎯 Best Practices

### **Development Workflow**

1. **Local Testing**: Always test locally first
2. **GitHub Integration**: Use git for version control
3. **Environment Variables**: Don't hardcode sensitive data
4. **Logging**: Implement proper logging
5. **Error Handling**: Graceful error handling

### **Deployment Checklist**

- [ ] All dependencies in `requirements_web.txt`
- [ ] `Procfile` correctly configured
- [ ] `runtime.txt` specifies correct Python version
- [ ] Environment variables set
- [ ] Application handles `PORT` environment variable
- [ ] Health check endpoint implemented
- [ ] Error handling in place

## 🎉 Success!

Once deployed, your VARIANT application will be:
- 🌐 **Publicly accessible** via Railway URL
- 🎨 **Professionally designed** with clean blue theme
- 📁 **Fully functional** with all features
- 🔒 **Secure** with HTTPS
- 📊 **Monitored** with Railway dashboard
- 💰 **Free** within usage limits

**Your VARIANT application is now live and ready for researchers worldwide!** 🚀
