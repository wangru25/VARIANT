# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-23 11:00:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-24 18:24:41
Email: rw3594@nyu.edu
FilePath: /VARIANT/web_app.py
Description: FastAPI web application for VARIANT virus mutation analysis tool with dynamic virus support.
'''

import os
import sys
import tempfile
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime
import uuid
import zipfile
import json

from fastapi import FastAPI, File, UploadFile, HTTPException, BackgroundTasks, Form, Request
from fastapi.responses import HTMLResponse, FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from pydantic import BaseModel
import uvicorn
import yaml
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from jinja2 import Environment, FileSystemLoader

# Add src to path for imports
sys.path.append(str(Path(__file__).parent / "src"))

try:
    from src.core.mutation_processor import MutationProcessor
    from scripts.setup_virus_dataset import VirusDatasetSetup
except ImportError as e:
    print(f"Warning: Could not import core modules: {e}")

# Set Plotly template
pio.templates.default = "simple_white"

# Initialize FastAPI app
app = FastAPI(
    title="VARIANT - Viral mutAtion trackeR aImed At GeNome and proTein-level",
    description="A comprehensive web-based tool for virus mutation analysis",
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc"
)

# Create necessary directories
os.makedirs("static", exist_ok=True)
os.makedirs("templates", exist_ok=True)
os.makedirs("uploads", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("data", exist_ok=True)
os.makedirs("result", exist_ok=True)
os.makedirs("user_sessions", exist_ok=True)

# Mount static files
app.mount("/static", StaticFiles(directory="static"), name="static")

# Templates
templates = Jinja2Templates(directory="templates")

# Pydantic models
class AnalysisRequest(BaseModel):
    virus_name: str
    genome_id: Optional[str] = ""
    msa_file: Optional[str] = ""
    process_all: bool = False
    detect_frameshifts: bool = False
    segment: Optional[str] = None
    reference_genome: Optional[str] = None
    proteome_file: Optional[str] = None
    is_custom_virus: bool = False
    session_id: Optional[str] = None

class AnalysisResult(BaseModel):
    job_id: str
    status: str
    virus_name: str
    created_at: str
    results_url: Optional[str] = None
    error_message: Optional[str] = None

class DataUploadRequest(BaseModel):
    virus_name: str
    file_type: str  # 'reference_genome', 'proteome', 'msa'
    segment: Optional[str] = None
    is_custom_virus: bool = False
    session_id: Optional[str] = None

class CustomVirusRequest(BaseModel):
    virus_name: str
    is_multi_segment: bool = False
    segments: Optional[List[str]] = None

# Global storage for analysis jobs and user sessions
analysis_jobs: Dict[str, Dict[str, Any]] = {}
user_sessions: Dict[str, Dict[str, Any]] = {}

def load_virus_config():
    '''Load virus configuration from YAML file.'''
    try:
        with open("virus_config.yaml", "r") as f:
            return yaml.safe_load(f)
    except Exception as e:
        print(f"Warning: Could not load virus config: {e}")
        return {"viruses": {}}

def save_virus_config(config):
    '''Save virus configuration to YAML file.'''
    try:
        with open("virus_config.yaml", "w") as f:
            yaml.dump(config, f, default_flow_style=False)
        return True
    except Exception as e:
        print(f"Error saving virus config: {e}")
        return False

def create_custom_virus_config(virus_name: str, is_multi_segment: bool = False, segments: Optional[List[str]] = None):
    '''Create a new virus configuration for custom viruses.'''
    config = load_virus_config()
    
    if virus_name in config.get("viruses", {}):
        raise ValueError(f"Virus '{virus_name}' already exists in configuration")
    
    virus_config = {
        "codon_table_id": 1,
        "description": f"Custom virus: {virus_name}",
    }
    
    if is_multi_segment and segments:
        virus_config["segments"] = {}
        for segment in segments:
            virus_config["segments"][segment] = {
                "default_msa_file": "",
                "proteome_file": "",
                "reference_genome": ""
            }
    else:
        virus_config.update({
            "default_msa_file": "",
            "proteome_file": "",
            "reference_genome": ""
        })
    
    if "viruses" not in config:
        config["viruses"] = {}
    
    config["viruses"][virus_name] = virus_config
    
    if save_virus_config(config):
        return virus_config
    else:
        raise Exception("Failed to save virus configuration")

def get_session_id(request: Request) -> str:
    '''Get or create a session ID for the user.'''
    session_id = request.cookies.get("session_id")
    if not session_id:
        session_id = str(uuid.uuid4())
    return session_id

@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    '''Render the main homepage.'''
    session_id = get_session_id(request)
    response = templates.TemplateResponse("index.html", {"request": request})
    
    # Set session cookie if not already set
    if not request.cookies.get("session_id"):
        response.set_cookie(key="session_id", value=session_id, httponly=True, max_age=3600*24*7)  # 7 days
    
    return response

@app.get("/api/viruses")
async def get_available_viruses(request: Request):
    '''Get list of available viruses from configuration.'''
    session_id = get_session_id(request)
    
    try:
        config = load_virus_config()
        
        viruses = []
        for virus_name, virus_config in config.get("viruses", {}).items():
            virus_info = {
                "name": virus_name,
                "description": virus_config.get("description", ""),
                "is_multi_segment": "segments" in virus_config,
                "segments": list(virus_config.get("segments", {}).keys()) if "segments" in virus_config else [],
                "reference_genome": virus_config.get("reference_genome", ""),
                "proteome_file": virus_config.get("proteome_file", ""),
                "default_msa_file": virus_config.get("default_msa_file", ""),
                "is_custom": virus_config.get("description", "").startswith("Custom virus:")
            }
            viruses.append(virus_info)
        
        return {"viruses": viruses}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error loading virus configuration: {str(e)}")

@app.post("/api/create-custom-virus")
async def create_custom_virus(custom_virus: CustomVirusRequest, request: Request):
    '''Create a new custom virus configuration.'''
    session_id = get_session_id(request)
    
    try:
        virus_config = create_custom_virus_config(
            custom_virus.virus_name,
            custom_virus.is_multi_segment,
            custom_virus.segments
        )
        
        return {
            "success": True,
            "virus_name": custom_virus.virus_name,
            "config": virus_config,
            "message": f"Custom virus '{custom_virus.virus_name}' created successfully"
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

@app.post("/api/upload-data")
async def upload_data_file(
    virus_name: str = Form(...),
    file_type: str = Form(...),
    segment: Optional[str] = Form(None),
    is_custom_virus: bool = Form(False),
    file: UploadFile = File(...),
    request: Request = None
):
    '''Upload reference genome, proteome, or MSA files.'''
    session_id = get_session_id(request) if request else "default"
    
    allowed_extensions = {
        'reference_genome': ['.fasta', '.fa', '.fna'],
        'proteome': ['.fasta', '.fa', '.faa'],
        'msa': ['.txt', '.fasta', '.fa', '.aln', '.clustal']
    }
    
    if file_type not in allowed_extensions:
        raise HTTPException(status_code=400, detail=f"Invalid file type: {file_type}")
    
    file_ext = Path(file.filename).suffix.lower()
    if file_ext not in allowed_extensions[file_type]:
        raise HTTPException(
            status_code=400, 
            detail=f"Invalid file extension for {file_type}. Allowed: {', '.join(allowed_extensions[file_type])}"
        )
    
    # Create virus-specific directory structure
    if segment:
        target_dir = f"data/{virus_name}/{segment}"
        if file_type in ['reference_genome', 'proteome']:
            target_dir += "/refs"
        elif file_type == 'msa':
            target_dir += "/clustalW"
    else:
        target_dir = f"data/{virus_name}"
        if file_type in ['reference_genome', 'proteome']:
            target_dir += "/refs"
        elif file_type == 'msa':
            target_dir += "/clustalW"
    
    os.makedirs(target_dir, exist_ok=True)
    
    # Save file with user-friendly naming
    original_filename = file.filename
    base_name = Path(original_filename).stem
    extension = Path(original_filename).suffix
    
    # Create a descriptive filename based on file type
    if file_type == "reference_genome":
        filename = f"{virus_name}_reference_genome{extension}"
    elif file_type == "proteome":
        filename = f"{virus_name}_proteome{extension}"
    elif file_type == "msa":
        filename = f"{virus_name}_msa{extension}"
    else:
        filename = original_filename
    
    # If file already exists, add a number suffix
    counter = 1
    original_filename = filename
    while os.path.exists(os.path.join(target_dir, filename)):
        name_without_ext = Path(original_filename).stem
        ext = Path(original_filename).suffix
        filename = f"{name_without_ext}_{counter}{ext}"
        counter += 1
    
    filepath = os.path.join(target_dir, filename)
    
    try:
        with open(filepath, "wb") as buffer:
            shutil.copyfileobj(file.file, buffer)
        
        # Update virus configuration with file path
        if is_custom_virus:
            config = load_virus_config()
            if virus_name in config.get("viruses", {}):
                virus_config = config["viruses"][virus_name]
                if segment and "segments" in virus_config:
                    if segment in virus_config["segments"]:
                        if file_type == "reference_genome":
                            virus_config["segments"][segment]["reference_genome"] = filename
                        elif file_type == "proteome":
                            virus_config["segments"][segment]["proteome_file"] = filename
                        elif file_type == "msa":
                            virus_config["segments"][segment]["default_msa_file"] = filename
                else:
                    if file_type == "reference_genome":
                        virus_config["reference_genome"] = filename
                    elif file_type == "proteome":
                        virus_config["proteome_file"] = filename
                    elif file_type == "msa":
                        virus_config["default_msa_file"] = filename
                
                save_virus_config(config)
        
        return {
            "success": True,
            "filename": filename,
            "filepath": filepath,
            "size": os.path.getsize(filepath),
            "virus_name": virus_name,
            "file_type": file_type,
            "segment": segment,
            "is_custom_virus": is_custom_virus
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error uploading file: {str(e)}")

@app.post("/api/analyze", response_model=AnalysisResult)
async def start_analysis(analysis_request: AnalysisRequest, background_tasks: BackgroundTasks, request: Request):
    '''Start a virus mutation analysis job.'''
    session_id = get_session_id(request)
    job_id = str(uuid.uuid4())
    
    # Create job record
    job_record = {
        "job_id": job_id,
        "status": "queued",
        "virus_name": analysis_request.virus_name,
        "created_at": datetime.now().isoformat(),
        "request": analysis_request.dict(),
        "results": None,
        "error": None,
        "session_id": session_id
    }
    
    analysis_jobs[job_id] = job_record
    
    # Add background task
    background_tasks.add_task(run_analysis_job, job_id, analysis_request)
    
    return AnalysisResult(
        job_id=job_id,
        status="queued",
        virus_name=analysis_request.virus_name,
        created_at=job_record["created_at"]
    )

@app.get("/api/job/{job_id}")
async def get_job_status(job_id: str, request: Request):
    '''Get the status and results of an analysis job for the current user session.'''
    session_id = get_session_id(request)
    
    if job_id not in analysis_jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = analysis_jobs[job_id]
    
    # Check if job belongs to current user session
    if job.get("session_id") != session_id:
        raise HTTPException(status_code=403, detail="Access denied")
    
    return {
        "job_id": job_id,
        "status": job["status"],
        "virus_name": job["virus_name"],
        "created_at": job["created_at"],
        "results": job.get("results"),
        "error": job.get("error")
    }

@app.get("/api/results/{job_id}/download")
async def download_results(job_id: str, request: Request):
    '''Download analysis results as a zip file for the current user session.'''
    session_id = get_session_id(request)
    
    if job_id not in analysis_jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job = analysis_jobs[job_id]
    
    # Check if job belongs to current user session
    if job.get("session_id") != session_id:
        raise HTTPException(status_code=403, detail="Access denied")
    
    if job["status"] != "completed":
        raise HTTPException(status_code=400, detail="Analysis not completed")
    
    # Create zip file with results
    zip_path = f"results/{job_id}_results.zip"
    
    with zipfile.ZipFile(zip_path, 'w') as zipf:
        if job.get("results"):
            # Add result files
            for result_file in job["results"].get("files", []):
                if os.path.exists(result_file):
                    zipf.write(result_file, f"results/{os.path.basename(result_file)}")
            
            # Add data files if they exist
            virus_name = job["virus_name"]
            data_dir = f"data/{virus_name}"
            if os.path.exists(data_dir):
                for root, dirs, files in os.walk(data_dir):
                    for file in files:
                        file_path = os.path.join(root, file)
                        arc_name = f"data/{os.path.relpath(file_path, 'data')}"
                        zipf.write(file_path, arc_name)
    
    return FileResponse(
        zip_path,
        media_type="application/zip",
        filename=f"variant_results_{job_id}.zip"
    )

@app.get("/api/data-files/{virus_name}")
async def get_data_files(virus_name: str, request: Request, segment: Optional[str] = None):
    '''Get list of uploaded data files for a virus (session-based access).'''
    session_id = get_session_id(request)
    
    # For now, allow access to all virus data files
    # In a more secure implementation, you might want to track which viruses each user has uploaded
    # and only allow access to those specific viruses
    base_dir = f"data/{virus_name}"
    if segment:
        base_dir += f"/{segment}"
    
    if not os.path.exists(base_dir):
        return {"files": []}
    
    files = []
    for root, dirs, filenames in os.walk(base_dir):
        for filename in filenames:
            file_path = os.path.join(root, filename)
            rel_path = os.path.relpath(file_path, base_dir)
            file_type = "unknown"
            
            # Determine file type based on filename pattern (our new naming convention)
            if filename.startswith(f"{virus_name}_reference_genome"):
                file_type = "reference_genome"
            elif filename.startswith(f"{virus_name}_proteome"):
                file_type = "proteome"
            elif filename.startswith(f"{virus_name}_msa"):
                file_type = "msa"
            # Handle multi-segment virus naming patterns
            elif segment and filename.startswith(f"{segment}_reference_genome"):
                file_type = "reference_genome"
            elif segment and filename.startswith(f"{segment}_proteome"):
                file_type = "proteome"
            elif segment and filename.startswith(f"{segment}_msa"):
                file_type = "msa"
            else:
                # Fallback to directory-based detection for legacy files
                if "refs" in file_path:
                    if filename.endswith(('.fasta', '.fa', '.fna')):
                        # Check if it's a proteome file by looking for protein sequences
                        if filename.endswith(('.faa')) or 'proteome' in filename.lower():
                            file_type = "proteome"
                        else:
                            file_type = "reference_genome"
                    elif filename.endswith(('.faa')):
                        file_type = "proteome"
                elif "clustalW" in file_path:
                    file_type = "msa"
                else:
                    # Additional fallback for files in main virus directory
                    if filename.endswith(('.fasta', '.fa', '.fna')):
                        if 'proteome' in filename.lower():
                            file_type = "proteome"
                        elif any(keyword in filename.lower() for keyword in ['genome', 'ref', 'reference']):
                            file_type = "reference_genome"
                        else:
                            # Default to reference genome for .fasta files in main directory
                            file_type = "reference_genome"
                    elif filename.endswith(('.txt', '.aln', '.clustal')):
                        file_type = "msa"
            
            # Debug logging
            print(f"File: {filename}, Path: {file_path}, Type: {file_type}")
            
            files.append({
                "filename": filename,
                "path": rel_path,
                "type": file_type,
                "size": os.path.getsize(file_path),
                "modified": datetime.fromtimestamp(os.path.getmtime(file_path)).isoformat()
            })
    
    return {"files": files}

@app.get("/api/result-files/{virus_name}")
async def get_result_files(virus_name: str, request: Request, segment: Optional[str] = None):
    '''Get list of result files for a virus.'''
    base_dir = f"result/{virus_name}"
    if segment:
        base_dir += f"/{segment}"
    
    if not os.path.exists(base_dir):
        return {"files": []}
    
    files = []
    for root, dirs, filenames in os.walk(base_dir):
        for filename in filenames:
            if filename.endswith(('.txt', '.csv', '.html', '.pdf')):
                file_path = os.path.join(root, filename)
                rel_path = os.path.relpath(file_path, base_dir)
                
                files.append({
                    "filename": filename,
                    "path": rel_path,
                    "size": os.path.getsize(file_path),
                    "modified": datetime.fromtimestamp(os.path.getmtime(file_path)).isoformat()
                })
    
    return {"files": files}

@app.get("/api/download-file")
async def download_file(file_path: str, request: Request):
    '''Download a specific file from the data or result directories (session-based access).'''
    session_id = get_session_id(request)
    
    # For now, allow access to all files within allowed directories
    # In a more secure implementation, you might want to track which files each user has access to
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="File not found")
    
    # Security check: ensure file is within allowed directories
    allowed_dirs = ["data/", "result/", "uploads/", "results/"]
    if not any(file_path.startswith(allowed_dir) for allowed_dir in allowed_dirs):
        raise HTTPException(status_code=403, detail="Access denied")
    
    return FileResponse(
        file_path,
        filename=os.path.basename(file_path)
    )

@app.get("/api/download-virus-files/{virus_name}")
async def download_virus_files(virus_name: str, request: Request, type: str = "all", segment: Optional[str] = None):
    '''Download all files for a virus as a zip file.
    
    Args:
        virus_name: Name of the virus
        type: Type of files to download ('data', 'results', 'all')
        segment: Segment name for multi-segment viruses
    '''
    # Security check: ensure virus name is valid
    if not virus_name or ".." in virus_name or "/" in virus_name:
        raise HTTPException(status_code=400, detail="Invalid virus name")
    
    # Create zip file
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    zip_filename = f"{virus_name}_{type}_files_{timestamp}.zip"
    zip_path = f"results/{zip_filename}"
    
    with zipfile.ZipFile(zip_path, 'w') as zipf:
        if type in ["data", "all"]:
            # Add data files
            data_dir = f"data/{virus_name}"
            if segment:
                data_dir += f"/{segment}"
            
            if os.path.exists(data_dir):
                for root, dirs, files in os.walk(data_dir):
                    for file in files:
                        file_path = os.path.join(root, file)
                        arc_name = f"data/{os.path.relpath(file_path, 'data')}"
                        zipf.write(file_path, arc_name)
        
        if type in ["results", "all"]:
            # Add result files
            result_dir = f"result/{virus_name}"
            if segment:
                result_dir += f"/{segment}"
            
            if os.path.exists(result_dir):
                for root, dirs, files in os.walk(result_dir):
                    for file in files:
                        if file.endswith(('.txt', '.csv', '.html', '.pdf')):
                            file_path = os.path.join(root, file)
                            arc_name = f"results/{os.path.relpath(file_path, 'result')}"
                            zipf.write(file_path, arc_name)
    
    return FileResponse(
        zip_path,
        media_type="application/zip",
        filename=zip_filename
    )

async def run_analysis_job(job_id: str, analysis_request: AnalysisRequest):
    '''Run the actual analysis in the background.'''
    try:
        # Update job status
        analysis_jobs[job_id]["status"] = "running"
        
        # Determine MSA file path for custom viruses
        msa_file_path = analysis_request.msa_file
        if not msa_file_path and analysis_request.is_custom_virus:
            # For custom viruses, find the MSA file in the data directory
            virus_data_dir = f"data/{analysis_request.virus_name}"
            if analysis_request.segment:
                virus_data_dir += f"/{analysis_request.segment}"
            
            clustalw_dir = f"{virus_data_dir}/clustalW"
            if os.path.exists(clustalw_dir):
                # Look for MSA files
                for file in os.listdir(clustalw_dir):
                    if file.endswith(('.txt', '.fasta', '.fa', '.aln', '.clustal')):
                        msa_file_path = os.path.join(clustalw_dir, file)
                        break
            
            if not msa_file_path:
                raise Exception(f"No MSA file found for virus {analysis_request.virus_name}")
        
        # Initialize processor
        processor = MutationProcessor(analysis_request.virus_name, "virus_config.yaml")
        
        # Create temporary args object
        class Args:
            def __init__(self, **kwargs):
                for key, value in kwargs.items():
                    setattr(self, key, value)
        
        # Run main mutation analysis (always run this first)
        args_main = Args(
            virus=analysis_request.virus_name,
            genome_id=analysis_request.genome_id,
            msa_file=msa_file_path,
            process_all=analysis_request.process_all,
            detect_frameshifts=False,  # Disable frameshift detection for main analysis
            segment=analysis_request.segment,
            config="virus_config.yaml"
        )
        
        # Run main analysis
        processor.process(args_main)
        
        # Run frameshift analysis separately if requested
        if analysis_request.detect_frameshifts:
            args_frameshift = Args(
                virus=analysis_request.virus_name,
                genome_id=analysis_request.genome_id,
                msa_file=msa_file_path,
                process_all=analysis_request.process_all,
                detect_frameshifts=True,  # Enable frameshift detection
                segment=analysis_request.segment,
                config="virus_config.yaml"
            )
            
            # Run frameshift analysis
            processor.process(args_frameshift)
        
        # Collect results
        result_files = []
        virus_result_dir = f"result/{analysis_request.virus_name}"
        if os.path.exists(virus_result_dir):
            for root, dirs, files in os.walk(virus_result_dir):
                for file in files:
                    if file.endswith(('.csv', '.txt', '.html', '.pdf')):
                        result_files.append(os.path.join(root, file))
        
        # Update job with results
        analysis_jobs[job_id]["status"] = "completed"
        analysis_jobs[job_id]["results"] = {
            "files": result_files,
            "summary": f"Analysis completed for {analysis_request.virus_name}"
        }
        
    except Exception as e:
        analysis_jobs[job_id]["status"] = "failed"
        analysis_jobs[job_id]["error"] = str(e)

@app.get("/api/jobs")
async def list_jobs(request: Request):
    '''List analysis jobs for the current user session.'''
    session_id = get_session_id(request)
    user_jobs = [
        {
            "job_id": job_id,
            "status": job["status"],
            "virus_name": job["virus_name"],
            "created_at": job["created_at"]
        }
        for job_id, job in analysis_jobs.items()
        if job.get("session_id") == session_id
    ]
    return {"jobs": user_jobs}

if __name__ == "__main__":
    import os
    port = int(os.environ.get("PORT", 8000))
    uvicorn.run(app, host="0.0.0.0", port=port, reload=False)
