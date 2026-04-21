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
import time
from pathlib import Path
from typing import Dict, List, Optional, Any
from datetime import datetime
import uuid
import zipfile
import json

from fastapi import FastAPI, File, UploadFile, HTTPException, BackgroundTasks, Form, Request
from fastapi.responses import HTMLResponse, FileResponse, JSONResponse
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

# Add Existing-Dual-Search to path
_DUAL_SEARCH_DIR = str(Path(__file__).parent / "Existing-Dual-Search")
sys.path.insert(0, _DUAL_SEARCH_DIR)

try:
    from src.core.mutation_processor import MutationProcessor
    from scripts.setup_virus_dataset import VirusDatasetSetup
except ImportError as e:
    print(f"Warning: Could not import core modules: {e}")

try:
    from convert_dbn_to_dssr import convert_dbn_to_dssr
    from PDBto2D import process_id as pdbto2d_process_id
    from Dual_Library import run_dual_library
    from plot_dual_graph import plot_dual_graph
    _DUAL_SEARCH_AVAILABLE = True
except ImportError as e:
    print(f"Warning: Could not import Existing-Dual-Search modules: {e}")
    _DUAL_SEARCH_AVAILABLE = False

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
os.makedirs("templates", exist_ok=True)
os.makedirs("uploads", exist_ok=True)
os.makedirs("results", exist_ok=True)
os.makedirs("data", exist_ok=True)
os.makedirs("result", exist_ok=True)
os.makedirs("user_sessions", exist_ok=True)

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
    visualization_type: Optional[str] = None

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
session_virus_access: Dict[str, List[str]] = {}

_JOBS_FILE = "results/analysis_jobs.json"
_ACCESS_FILE = "results/session_virus_access.json"

def _load_jobs_from_disk() -> None:
    """Load persisted jobs from disk into the in-memory dict on startup."""
    global analysis_jobs
    if os.path.exists(_JOBS_FILE):
        try:
            with open(_JOBS_FILE, "r") as f:
                analysis_jobs = json.load(f)
            print(f"Loaded {len(analysis_jobs)} persisted jobs from {_JOBS_FILE}")
        except Exception as e:
            print(f"Warning: could not load persisted jobs: {e}")

def _save_jobs_to_disk() -> None:
    """Persist in-memory jobs to disk (best-effort, non-blocking)."""
    try:
        os.makedirs(os.path.dirname(_JOBS_FILE), exist_ok=True)
        # Write via a tmp file to avoid corruption on crash
        tmp = _JOBS_FILE + ".tmp"
        with open(tmp, "w") as f:
            json.dump(analysis_jobs, f, default=str)
        os.replace(tmp, _JOBS_FILE)
    except Exception as e:
        print(f"Warning: could not persist jobs: {e}")


def _load_access_from_disk() -> None:
    """Load persisted session virus access map."""
    global session_virus_access
    if os.path.exists(_ACCESS_FILE):
        try:
            with open(_ACCESS_FILE, "r") as f:
                raw = json.load(f)
            # Normalize to {session_id: [virus_name, ...]}
            if isinstance(raw, dict):
                session_virus_access = {
                    str(k): [str(v) for v in vals] if isinstance(vals, list) else []
                    for k, vals in raw.items()
                }
            print(f"Loaded {len(session_virus_access)} session access records from {_ACCESS_FILE}")
        except Exception as e:
            print(f"Warning: could not load session access records: {e}")


def _save_access_to_disk() -> None:
    """Persist session virus access map."""
    try:
        os.makedirs(os.path.dirname(_ACCESS_FILE), exist_ok=True)
        tmp = _ACCESS_FILE + ".tmp"
        with open(tmp, "w") as f:
            json.dump(session_virus_access, f, default=str)
        os.replace(tmp, _ACCESS_FILE)
    except Exception as e:
        print(f"Warning: could not persist session access records: {e}")

# Load any previously persisted jobs so downloads still work after a server restart
_load_jobs_from_disk()
_load_access_from_disk()

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


def _ensure_custom_virus_session_access(virus_name: str, session_id: str) -> None:
    '''Ensure a custom virus is accessible to the current session.

    Custom-virus file/result endpoints authorize via `analysis_jobs`.
    Upload/create flows can happen before analysis, so we register a lightweight
    sentinel job record here to grant session-level access immediately.
    '''
    sentinel_id = f"custom_{virus_name}_{session_id[:8]}"
    if sentinel_id in analysis_jobs:
        return

    analysis_jobs[sentinel_id] = {
        "job_id": sentinel_id,
        "virus_name": virus_name,
        "session_id": session_id,
        "status": "custom_virus_ready",
        "created_at": datetime.utcnow().isoformat() + "Z",
    }
    # Grant dedicated access ownership (separate from job tracking)
    owned = session_virus_access.get(session_id, [])
    if virus_name not in owned:
        owned.append(virus_name)
    session_virus_access[session_id] = owned
    _save_jobs_to_disk()
    _save_access_to_disk()


def _grant_session_access(virus_name: str, session_id: str) -> None:
    """Grant a session access to a virus resource."""
    owned = session_virus_access.get(session_id, [])
    if virus_name not in owned:
        owned.append(virus_name)
        session_virus_access[session_id] = owned
        _save_access_to_disk()


def _session_has_virus_access(virus_name: str, session_id: str) -> bool:
    """Check whether the current session owns/accesses this virus resources."""
    if virus_name in session_virus_access.get(session_id, []):
        return True
    # Backward compatibility: older sessions tracked only in analysis_jobs.
    for _, job in analysis_jobs.items():
        if job.get("session_id") == session_id and job.get("virus_name") == virus_name:
            return True
    return False

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
        _ensure_custom_virus_session_access(custom_virus.virus_name, session_id)
        
        return {
            "success": True,
            "virus_name": custom_virus.virus_name,
            "config": virus_config,
            "message": f"Custom virus '{custom_virus.virus_name}' created successfully"
        }
    except Exception as e:
        raise HTTPException(status_code=400, detail=str(e))

import re as _re

_VALID_NUCLEOTIDES = _re.compile(r'^[ACGTURYSWKMBDHVNacgturyswkmbdhvn\-\n\r ]+$')
_VALID_AMINO_ACIDS  = _re.compile(r'^[ACDEFGHIKLMNPQRSTVWYXacdefghiklmnpqrstvwyx\*\-\n\r ]+$')

def _validate_fasta_content(content: bytes, file_type: str) -> Optional[str]:
    """Return an error string if FASTA content is invalid, else None."""
    try:
        text = content.decode('utf-8', errors='replace')
    except Exception:
        return "File could not be decoded as UTF-8 text."

    lines = [l.strip() for l in text.splitlines()]
    if not lines or not lines[0].startswith('>'):
        if file_type == 'msa':
            # Allow CLUSTAL format for MSA
            if not (lines and lines[0].upper().startswith('CLUSTAL')):
                return "MSA file must be in FASTA (starts with '>') or CLUSTAL format."
            return None
        return "File does not appear to be valid FASTA format (must start with '>')."

    seq_chars = []
    for line in lines:
        if line.startswith('>') or not line:
            continue
        seq_chars.append(line)
    combined = ''.join(seq_chars)

    if not combined:
        return "File contains no sequence data."

    if file_type == 'reference_genome':
        invalid = _re.findall(r'[^ACGTURYSWKMBDHVNacgturyswkmbdhvn\-]', combined)
        if invalid:
            bad = ', '.join(sorted(set(invalid))[:10])
            return (f"Reference genome contains non-nucleotide characters: {bad}. "
                    "Only IUPAC nucleotide codes (A, C, G, T, U, R, Y, …) are allowed.")

    elif file_type == 'proteome':
        invalid = _re.findall(r'[^ACDEFGHIKLMNPQRSTVWYXacdefghiklmnpqrstvwyx\*\-]', combined)
        if invalid:
            bad = ', '.join(sorted(set(invalid))[:10])
            return (f"Proteome file contains non-amino-acid characters: {bad}. "
                    "Only standard one-letter amino acid codes are allowed.")

    return None


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

    # Read file content once for validation (avoids reading twice)
    file_content = await file.read()
    if not file_content:
        raise HTTPException(status_code=400, detail="Uploaded file is empty.")

    # Validate file content (nucleotide/protein characters)
    if file_type in ('reference_genome', 'proteome', 'msa'):
        validation_error = _validate_fasta_content(file_content, file_type)
        if validation_error:
            raise HTTPException(status_code=422, detail=validation_error)

    # Block uploading reference genome / proteome for built-in (non-custom) viruses.
    # Built-in viruses ship with their own curated reference files; overwriting them
    # silently would produce nonsensical results. Users should create a Custom Virus instead.
    if not is_custom_virus and file_type in ('reference_genome', 'proteome'):
        config = load_virus_config()
        virus_cfg = config.get("viruses", {}).get(virus_name, {})
        is_builtin = virus_cfg and not virus_cfg.get("description", "").startswith("Custom virus:")
        if is_builtin:
            raise HTTPException(
                status_code=400,
                detail=(
                    f"'{virus_name}' is a built-in virus with a curated reference genome and proteome. "
                    "Uploading your own reference or proteome for a built-in virus is not supported — "
                    "it would silently replace the shared reference and break other users' analyses. "
                    "To use your own reference files, please create a Custom Virus instead."
                )
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
        # Keep original filename for MSA files to preserve user's naming
        filename = original_filename
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
            buffer.write(file_content)  # already read for validation
        
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
                _ensure_custom_virus_session_access(virus_name, session_id)
        
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

@app.post("/api/load-sample-data")
async def load_sample_data(virus_name: str = "SARS-CoV-2", request: Request = None):
    """Grant the current session access to a built-in virus's sample data without an upload."""
    config = load_virus_config()
    if virus_name not in config.get("viruses", {}):
        raise HTTPException(status_code=404, detail=f"Virus '{virus_name}' not found")
    virus_cfg = config["viruses"][virus_name]
    if virus_cfg.get("description", "").startswith("Custom virus:"):
        raise HTTPException(status_code=400, detail="load-sample-data only works for built-in viruses")

    session_id = get_session_id(request) if request else "default"
    # Register a sentinel job so this session can access the virus data files
    sentinel_id = f"sample_{virus_name}_{session_id[:8]}"
    if sentinel_id not in analysis_jobs:
        analysis_jobs[sentinel_id] = {
            "job_id": sentinel_id,
            "virus_name": virus_name,
            "session_id": session_id,
            "status": "sample_loaded",
            "created_at": datetime.utcnow().isoformat() + "Z",
        }
        _save_jobs_to_disk()

    # Return available MSA files so the frontend can pre-select one
    msa_dir = f"data/{virus_name}/clustalW"
    msa_files = []
    if os.path.isdir(msa_dir):
        msa_files = [f for f in os.listdir(msa_dir)
                     if os.path.isfile(os.path.join(msa_dir, f))]
    return {"success": True, "virus_name": virus_name, "msa_files": sorted(msa_files)}


@app.get("/api/msa-files/{virus_name}")
async def get_msa_files(virus_name: str, segment: Optional[str] = None, request: Request = None):
    '''Get all available MSA files for a virus.'''
    try:
        config = load_virus_config()
        
        if virus_name not in config.get("viruses", {}):
            raise HTTPException(status_code=404, detail=f"Virus '{virus_name}' not found")
        
        virus_config = config["viruses"][virus_name]
        msa_files = []
        
        if segment and "segments" in virus_config:
            # Multi-segment virus
            if segment in virus_config["segments"]:
                segment_config = virus_config["segments"][segment]
                msa_dir = f"data/{virus_name}/{segment}/clustalW"
                if os.path.exists(msa_dir):
                    for file in os.listdir(msa_dir):
                        if file.endswith(('.txt', '.fasta', '.fa', '.aln', '.clustal')):
                            msa_files.append({
                                "filename": file,
                                "path": f"{msa_dir}/{file}",
                                "size": os.path.getsize(f"{msa_dir}/{file}"),
                                "is_default": file == segment_config.get("default_msa_file", "")
                            })
        else:
            # Single-segment virus
            msa_dir = f"data/{virus_name}/clustalW"
            if os.path.exists(msa_dir):
                for file in os.listdir(msa_dir):
                    if file.endswith(('.txt', '.fasta', '.fa', '.aln', '.clustal')):
                        msa_files.append({
                            "filename": file,
                            "path": f"{msa_dir}/{file}",
                            "size": os.path.getsize(f"{msa_dir}/{file}"),
                            "is_default": file == virus_config.get("default_msa_file", "")
                        })
        
        return {
            "success": True,
            "virus_name": virus_name,
            "segment": segment,
            "msa_files": msa_files
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error getting MSA files: {str(e)}")

@app.post("/api/set-default-msa")
async def set_default_msa_file(
    virus_name: str = Form(...),
    msa_filename: str = Form(...),
    segment: Optional[str] = Form(None),
    request: Request = None
):
    '''Set the default MSA file for a virus.'''
    try:
        config = load_virus_config()
        
        if virus_name not in config.get("viruses", {}):
            raise HTTPException(status_code=404, detail=f"Virus '{virus_name}' not found")
        
        virus_config = config["viruses"][virus_name]
        
        if segment and "segments" in virus_config:
            # Multi-segment virus
            if segment not in virus_config["segments"]:
                raise HTTPException(status_code=404, detail=f"Segment '{segment}' not found")
            virus_config["segments"][segment]["default_msa_file"] = msa_filename
        else:
            # Single-segment virus
            virus_config["default_msa_file"] = msa_filename
        
        save_virus_config(config)
        
        return {
            "success": True,
            "message": f"Default MSA file set to '{msa_filename}' for {virus_name}"
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error setting default MSA file: {str(e)}")

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
        "created_at": datetime.utcnow().isoformat() + "Z",
        "request": analysis_request.dict(),
        "results": None,
        "error": None,
        "session_id": session_id
    }
    
    analysis_jobs[job_id] = job_record
    _grant_session_access(analysis_request.virus_name, session_id)
    _save_jobs_to_disk()

    # Add background task
    background_tasks.add_task(run_analysis_job, job_id, analysis_request)
    
    return AnalysisResult(
        job_id=job_id,
        status="queued",
        virus_name=analysis_request.virus_name,
        created_at=job_record["created_at"]
    )

@app.get("/api/jobs")
async def get_user_jobs(request: Request):
    '''Get all analysis jobs for the current user session.'''
    session_id = get_session_id(request)
    
    # Filter jobs by session ID
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
    
    # Sort by creation time (newest first)
    user_jobs.sort(key=lambda x: x["created_at"], reverse=True)
    
    return {"jobs": user_jobs}

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
    
    # User isolation: pre-built viruses are always accessible; custom viruses
    # require the user to have uploaded or analyzed them in this session.
    config = load_virus_config()
    virus_cfg = config.get("viruses", {}).get(virus_name, {})
    is_builtin = virus_cfg and not virus_cfg.get("description", "").startswith("Custom virus:")

    if not is_builtin:
        if not _session_has_virus_access(virus_name, session_id):
            raise HTTPException(status_code=403, detail="Access denied: You can only access data for viruses you have uploaded or analyzed")
    
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
                elif "fasta" in file_path:
                    # Files in fasta/ folder are FASTA sequences for MSA generation
                    file_type = "fasta_sequences"
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
                "modified": datetime.utcfromtimestamp(os.path.getmtime(file_path)).isoformat() + "Z"
            })
    
    return {"files": files}

@app.get("/api/result-files/{virus_name}")
async def get_result_files(virus_name: str, request: Request, segment: Optional[str] = None):
    '''Get list of result files for a virus.'''
    session_id = get_session_id(request)

    # Pre-built viruses always have accessible results; only custom viruses
    # are restricted to the session that ran the analysis.
    config = load_virus_config()
    virus_cfg = config.get("viruses", {}).get(virus_name, {})
    is_builtin = virus_cfg and not virus_cfg.get("description", "").startswith("Custom virus:")

    if not is_builtin:
        if not _session_has_virus_access(virus_name, session_id):
            raise HTTPException(status_code=403, detail="Access denied: You can only access results for viruses you have analyzed")
    
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
                    "modified": datetime.utcfromtimestamp(os.path.getmtime(file_path)).isoformat() + "Z"
                })
    
    return {"files": files}

@app.get("/api/download-file")
async def download_file(file_path: str, request: Request):
    '''Download a specific file from the data or result directories (session-based access).'''
    session_id = get_session_id(request)
    
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="File not found")
    
    # Security check: ensure file is within allowed directories
    allowed_dirs = ["data/", "result/", "uploads/", "results/"]
    if not any(file_path.startswith(allowed_dir) for allowed_dir in allowed_dirs):
        raise HTTPException(status_code=403, detail="Access denied")
    
    # User isolation: Check if user has access to this file
    config = load_virus_config()

    if file_path.startswith("result/"):
        # For result files, check if they belong to user's analysis jobs.
        # Pre-built virus results are always accessible.
        virus_name = file_path.split("/")[1] if len(file_path.split("/")) > 1 else None
        if virus_name:
            virus_cfg = config.get("viruses", {}).get(virus_name, {})
            is_builtin = virus_cfg and not virus_cfg.get("description", "").startswith("Custom virus:")
            if not is_builtin:
                if not _session_has_virus_access(virus_name, session_id):
                    raise HTTPException(status_code=403, detail="Access denied: You can only download your own analysis results")

    elif file_path.startswith("data/"):
        # For data files, check if user has uploaded files for this virus.
        # Pre-built virus data is always accessible.
        virus_name = file_path.split("/")[1] if len(file_path.split("/")) > 1 else None
        if virus_name:
            virus_cfg = config.get("viruses", {}).get(virus_name, {})
            is_builtin = virus_cfg and not virus_cfg.get("description", "").startswith("Custom virus:")
            if not is_builtin:
                if not _session_has_virus_access(virus_name, session_id):
                    raise HTTPException(status_code=403, detail="Access denied: You can only download data for viruses you have uploaded")
    
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
    session_id = get_session_id(request)
    
    # Security check: ensure virus name is valid
    if not virus_name or ".." in virus_name or "/" in virus_name:
        raise HTTPException(status_code=400, detail="Invalid virus name")
    
    # User isolation: pre-built viruses are always downloadable;
    # custom viruses require a session job.
    config = load_virus_config()
    virus_cfg = config.get("viruses", {}).get(virus_name, {})
    is_builtin = virus_cfg and not virus_cfg.get("description", "").startswith("Custom virus:")

    if not is_builtin:
        if not _session_has_virus_access(virus_name, session_id):
            raise HTTPException(status_code=403, detail="Access denied: You can only download files for viruses you have uploaded or analyzed")

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

def _user_friendly_error(e: Exception, context: str = "") -> str:
    """Convert an internal exception into a readable message for the job record."""
    msg = str(e)
    prefix = f"[{context}] " if context else ""

    # Common, actionable patterns
    if "No MSA file" in msg or "default_msa_file" in msg.lower():
        return (prefix + "No MSA (Multiple Sequence Alignment) file is configured for this virus. "
                "Please upload an MSA file in the Data Upload tab before running analysis.")
    if "not found" in msg.lower() and ("msa" in msg.lower() or "clustalw" in msg.lower()):
        return (prefix + f"MSA file could not be found on the server: {msg}. "
                "Try re-uploading your MSA file.")
    if "reference genome" in msg.lower() or "ref_genome" in msg.lower():
        return (prefix + "Reference genome file is missing or misconfigured. "
                "Ensure you have uploaded a valid reference genome for this virus.")
    if "proteome" in msg.lower():
        return (prefix + "Proteome file is missing or misconfigured. "
                "Ensure you have uploaded a valid proteome file for this virus.")
    if "FileNotFoundError" in type(e).__name__ or "No such file" in msg:
        return (prefix + f"A required file was not found: {msg}. "
                "Check that all required data files have been uploaded.")
    if "KeyError" in type(e).__name__:
        return (prefix + f"Configuration key missing: {msg}. "
                "This may indicate an incomplete virus configuration.")
    if "ValueError" in type(e).__name__:
        return (prefix + f"Invalid value encountered during analysis: {msg}. "
                "Check that your input files are correctly formatted.")
    # Generic fallback — include the raw message but add guidance
    return prefix + msg + (" — Check server logs for details." if len(msg) < 80 else "")


async def run_analysis_job(job_id: str, analysis_request: AnalysisRequest):
    '''Run the actual analysis in the background.'''
    try:
        analysis_jobs[job_id]["status"] = "running"

        # ── Resolve MSA file path ────────────────────────────────────────────
        msa_file_path = None

        if analysis_request.msa_file:
            virus_data_dir = f"data/{analysis_request.virus_name}"
            if analysis_request.segment:
                virus_data_dir += f"/{analysis_request.segment}"
            msa_file_path = f"{virus_data_dir}/clustalW/{analysis_request.msa_file}"
            if not os.path.exists(msa_file_path):
                raise FileNotFoundError(
                    f"Specified MSA file '{analysis_request.msa_file}' not found at {msa_file_path}."
                )
        else:
            config = load_virus_config()
            virus_config = config.get("viruses", {}).get(analysis_request.virus_name, {})

            if not virus_config:
                raise KeyError(
                    f"Virus '{analysis_request.virus_name}' not found in configuration. "
                    "Create a custom virus entry before running analysis."
                )

            if analysis_request.segment and "segments" in virus_config:
                default_msa = virus_config["segments"].get(analysis_request.segment, {}).get("default_msa_file", "")
            else:
                default_msa = virus_config.get("default_msa_file", "")

            if default_msa:
                virus_data_dir = f"data/{analysis_request.virus_name}"
                if analysis_request.segment:
                    virus_data_dir += f"/{analysis_request.segment}"
                msa_file_path = f"{virus_data_dir}/clustalW/{default_msa}"
                if not os.path.exists(msa_file_path):
                    raise FileNotFoundError(
                        f"Configured MSA file '{default_msa}' not found. "
                        "Please re-upload the MSA file."
                    )
            else:
                raise ValueError(
                    f"No MSA file is configured for '{analysis_request.virus_name}'. "
                    "Upload an MSA file in the Data Upload tab first."
                )
        
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
        
        # Generate visualization(s): requested type plus PRF when frameshift detection is enabled.
        visualization_files = []
        visualization_types_to_generate = []
        if analysis_request.visualization_type:
            visualization_types_to_generate.append(analysis_request.visualization_type)
        if analysis_request.detect_frameshifts and "prf" not in visualization_types_to_generate:
            visualization_types_to_generate.append("prf")

        if visualization_types_to_generate:
            import subprocess
            import sys

            output_dir = f"imgs/visualizations/{analysis_request.virus_name}"
            os.makedirs(output_dir, exist_ok=True)

            type_patterns = {
                "mutation": ["combined_analysis", "mutation"],
                "row-hot": ["row_hot_mutations", "row-hot"],
                "prf": ["prf_regions", "prf"]
            }

            for viz_type in visualization_types_to_generate:
                try:
                    cmd = [sys.executable, "plot.py", "--type", viz_type, "--virus", analysis_request.virus_name]
                    if analysis_request.genome_id:
                        cmd.extend(["--genome-id", analysis_request.genome_id])

                    result = subprocess.run(cmd, capture_output=True, text=True, cwd=".")
                    if result.returncode != 0:
                        print(f"Visualization '{viz_type}' failed: {result.stderr}")
                        continue

                    current_time = time.time()
                    for file in os.listdir(output_dir):
                        if not file.endswith(".html"):
                            continue
                        file_lower = file.lower()
                        if not any(pattern in file_lower for pattern in type_patterns.get(viz_type, [])):
                            continue

                        file_path = os.path.join(output_dir, file)
                        if current_time - os.path.getmtime(file_path) < 300 and file_path not in visualization_files:
                            visualization_files.append(file_path)

                    print(f"Visualization '{viz_type}' generated successfully")
                except Exception as viz_error:
                    print(f"Error generating visualization '{viz_type}': {viz_error}")
        
        # Update job with results
        analysis_jobs[job_id]["status"] = "completed"
        analysis_jobs[job_id]["results"] = {
            "files": result_files,
            "visualization_files": visualization_files,
            "summary": f"Analysis completed for {analysis_request.virus_name}"
        }
        _save_jobs_to_disk()

    except Exception as e:
        analysis_jobs[job_id]["status"] = "failed"
        analysis_jobs[job_id]["error"] = _user_friendly_error(e, "Analysis")
        print(f"Analysis job {job_id} failed: {e}")
        _save_jobs_to_disk()

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

@app.get("/api/genome-ids/{virus_name}")
async def get_genome_ids(virus_name: str, segment: Optional[str] = None, msa_file: Optional[str] = None):
    '''Get genome IDs from MSA file for a specific virus.'''
    try:
        # Load virus configuration
        config = load_virus_config()
        virus_config = config.get("viruses", {}).get(virus_name, {})

        if not virus_config:
            raise HTTPException(status_code=404, detail=f"Virus '{virus_name}' not found in configuration")

        # Determine MSA file path
        if segment and "segments" in virus_config:
            default_msa = virus_config["segments"].get(segment, {}).get("default_msa_file", "")
            virus_data_dir = f"data/{virus_name}/{segment}"
        else:
            default_msa = virus_config.get("default_msa_file", "")
            virus_data_dir = f"data/{virus_name}"

        if not default_msa and not msa_file:
            raise HTTPException(status_code=404, detail=f"No default MSA file configured for virus '{virus_name}'")

        # Allow caller to override the MSA file (e.g. user selected a non-default MSA)
        chosen_msa = msa_file if msa_file else default_msa
        msa_file_path = f"{virus_data_dir}/clustalW/{chosen_msa}"
        
        if not os.path.exists(msa_file_path):
            raise HTTPException(status_code=404, detail=f"MSA file not found: {msa_file_path}")
        
        # Get reference genome ID from virus configuration
        reference_genome_id = None
        if segment and "segments" in virus_config:
            reference_genome_id = virus_config["segments"].get(segment, {}).get("reference_genome", "")
        else:
            reference_genome_id = virus_config.get("reference_genome", "")
        
        # Remove file extension from reference genome ID for comparison
        if reference_genome_id:
            reference_genome_id = reference_genome_id.replace('.fasta', '').replace('.fa', '')
            print(f"DEBUG: Reference genome ID for {virus_name}: {reference_genome_id}")
        
        # Extract genome IDs from MSA file
        genome_ids = []
        
        try:
            with open(msa_file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    # Skip empty lines and CLUSTAL header lines
                    if not line or line.startswith('CLUSTAL') or line.startswith('*') or line.startswith(' '):
                        continue
                    
                    # Extract genome ID from lines that contain sequence data
                    # Format: "GENOME_ID    SEQUENCE_DATA"
                    parts = line.split()
                    if len(parts) >= 2:
                        full_genome_id = parts[0]
                        
                        # Handle different genome ID formats
                        if '|' in full_genome_id:
                            # Format: "NAME|EPI_ISL_XXXXX|DATE" or "NAME|DESCRIPTION"
                            pipe_parts = full_genome_id.split('|')
                            if len(pipe_parts) >= 2:
                                # Check if second part looks like an EPI_ISL ID
                                if pipe_parts[1].startswith('EPI_ISL_'):
                                    genome_id = pipe_parts[1]  # Use EPI_ISL ID
                                elif len(pipe_parts) >= 3 and pipe_parts[2].startswith('NC_'):
                                    genome_id = pipe_parts[2]  # Use NC_ ID (reference genome)
                                else:
                                    genome_id = pipe_parts[0]  # Use first part
                            else:
                                genome_id = pipe_parts[0]
                        else:
                            genome_id = full_genome_id
                        
                        # Skip the reference genome (using configuration, not hardcoded patterns)
                        if reference_genome_id and genome_id == reference_genome_id:
                            print(f"DEBUG: Skipping reference genome: {genome_id}")
                            continue
                        
                        # Add non-reference genome IDs
                        if genome_id not in genome_ids:
                            genome_ids.append(genome_id)
        except Exception as e:
            raise HTTPException(status_code=500, detail=f"Error reading MSA file: {str(e)}")
        
        return {"genome_ids": genome_ids}
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error extracting genome IDs: {str(e)}")

# Visualization endpoints
@app.post("/api/visualize")
async def create_visualization(
    virus_name: str = Form(...),
    visualization_type: str = Form(...),  # 'mutation', 'row-hot', 'prf'
    genome_id: Optional[str] = Form(None),
    output_path: Optional[str] = Form(None),
    request: Request = None
):
    '''Create visualization plots for a virus.'''
    try:
        # Validate visualization type
        valid_types = ['mutation', 'row-hot', 'prf']
        if visualization_type not in valid_types:
            raise HTTPException(status_code=400, detail=f"Invalid visualization type. Must be one of: {valid_types}")
        
        # Check if virus exists
        virus_data_dir = f"data/{virus_name}"
        if not os.path.exists(virus_data_dir):
            raise HTTPException(status_code=404, detail=f"Virus '{virus_name}' not found")
        
        # Create output directory
        output_dir = f"imgs/visualizations/{virus_name}"
        os.makedirs(output_dir, exist_ok=True)
        
        # Run the visualization using the unified plot script
        import subprocess
        import sys
        
        cmd = [
            sys.executable, "plot.py",
            "--type", visualization_type,
            "--virus", virus_name
        ]
        
        if genome_id:
            cmd.extend(["--genome-id", genome_id])
        if output_path:
            cmd.extend(["--output", output_path])
        
        # Run the command
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=os.getcwd()
        )
        
        # Check if visualization failed but still try to find existing files
        print(f"DEBUG: returncode={result.returncode}, stderr contains 'Error:': {'Error:' in result.stderr}")
        print(f"DEBUG: stderr content: {repr(result.stderr)}")
        
        # Check for specific errors that indicate complete failure
        kaleido_error = "Kaleido requires Google Chrome" in result.stderr
        general_error = "Error:" in result.stderr
        command_failed = result.returncode != 0
        
        visualization_failed = command_failed or general_error
        if visualization_failed:
            print(f"Warning: Visualization command failed: {result.stderr}")
            # Continue to try to find existing files
        else:
            # Check if there was an error in stderr even if return code was 0
            if general_error:
                visualization_failed = True
        
        # Find the generated files
        generated_files = []
        html_file = None
        pdf_file = None
        
        if os.path.exists(output_dir):
            # Get the most recent files that match the visualization type
            files = []
            for file in os.listdir(output_dir):
                if file.endswith(('.html', '.pdf')):
                    file_path = os.path.join(output_dir, file)
                    files.append((file, os.path.getmtime(file_path)))
            
            # Sort by modification time (most recent first)
            files.sort(key=lambda x: x[1], reverse=True)
            
            # Map visualization types to file name patterns
            type_patterns = {
                'mutation': ['combined_analysis', 'mutation'],
                'row-hot': ['row_hot_mutations', 'row-hot'],
                'prf': ['prf_regions', 'prf']
            }
            
            # Find the most recent files that match the visualization type
            print(f"DEBUG: Checking files in {output_dir}: {[f[0] for f in files]}")
            
            # Only consider files created in the last 5 minutes (300 seconds)
            current_time = time.time()
            recent_files = [(file, mtime) for file, mtime in files if current_time - mtime < 300]
            print(f"DEBUG: Recent files (last 5 minutes): {[f[0] for f in recent_files]}")
            
            for file, mtime in recent_files:
                file_matches = False
                if visualization_type in type_patterns:
                    for pattern in type_patterns[visualization_type]:
                        if pattern in file.lower():
                            file_matches = True
                            break
                
                if file_matches:
                    if file.endswith('.html') and html_file is None:
                        html_file = {
                            "filename": file,
                            "path": f"/api/visualization-files/{virus_name}/{file}",
                            "type": "html"
                        }
                        print(f"DEBUG: Found HTML file: {file}")
                    elif file.endswith('.pdf') and pdf_file is None:
                        pdf_file = {
                            "filename": file,
                            "path": f"/api/visualization-files/{virus_name}/{file}",
                            "type": "pdf"
                        }
        
        # Don't include PDF files - only HTML for embedding
        # if pdf_file:
        #     generated_files.append(pdf_file)
        
        # Prepare response message and success status
        success = True
        if html_file:
            # Files were found (either newly created or existing)
            if visualization_failed:
                # HTML was generated despite errors (likely just PDF failed)
                message = f"Visualization '{visualization_type}' HTML generated successfully, but PDF generation failed (Chrome required)"
            else:
                message = f"Visualization '{visualization_type}' created successfully"
        else:
            # No files were found at all
            if kaleido_error:
                message = f"Visualization '{visualization_type}' failed to generate files (Chrome required for PDF)"
            elif general_error:
                message = f"Visualization '{visualization_type}' failed to generate files"
            else:
                message = f"Visualization '{visualization_type}' completed but no files were generated"
            success = False
        
        # Update visualization_failed based on stderr content
        if kaleido_error:
            if html_file:
                message = f"Visualization '{visualization_type}' HTML generated successfully, but PDF generation failed (Chrome required)"
            else:
                message = f"Visualization '{visualization_type}' failed to generate files (Chrome required for PDF)"
                success = False
        elif general_error:
            if html_file:
                message = f"Visualization '{visualization_type}' HTML generated successfully, but some errors occurred"
            else:
                message = f"Visualization '{visualization_type}' failed to generate files"
                success = False
        
        return {
            "success": success,
            "message": message,
            "virus_name": virus_name,
            "visualization_type": visualization_type,
            "files": [],  # No files to download
            "html_file": html_file,
            "output": result.stdout,
            "visualization_failed": visualization_failed
        }
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error creating visualization: {str(e)}")

@app.get("/api/visualization-files/{virus_name}/{filename}")
async def get_visualization_file(virus_name: str, filename: str):
    '''Serve visualization files.'''
    file_path = f"imgs/visualizations/{virus_name}/{filename}"
    
    if not os.path.exists(file_path):
        raise HTTPException(status_code=404, detail="Visualization file not found")
    
    # Set appropriate media type based on file extension
    if filename.endswith('.html'):
        media_type = "text/html"
    elif filename.endswith('.pdf'):
        media_type = "application/pdf"
    else:
        media_type = "application/octet-stream"
    
    return FileResponse(
        file_path,
        media_type=media_type
    )

@app.get("/api/viz-html/{virus_name}")
async def get_viz_html(virus_name: str):
    """Return the most-recent HTML visualization file per type for a virus.

    Returns a dict like:
      { "mutation": "MT526800_combined_analysis.html",
        "row_hot":  "EPI_ISL_20090675_row_hot_mutations.html",
        "prf":      "potential_PRF_prf_regions.html" }
    Missing types have a null value.
    """
    viz_dir = f"imgs/visualizations/{virus_name}"
    if not os.path.exists(viz_dir):
        return {"mutation": None, "row_hot": None, "prf": None}

    # Collect HTML files grouped by type (newest first)
    mutation_files, row_hot_files, prf_files = [], [], []
    for fname in os.listdir(viz_dir):
        if not fname.endswith(".html"):
            continue
        fpath = os.path.join(viz_dir, fname)
        mtime = os.path.getmtime(fpath)
        if "combined_analysis" in fname:
            mutation_files.append((mtime, fname))
        elif "row_hot" in fname:
            row_hot_files.append((mtime, fname))
        elif "prf" in fname:
            prf_files.append((mtime, fname))

    pick = lambda lst: sorted(lst, reverse=True)[0][1] if lst else None
    return {
        "mutation": pick(mutation_files),
        "row_hot":  pick(row_hot_files),
        "prf":      pick(prf_files),
    }


@app.get("/api/visualization-types")
async def get_visualization_types():
    '''Get available visualization types.'''
    return {
        "visualization_types": [
            {
                "id": "mutation",
                "name": "Mutation Analysis",
                "description": "Combined genome organization and protein mutation analysis"
            },
            {
                "id": "row-hot", 
                "name": "Row/Hot Mutations",
                "description": "Visualization of row and hot mutations on protein bars"
            },
            {
                "id": "prf",
                "name": "PRF Regions", 
                "description": "Programmed Ribosomal Frameshifting regions visualization"
            }
        ]
    }

# ─────────────────────────────────────────────────────────────────────────────
# RNA Dual Graph routes
# ─────────────────────────────────────────────────────────────────────────────

def _get_dual_session_dir(session_id: str) -> str:
    """Return a session-specific scratch directory for the dual-graph pipeline."""
    path = os.path.join("uploads", "dual_graph", session_id)
    for sub in ["PDB_DBN", "PDB_DSSR", "PDB_DSSR_2D", "PDB_DSSR_CT", "PDB_DSSR_Dual"]:
        os.makedirs(os.path.join(path, sub), exist_ok=True)
    return path


async def run_dual_graph_job(job_id: str, dbn_content: str, structure_id: str, session_id: str):
    """Background task: run the 3-step dual graph pipeline for one DBN input."""
    analysis_jobs[job_id]["status"] = "running"
    base_dir = _get_dual_session_dir(session_id)

    try:
        if not _DUAL_SEARCH_AVAILABLE:
            raise RuntimeError("Existing-Dual-Search modules could not be imported. Check server logs.")

        # Step 1 — write .dbn and convert to DSSR format
        dbn_path = os.path.join(base_dir, "PDB_DBN", f"{structure_id}-2D.dbn")
        dssr_path = os.path.join(base_dir, "PDB_DSSR", f"{structure_id}.out")
        with open(dbn_path, "w") as f:
            f.write(dbn_content)
        convert_dbn_to_dssr(dbn_path, dssr_path, structure_id)

        # Step 2 — DSSR → 2D summary + CT files
        pdbto2d_process_id(structure_id, base_dir)

        # Step 3 — CT → dual graph IDs
        results = run_dual_library([structure_id], base_dir)

        analysis_jobs[job_id]["status"] = "completed"
        analysis_jobs[job_id]["results"] = {
            "graph_assignments": results["graph_assignments"],
            "no_ct": results["no_ct"],
            "base_dir": base_dir,
        }
        _save_jobs_to_disk()

    except Exception as e:
        analysis_jobs[job_id]["status"] = "failed"
        analysis_jobs[job_id]["error"] = str(e)
        _save_jobs_to_disk()


@app.post("/api/dual-graph")
async def submit_dual_graph(
    background_tasks: BackgroundTasks,
    request: Request,
    dbn_text: Optional[str] = Form(None),
    dbn_file: Optional[UploadFile] = File(None),
    structure_id: str = Form("USER_INPUT"),
):
    """Submit an RNA secondary structure (dot-bracket) for dual graph topology assignment."""
    if not dbn_text and not dbn_file:
        raise HTTPException(status_code=400, detail="Provide either dbn_text or dbn_file.")

    session_id = get_session_id(request)
    job_id = str(uuid.uuid4())

    # Read file upload if provided
    if dbn_file:
        dbn_content = (await dbn_file.read()).decode("utf-8")
        if structure_id == "USER_INPUT":
            structure_id = Path(dbn_file.filename).stem.split("-")[0].upper() or "INPUT"
    else:
        dbn_content = dbn_text

    analysis_jobs[job_id] = {
        "job_id": job_id,
        "type": "dual_graph",
        "status": "queued",
        "virus_name": f"dual_graph:{structure_id}",
        "structure_id": structure_id,
        "created_at": datetime.utcnow().isoformat() + "Z",
        "session_id": session_id,
        "results": None,
        "error": None,
    }

    background_tasks.add_task(run_dual_graph_job, job_id, dbn_content, structure_id, session_id)

    response = JSONResponse({"job_id": job_id, "status": "queued", "structure_id": structure_id})
    if not request.cookies.get("session_id"):
        response.set_cookie(key="session_id", value=session_id, httponly=True, max_age=3600 * 24 * 7)
    return response


@app.get("/api/dual-graph/plot/{graph_id}")
async def get_dual_graph_plot(graph_id: str):
    """Return a PNG visualization for a dual graph topology ID (e.g. '3_4').

    PNGs are cached in uploads/dual_graph/plots/ so repeated requests are instant.
    """
    # Sanitize: only allow digits and underscores
    import re
    if not re.fullmatch(r"\d+_\d+", graph_id):
        raise HTTPException(status_code=400, detail="Invalid graph_id format. Expected e.g. '3_4'.")

    cache_dir = os.path.join("uploads", "dual_graph", "plots")
    png_path  = os.path.join(cache_dir, f"{graph_id}.png")

    if not os.path.exists(png_path):
        if not _DUAL_SEARCH_AVAILABLE:
            raise HTTPException(status_code=503, detail="Dual graph modules not available.")
        ok = plot_dual_graph(graph_id, png_path)
        if not ok:
            raise HTTPException(status_code=404, detail=f"Graph '{graph_id}' not found in library.")

    return FileResponse(png_path, media_type="image/png")


@app.get("/api/dual-graph/result-files/{session_id}/{structure_id}")
async def get_dual_graph_files(session_id: str, structure_id: str, request: Request):
    """List result files produced by the dual graph pipeline for a given session+structure."""
    req_session = get_session_id(request)
    if req_session != session_id:
        raise HTTPException(status_code=403, detail="Access denied")

    base_dir = os.path.join("uploads", "dual_graph", session_id)
    dual_dir = os.path.join(base_dir, "PDB_DSSR_Dual")

    if not os.path.exists(dual_dir):
        return {"files": []}

    files = []
    for fname in os.listdir(dual_dir):
        if fname.endswith(".txt") and fname.startswith(structure_id):
            fpath = os.path.join(dual_dir, fname)
            with open(fpath) as f:
                content = f.read().strip()
            files.append({"filename": fname, "graph_id": content})
    return {"files": files}


# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    import os
    port = int(os.environ.get("PORT", 8000))
    uvicorn.run(app, host="0.0.0.0", port=port, reload=False)
