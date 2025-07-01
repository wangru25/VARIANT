#!/usr/bin/env python3
"""
Setup script for MutParser package.

This file provides backward compatibility for pip installation.
For modern Python packaging, use pyproject.toml instead.
"""

from setuptools import find_packages, setup

if __name__ == "__main__":
    setup(
        name="viralytics-mut",
        version="1.0.0",
        description="A Comprehensive Framework for Multi-Scale Viral Mutation Analysis with Integrated Programmed Ribosomal Frameshifting Detection",
        author="Rui Wang",
        author_email="rw3594@nyu.edu",
        packages=find_packages(where="src"),
        package_dir={"": "src"},
        python_requires=">=3.8",
        install_requires=[
            "biopython>=1.79",
            "numpy>=1.21.0",
            "pandas>=1.3.0",
            "scikit-learn>=1.0.0",
            "fuzzysearch==0.7.3",
            "more-itertools>=8.12.0",
            "pyyaml>=6.0",
        ],
        extras_require={
            "dev": [
                "pytest>=7.0.0",
                "pytest-cov>=4.0.0",
                "black>=22.0.0",
                "isort>=5.10.0",
                "flake8>=5.0.0",
                "mypy>=1.0.0",
                "pre-commit>=2.20.0",
            ]
        },
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
    )
