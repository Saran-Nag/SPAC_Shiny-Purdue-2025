# Quickstart Guide: SPAC Visualization on Windows with WSL
This guide provides detailed instructions for setting up and launching the Spacy Visualization dashboard on a Windows machine using WSL (Windows Subsystem for Linux). It supports local testing and development workflows including environment setup, dependency resolution, and launching the app for analysis.

## Prerequisites
Before starting, ensure you’ve installed the following:
- Git for Windows
- Miniconda (Python=3.9.13)
- Ubuntu on WSL
- Visual Studio Code (WSL extension recommended)
Optional: GitHub CLI or access tokens for authentication.

## Windows + WSL Setup Instructions
### 1. Setup your Workspace 

```bash
mkdir -p ~/summer2025
cd ~/summer2025
```

### 2. Cloning the SPAC Shiny repository

After installing Git for Windows, clone the modularize branch to ensure the latest version of the app is available. Then switch the updated branch 

```bash
git clone https://github.com/Summer2025-SPAC/SPAC_Shiny 
cd SPAC_Shiny 
git checkout keyerror_fix 
```

### 3. Creating the Conda environment

Using miniconda, stop any current environment from running and create a virutal environment called spac_env_3119 using python version 3.11.9. Use the next command to activate the environment.

```bash
conda deactivate
conda create -n spacy_viz_py39_wsl python=3.9.13 -y
conda activate spacy_viz_py39_wsl
```

### 4. Installling dependencies

Install the dependencies from the requirements.txt file. Ensure python=3.9.13

```bash
pip install -r requirements.txt
```

### 5. Open VS Code in WSL

```bash
code ./
```

### 6. Run the app.
Use this command to run the SPAC app in the local browser.

```bash
shiny run --reload app.py
```

Visit http://127.0.0.1:8000 in your browser. Use this sample [dataset](https://zenodo.org/records/6376767/files/healthy_lung_adata.h5ad?download=1) set to test out the app. 

