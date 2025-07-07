# SPAC Shiny: Windows Installation Guide

The SPAC Shiny App is a web-based dashboard that allows researchers to analyze single-cell data. This dahboard contains a user friendly interface with a customizable graphs, along with a modular pipeline for examining images. Use this installation guide for WSL environments to setup the SPAC Shiny dashboard.

This setup guide is for the `modularize` branch.


## Windows Installation using Ubuntu

### Prerequisites
Git for Windows, Miniconda, Vscode, and Ubuntu are all used to setup and launch this app.

- [Git for Windows](https://git-scm.com/download/win)
- [Miniconda (Python 3.11+)](https://docs.conda.io/en/latest/miniconda.html)  
- [Ubuntu](https://ubuntu.com/download)
- [Visual Studio Code](https://code.visualstudio.com/) (optional)

### Steps

1. **Access Ubuntu**

    In the terminal, open Ubuntu usin this command.
    ```bash
    wsl --install
2. **Clone the repository**
   
   After installing Git for Windows, clone the modularize branch to ensure the latest version of the app is available.
   ```bash
   git clone --branch modularize https://github.com/FNLCR-DMAP/SPAC_Shiny.git
   cd SPAC_Shiny
3. **Create and activate environment**

    Using miniconda, create a virutal environment called spac_env_3119 using python version 3.11.9. Use the next command to activate the environment.
    ```bash
    conda create -n spac_env_3119 python=3.11.9
    conda activate spac_env_3119
4. **Install dependencies**

    Use the upgrade command to ensure you are using the latest version of pip. Then install scikit learn and llvmlite with the command. This forces pip to use precompiled binary wheels instead of building from source — which avoids messy compiler errors, especially on Windows or WSL setups.
    
    ```bash
    pip install --upgrade pip wheel
    pip install scikit-learn==1.3.2 --only-binary :all:
    pip install llvmlite==0.41.1 --only-binary :all:
    pip install shiny scanpy
    pip install --no-deps -r requirements.txt
5. **Launch the app**

    Launch the app using python app.py. 
    ```bash
    python app.py
    ```

    Visit http://127.0.0.1:8000 in order to fully access the local web server. If the local server is not working,  try  
    ```bash
    shiny run --reload app.py
    ```
    This starts the Shiny server using the CLI tool. 
