# SPAC_Shiny – Windows Installation for Modularize Branch


## Branch

This setup guide is for the `modularize` branch.

---

## Windows Installation using Ubuntu

### Prerequisites

- [Git for Windows](https://git-scm.com/download/win)
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- [Visual Studio Code](https://code.visualstudio.com/) (optional)
- [wsl --install](https://learn.microsoft.com/en-us/windows/wsl/install)

### Steps

1. **Access Ubuntu**
    In terminal
    ```bash
    wsl --install

2. **Clone the repository**
   ```bash
   git clone --branch modularize https://github.com/<your-username>/SPAC_Shiny.git
   cd SPAC_Shiny

3. **Create and activate environment**
    ```bash
    conda create -n spac_env_3119 python=3.11
    conda activate spac_env_3119

4. **Install dependencies**
    ```bash
    pip install --upgrade pip wheel
    pip install scikit-learn==1.3.2 --only-binary :all:
    pip install llvmlite==0.41.1 --only-binary :all:
    pip install shiny scanpy
    pip install --no-deps -r requirements.txt


5. **Launch the app**
    ```bash
    python app.py

    Visit http://127.0.0.1:8000 in order to fully access the local web server
    
