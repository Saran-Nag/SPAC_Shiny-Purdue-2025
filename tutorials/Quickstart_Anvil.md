This is a guide for Frederick-Analysis 2025-2026 Team to run shiny on Anvil, which may also be helpful for other HPC users. 

This guide uses a shared conda environment. You may also create your own one. Check other tutorial files or conda's official site for reference.

Warning: The Data Mine Team does not recommend using conda environments due to debugging issues. They may not help debug issues that arise from personal conda environment. For future Data Mine Teams, please think carefully before using it.


Author: Boqiang Zhang

Date: 25.10.16


# Step 1: Connect to Anvil through VS Code

In order to run shiny, we have to access Anvil through an interactive app, and the only choice I find that works is VS Code. (If there are another ways, please let me know!) For Data Mine students, here are the steps to connect to Anvil through VS Code:

1. Go to https://ondemand.anvil.rcac.purdue.edu/, log in with your ACCESS account.

2. Go to "The Data Mine" - "Visual Studio Code".

3. In the new page, select cores and time if needed. Click "Launch". 

4. Wait until your session is ready, then click "Connect to VS Code". It may take a few minutes.

You will start a VS Code session on Anvil in a new tab.



# Step 2: Set up environment for shiny

## 1. Open a terminal and go to the shiny repo

If you haven't cloned the shiny repo to your own folder, please do it first. The link to our fork is: https://github.com/Saran-Nag/SPAC_Shiny-Purdue-2025. Or use the SPAC_Shiny repo with this command:

```bash
git clone https://github.com/FNLCR-DMAP/SPAC_Shiny.git
```

Please also make sure you are under the shiny repo's directory before any further steps. Run the following command if necessary:

```bash
cd <home directory of SPAC_Shiny folder>
```


## 2. Load the shared conda environment for shiny

To load the shared conda environment, run the following commands each time you open a new terminal:

```bash
module load conda
module use $PROJECT/corporate/frederick-analysis/2025-2026/etc/modules
module load conda-env/shiny_shared-py3.12.8
```

Try running "module purge" first to unload unnecessary modules if needed.

For those who want to know how to create this shared environment module, please refer to **Appendix A** for more information

You may also activate the conda environment using the standard "conda activate" command. See **Appendix B** for more information.



# Step 3: Run shiny

Now we are able to run shiny. Run the following commands in the same terminal:

```bash
shiny run --reload app.py
```

Wait for a while until a pop-up window appears. Click "Open in Browser". You may also find the forwarded address in "PORTS". The default port # is 8000. 

You have sucessfully run the shiny app on Anvil!

Use this dataset to test: https://zenodo.org/records/6376767/files/healthy_lung_adata.h5ad?download=1 



# Cheat Sheet to run shiny on Anvil

- Connect to Anvil through VS Code
- Set up environment
```bash
module load conda
module use $PROJECT/corporate/frederick-analysis/2025-2026/etc/modules
module load conda-env/shiny_shared-py3.12.8
```
- Run Shiny
```bash
cd /path/to/your/repo/directory
shiny run --reload app.py
```



# Appendix A. How to create the shared conda environment from environment.yml

Reference: https://www.rcac.purdue.edu/knowledge/anvil/run/examples/apps/python/packages. This link also teaches you how to generate a Jupyter kernel from a conda environment. 

If you are using VS Code on Anvil, please make sure to select at least 3 cores (that's 5~6 GB RAM) when you start your session. Otherwise the creating process may crash.

The following are the commands I use to create the shared environment. Replace the location by your designated one. 

```bash
# Load conda module
module load conda

# Create an conda environment from environment.yml. Use --prefix plus the location you want to store your environment
conda env create -f environment.yml --prefix $PROJECT/corporate/frederick-analysis/2025-2026/.conda/envs/shiny_shared

# Run this to initialize conda if you have not done so
conda init --all

# Activate the environment
conda activate $PROJECT/x-cis220051/corporate/frederick-analysis/2025-2026/.conda/envs/shiny_shared

# Use the conda-env-mod script to create the module. Here --local-py means we will use the local python
conda-env-mod module -p $PROJECT/corporate/frederick-analysis/2025-2026/.conda/envs/shiny_shared -m $PROJECT/corporate/frederick-analysis/2025-2026/etc/modules --local-py
```

This will give an output:
```bash
+-----------------------------------------------------------------------------------------------+
| To use this environment, load the following modules:                                          |
|     module use $PROJECT/corporate/frederick-analysis/2025-2026/etc/modules |
|     module load conda-env/shiny_shared-py3.12.8                                               |
| (then standard 'conda install' / 'pip install' / run scripts)                                 |
+-----------------------------------------------------------------------------------------------+
```


# Appendix B. Another way to activate the shared environment

This is an alternative way to **2.2 Load the shared conda environment for shiny**

In Terminal, run these commands instead to get the same environment:

```bash
module load conda
conda activate $PROJECT/corporate/frederick-analysis/2025-2026/.conda/envs/shiny_shared
```

Run "conda init --all" if it appears in your terminal. 

This is the standard way to activate a shiny environment. It will attach a prefix "(shiny_shared)" in your terminal, while the one in **Step 2** does not. Whichever way is fine. You may run "conda info --envs" to see the list of all available conda environments and your current status.