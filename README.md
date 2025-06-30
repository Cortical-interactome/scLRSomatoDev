# scLRSomatoDev

<p align="center">
  <br>
    <img width="600" src="images/logo_scLRSomatoDev.svg" alt="logo_scLRSomatoDev">
  <br>
  <br>
</p>

## Introduction

scLRSomatoDev is a powerful Shiny application designed to explore and analyze ligand-receptor (LR) interactions in cortical development. This tool represents a significant advancement in understanding the molecular mechanisms underlying brain circuit assembly through large-scale single-cell transcriptomic analyses. With its intuitive design, scLRSomatoDev enables seamless exploration of complex biological data.

### Key Features

* Access detailed spatiotemporal transcriptional landscapes across neuronal subtypes
* Explore gene expression patterns and LR interactions at any developmental stage
* Predict and analyze LR interactions between different neuronal cell types

### Citation

**If you use scLRSomatoDev in your research, please cite our [paper](https://www.biorxiv.org/content/10.1101/2024.09.02.610245v3):** 
Rémi Mathieu, Léa Corbières, Tangra Draia-Nicolau, Annousha Govindan, Vianney Bensa, Emilie Pallesi-Pocachard, Lucas Silvagnoli, Alfonso Represa, Carlos Cardoso, Ludovic Telley, Antoine de Chevigny. ***Inferring Ligand-Receptor Interactions between neuronal subtypes during mouse cortical development.*** <em>bioRxiv</em>

### Try the Lite Version Online (TODO UPDATE THIS MESSAGE)

For discovery purpose, we have set up a dedicated server running scLRSomatoDev-lite, an lightened version of the app with only a subset of the features,  that you can access online from anywhere in the world.

1. Visit our website: [sclrsomatodev.online](http://sclrsomatodev.online/)

2. Click on the "Try the Lite version!" button

3. Please read the following important message:

> [!CAUTION]
> Please note that the application page may take a few minutes to load. Your browser may display an error page at first ("no response from server"), but the app will load after a few moments. We are currently hosted on a server with limited resources, and appreciate your patience.

4. Wait for the server to load the app (this process can take several minutes as mentioned above)

5. Once the app is loaded and the "Overview" page appears, you're ready to start using scLRSomatoDev!

## Table of Contents

1. [Getting Started](#getting-started)
   * [Minimum System Requirements](#minimum-system-requirements)
   * [Download Required Files](#download-required-files)
2. [Run scLRSomatoDev](#run-sclrsomatodev)
   * [Using Docker](#using-docker)
   * [Using RStudio](#using-rstudio)
  
4. [Contact us](#contact-us)

## Prerequisites

### Minimum System Requirements

  * Processor: 8+ cores (Intel/AMD)
  * Memory: 32GB RAM (64GB recommended for Docker)
  * Storage: 50GB+ free space
  * Display: 1920 x 1080 resolution or Higher
  * Operating System: Windows 10/11 (64-bit), macOS 10.15+, or Ubuntu 20.04+

### Download Required Files

1. **Download the project files**:
   * Click on the "Code" tab
   * Select "Download ZIP"
   * Extract the contents of the ZIP file to your desired location

2. **Download the required data**:
   * Visit our [Zenodo repository](https://github.com/Cortical-interactome/scLRSomatoDev/) (TODO: Finish Zenodo Repo + add URL)
   * Download the `data.zip` file
   * Extract the contents into your scLRSomatoDev folder
   * Verify your folder structure matches the following :

    <p align="center">
    <br>
      <img width="600" src="images/Folder_structure.png" alt="Folder_structure">
    </p>

## Run scLRSomatoDev

### Using Docker (Recommanded)

#### Prerequisites

* **Docker**: If you don't have Docker installed, please [download and install it](https://www.docker.com/products/docker-desktop/) before proceeding.

### Building the Docker image

To build the docker image, open your terminal application in the directory where the copied GitHub folder is located (the root must be the folder containing scLRSomatoDev folder and the Dockerfile) and run:

```bash
docker build -t sclrshiny .
```
The initial image build will take some time as it needs to download and install all the required R packages (about 20 min). Subsequent builds will be much faster.

### Running the Docker container

To create a local container from the previous built image run:

To create a local container from the image you just built, run the appropriate command for your system from the root directory of the project:

**On Linux or macOS:**
```bash
docker run --rm -p 3838:3838 -v "$(pwd)/scLRSomatoDev/Data:/app/Data" sclrshiny
```

**On Windows (Command Prompt):**
```bash
docker run --rm -p 3838:3838 -v "%cd%\scLRSomatoDev\Data:/app/Data" sclrshiny
```

**On Windows (PowerShell):**
```bash
docker run --rm -p 3838:3838 -v "${PWD}\scLRSomatoDev\Data:/app/Data" sclrshiny
```

This command mounts your local `scLRSomatoDev/Data` directory into the `/app/Data` directory inside the container, where the Shiny app expects to find it.

### Launching the app in a web browser

After running the command, your terminal will display:

```bash
Listening on http://0.0.0.0:3838
```

The app will be up and running at http://localhost:3838. The shiny app scLRSomatoDev should be displayed after few minutes.
> [!NOTE]
> You do not have to go through all these steps each time to launch the scLRSomatoDev shiny app. **You have to build the image only the first time**.
>
> For all subsequent times, you just have to **run the Docker image and copy the link to your web browser**.

### Using RStudio

#### Prerequisites
* **R** : The app has been tested with R versions 4.3.2 and 4.2. If you do not have R installed, please [download and install it](https://pbil.univ-lyon1.fr/CRAN/) before proceeding.
* **RStudio**: If RStudio is not installed, please [download and install it](https://posit.co/download/rstudio-desktop/) before proceeding.
* **Conda**: Ensure Conda is installed on your system. If not, please [download and install Anaconda Or Miniconda](https://www.anaconda.com/download/success) before proceeding.

#### Setup Instructions

1. Create a conda environment using the provided environment file:
    ```{bash}
    conda env create -f environment_scLRSomatoDev.yml
    ```

2. Activate the environment:
    ```{bash}
    conda activate r_env
    ```

3. Open RStudio and set the working directory to the scLRSomatoDev folder.

4. Install the required R packages if they are not already installed:
    ```{r}
    install.packages(c("shiny", "tidyverse", "Seurat", "ggplot2", "dplyr", "tidyr"))
    ```
  
5. To run the Shiny app, click the 'Run' or 'Play' button in either the ui.r or server.r file or run the following command : 
    ```{r}
    shiny::runApp("app.R")
    ```

6. The app will launch either in your default web browser or in the RStudio dedicated window.

## Contact us

* **Technical Support & Bug Reports**
  * Create an issue on our [GitHub repository](https://github.com/yourusername/scLRSomatoDev/issues)
  * Include detailed steps to reproduce any problems

* **General Inquiries**
  * Email: [support@sclrsomatodev.online](mailto:support@sclrsomatodev.online)
  * Response time: Usually within 2 business days

> [!NOTE]
> For faster response, please include your operating system, Usage method (Online/RStudio/Docker), and any relevant error messages.


