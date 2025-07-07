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

*   Access detailed spatiotemporal transcriptional landscapes across neuronal subtypes.
*   Explore gene expression patterns and LR interactions at any developmental stage.
*   Predict and analyze LR interactions between different neuronal cell types.

## Citation

**If you use scLRSomatoDev in your research, please cite our [paper](https://www.biorxiv.org/content/10.1101/2024.09.02.610245v3):**

Rémi Mathieu, Léa Corbières, Tangra Draia-Nicolau, Annousha Govindan, Vianney Bensa, Emilie Pallesi-Pocachard, Lucas Silvagnoli, Alfonso Represa, Carlos Cardoso, Ludovic Telley, Antoine de Chevigny. ***Inferring Ligand-Receptor Interactions between neuronal subtypes during mouse cortical development.*** *bioRxiv*.

<br>

### 🚀 Try a lite version of scLRSomatoDev Online!

Explore a lite, feature-focused version of our app directly in your browser—no installation required.

1.  Go to [sclrsomatodev.online](http://sclrsomatodev.online/) and click the **"Try the Lite version!"** button.
2.  **Please be patient while the app loads.** This process can take several minutes, and your browser may temporarily display an error message.
3.  Once the "Overview" page appears, you're ready to start exploring!

> [!NOTE]
> We are currently hosted on a server with limited resources. We appreciate your understanding as we work to improve performance.

<br>

## Getting Started

### System Requirements

*   **Processor**: 8+ cores (Intel/AMD)
*   **Memory**: 32GB RAM minimum
*   **Storage**: 50GB+ free space
*   **Display**: 1920 x 1080 resolution or higher
*   **Operating System**: Windows 10/11 (64-bit), macOS 10.15+, or Ubuntu 20.04+

### Installation

1.  **Get the project files**:

    You can get the files using one of the following methods:

    *   **Option A: Download ZIP**
        *   Navigate to the "Code" tab on the GitHub repository page.
        *   Select "Download ZIP".
        *   Extract the contents of the ZIP file to your desired location.

    *   **Option B: Git Clone**
        *   Open your terminal and run the following command:
            ```bash
            git clone https://github.com/Cortical-interactome/scLRSomatoDev.git
            ```

2.  **Download the required data**:
    *   Visit our [Zenodo repository](https://github.com/Cortical-interactome/scLRSomatoDev/) (TODO: Finish Zenodo Repo + add URL).
    *   Download the `data.zip` file.
    *   Extract its contents into your `scLRSomatoDev` folder.
    *   Verify that your folder structure matches the following:

    <p align="center">
    <br>
      <img width="600" src="images/Folder_structure.png" alt="Folder_structure">
    </p>

## How to Run

### Option 1: Docker (Recommended)

#### Prerequisites

*   **Docker**: If you don't have Docker installed, please [download and install it](https://www.docker.com/products/docker-desktop/) before proceeding.

#### 1. Build the Docker Image

Open your terminal, navigate to the project's root directory (the one containing the `Dockerfile`), and run:

**On Linux, Windows and macOS (Intel/AMD):**
```bash
docker build -t sclrshiny .
```

**On macOS (Apple Silicon):**
```bash
docker buildx build --platform linux/amd64 -t sclrshiny .
```

The initial build may take around 20 minutes to download and install the required R packages. Subsequent builds will be much faster.

#### 2. Run the Docker Container

To create a local container from the image you just built, run the appropriate command for your system from the project's root directory:

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

This command mounts your local `scLRSomatoDev/Data` directory to the `/app/Data` directory inside the container, which is where the app expects to find the data.

#### 3. Launch the App

After running the command, your terminal will display:
```
Listening on http://0.0.0.0:3838
```
The app will be available at http://localhost:3838. It may take a few moments to load.

> [!NOTE]
> You only need to **build the image once**. For future use, you can simply **run the container** and open the link in your web browser.

### Option 2: RStudio

#### Prerequisites
*   **R**: The app has been tested with R versions 4.2 and 4.3. If you don't have R, please [download and install it](https://pbil.univ-lyon1.fr/CRAN/) first.
*   **RStudio**: If RStudio is not installed, please [download and install it](https://posit.co/download/rstudio-desktop/) before proceeding.

#### Setup Instructions

#### 1. Create a conda environment using the provided environment file:

```{bash}
conda env create -f environment_scLRSomatoDev.yml
```

#### 2. Activate the environment:

```{bash}
conda activate r_env
```

#### 3. Open RStudio and set the working directory to the scLRSomatoDev folder.

**TODO : ADD FULL COMMAND FOR MISSING PACKAGES + ADD RUN indication**

#### 4. Install the required R packages if they are not already installed:

Since R can be tricky when managing dependencies, we have included a list of packages in the `environment_scLRSomatoDev.yml` file.

* If you encounter an error, you may need to install the missing R packages. You can do this by running the following commands in the console:

for a single specific package:

```r
install.packages("package_name")
```
for multiple packages:

```r
install.packages(c("package_name1", "package_name2", "package_name3"))
```
some packages might not be available on CRAN, you can install them from Bioconductor or github using the following command:

for Bioconductor packages:

```r
BiocManager::install("package_name")
```

for github packages:

```r
remotes::install_github("username/repo")
```

We are aware that R packages management and dependencies can be difficult. If you cannot find the package you need, you can ask for help on the [issues](https://github.com/Cortical-interactome/scLRSomatoDev/issues) page.

If you encounter any issues, you can ask for help on the [issues](https://github.com/Cortical-interactome/scLRSomatoDev/issues) page

#### 5. Run the app  :

To launch the application from RStudio, open either the `ui.R` or `server.R` file and click the "Run App" button located at the top of the file.

Alternatively, you can execute the following command in the R console:

```r
shiny::runApp('Path/to/your/App/file')
```

<br>

## Contact Us

*   **Technical Support & Bug Reports**
    *   Create an issue on our [GitHub repository](https://github.com/Cortical-interactome/scLRSomatoDev/issues).
    *   Please include detailed steps to reproduce the problem.

*   **General Inquiries**
    *   Email: [support@sclrsomatodev.online](mailto:support@sclrsomatodev.online)
    *   We typically respond within two business days.

> [!NOTE]
> For a faster response, please include your operating system, usage method (Online/RStudio/Docker), and any relevant error messages.
