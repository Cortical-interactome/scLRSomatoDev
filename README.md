
# scLRSomatoDev

<p align="center">
  <br>
    <img width="600" src="images/logo_scLRSomatoDev.svg" alt="logo_scLRSomatoDev">
 <br>
 <br>
</p>

## Introduction

Lorem ipsum dolor sit amet, consectetur adipiscing elit. Sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.

Furthermore, an open-access journal article describing scLRSomatoDev has been published in Science, available [here](test).

## Run scLRSomatoDev

### Prerequisites

* **You will need to have Docker installed and running**. If you don't have Docker, please [download and install it](https://www.docker.com/products/docker-desktop/) before continuing.
  
* Download the files from this github repository.
    + To Download the project files, navigate to this GitHub repository and click on the "Code" tab, then select "Download ZIP".
    + Once downloaded, extract the contents of the .zip
    
* Download the required Data available on [Zenodo](inc).
    + Ddownload the data.zip and "extract here" within the the scLRSomatoDev folder.

* Before started files in your scLRSomatoDev folder shoud look like this :

  <p align="center">
  <br>
    <img width="600" src="Folder_structure.png" alt="Folder_structure">
</p>
  

### Building the Docker image

To build the docker image, open your terminal application in the directory where the copied GitHub folder is located (the root must be the folder containing scLRSomatoDev folder and the Dockerfile) and run:

```{bash}
docker build -t sclrshiny .
```
The image will take about 50 min to build.

### Running the Docker image

To create a local container from the previous built image run:

```{bash}
docker run --rm -p 3838:3838 -v "/path/to/folder:/app" sclrshiny
```
You have to replace "path/to/folder" by the path where the Data and www folders are located. For example if the Data and www folders are located in "/home/user/Documents/Data_scLRSomatoDev", you have to run:

```{bash}
docker run --rm -p 3838:3838 -v "/home/user/Documents/Data_scLRSomatoDev:/app" sclrshiny
```
### Launching the app in a web browser

On your terminal application, it will be displayed:

```{bash}
Listening on http://0.0.0.0:3838
```

Copy the link written after "Listening on" on a web browser. The shiny app scLRSomatoDev should be displayed after few minutes.

> [!NOTE]
>  You do not have to go through all these steps each time to launch the scLRSomatoDev shiny App. **You have to build the image only the first time**.
>
> For all subsequent times, you just have to **run the Docker image and copy the link to your web browser**.

<br>
<br>
<br>
