# Example for running a GLM on Tiled AI Data for AI4AD Consortium
This example requires access to tiled data available to members of the AI4AD consortium. If you require access, please let the folks at Curii know.

If you have access, data can be downloaded from here:
https://workbench.2xpu4.arvadosapi.com/projects/2xpu4-j7d0g-50i3pwf3a16ubbf

This example has a docker container to provide the resources to run the example.

To build the docker container, change to Dockerfile directory and run:
docker build -t example .

Below is the syntax on how to run a docker container (in this case called example) with 2 mounted volumes (using -v flag) containing code and downloaded data. You will have to replace these with correct paths for your directories.

docker run -v /YOURPATH/AI4AD:/AI4AD -v /YOURPATH/data:/data -ti example

To run the example:
* Adjust the file paths to reflect where you have your downloaded tile data (e.g., Xdata_file, Xrdata_file, etc)
* Run Rscript TiledData_Example.r
