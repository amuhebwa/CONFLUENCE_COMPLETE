Confluence Machine (BAM Algorithm).
This is the repo for the BAM algorithm part of the confluence machine. 
To run the code;
1. make sure that you have docker installed(https://docs.docker.com/v17.09/engine/installation/)and check that it is working.Docker doesn't work directly on Windows. Run the command docker --version in the powershell or cmd
2. The original netcdf files should be in the "original_data" folder. 
3. Then simply run the docker_driver.sh script. Note that the first time you run it, you need a stable internet
connection as the images will be downloaded locally for future use.The end to end machine takes atleast 10-20 minutes to run, depending on your processing power and the number of initial results.  
