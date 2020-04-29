Confluence Machine (BAM Algorithm).
This is the repo for the BAM algorithm part of the confluence machine. 
To run the code;
1. create two folders; data and output_files and grant them read and write permissions. If you skip this, your code will only excute the first stage. We are working on making sure that in future versions, you don't have to manually create these two folders
2. make sure that you have docker installed(https://docs.docker.com/v17.09/engine/installation/)and check that it is working.Docker doesn't work directly on Windows. Run the command docker --version in the powershell or cmd
3. The original netcdf files should be in the "original_data" folder. 
4. Then simply run the docker_driver.sh script. Note that the first time you run it, you need a stable internet
connection as the images will be downloaded locally for future use.The end to end machine takes atleast 10-20 minutes to run, depending on your processing power and the number of initial results.  
