import yaml #for writing yamls
import sys #for seeing if yaml is available 
import subprocess #allows for running bash in python
import os #for connecting to directories 

if 'yaml' in sys.modules:
    print(" Yaml Module already imported.")
else:
    print("Yaml Module not imported yet. Please insure the Yaml Module is installed")
    
# Reading data from the YAML file
with open('YourSamplesHere.yaml', 'r') as file:
    loaded_data = yaml.safe_load(file)

print("Data read from 'data.yaml':")
print(loaded_data)

def download_sra(accession, target_dir="initial_data"): #function to do SRA download
    command = ["prefetch", accession, "--output-directory", target_dir] #command is prefetch to get SRA + the accession
    if not os.path.exists(target_dir):
        os.makedirs(target_dir) 
    try: #try this section
        # Runs the command on the given accession
        result = subprocess.run(command, check=True, capture_output=True, text=True) 
        #check breaks is if no connection
        #capture_output means python will save what you get
        #text=true means python will decode bytes to text of the output
        print(f"Downloaded {accession}") #for comfort 
    except subprocess.CalledProcessError as e: #if this fails pass and throw error
        print(f"Error Downloading {accession}: {e.stderr}") #what failed, error



for a in loaded_data:
    download_sra(a)
    



