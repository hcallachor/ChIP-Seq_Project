import subprocess #allows for running bash in python
import os 

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


list_accessions = ['SRR057262', 'SRR057268', 'SRR6055333']
for a in list_accessions:
    download_sra(a)


