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

print("Data read from 'YourSamplesHere.yaml':")
print(loaded_data)

# This dictionary will hold the data for our output YAML
# You can update the Reference Genome string as needed
output_results = {
    "SRAs": [],
    "Single End": [],
    "Paired End": [],
    "Reference Genome": ["GCF_000002765.6"]
}

def get_library_layout(accession, target_dir):
    try:
        # -X 1: Only look at the first spot
        # -Z: Print to stdout (no files created)
        # --split-spot: If it's paired, it splits 1 spot into 2 reads
        cmd = ["fastq-dump", "-X", "1", "-Z", "--split-spot", accession]
        
        # Run the command from the target directory where the data was downloaded
        result = subprocess.run(cmd, capture_output=True, text=True, check=True, cwd=target_dir)
        
        # A single-end FASTQ record has 4 lines.
        # A paired-end FASTQ record (split) will have 8 lines.
        lines = result.stdout.strip().splitlines()
        line_count = len(lines)
        
        if line_count == 8:
            return "PAIRED"
        elif line_count == 4:
            return "SINGLE"
        elif line_count > 0 and line_count % 4 == 0:
            return "PAIRED"
        else:
            return "UNKNOWN"
            
    except Exception as e:
        # If fastq-dump also fails, we check the directory directly for help
        print(f"DEBUG: Detection failed for {accession}: {e}")
        return "UNKNOWN"

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

        #this calls get_library_layout to the download code
        layout = get_library_layout(accession, target_dir)
        print(f"Layout detected for {accession}: {layout}")
        
        # Add results dictionary for the final YAML as paired or single
        output_results["SRAs"].append(accession)
        if layout == "PAIRED":
            output_results["Paired End"].append(accession)
        elif layout == "SINGLE":
            output_results["Single End"].append(accession)

    except subprocess.CalledProcessError as e: #if this fails pass and throw error
        print(f"Error Downloading {accession}: {e.stderr}") #what failed, error

# Main execution loop
for a in loaded_data:
    download_sra(a)

# Writing the results out to a new YAML file
with open('CompProjectconfig.yaml', 'w') as outfile:
    #no sorting keys to maintain order made by single/paired
    yaml.dump(output_results, outfile, default_flow_style=False, sort_keys=False)

print("\nProcessing complete. Layout results written to 'CompProjectconfig.yaml'.")