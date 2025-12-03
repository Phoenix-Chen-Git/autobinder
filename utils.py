import pymol
from pymol import cmd
import os
import subprocess
import re
import yaml
import shutil
import sys
import pandas as pd

def get_mutable_residues(pdb_path, target_chain, binder_chain, hotspots, distance_threshold):
    """
    Identifies mutable positions on the binder chain that are within a certain distance
    of the specified hotspots on the target chain using PyMOL.

    Args:
        pdb_path (str): Path to the PDB file.
        target_chain (str): Chain ID of the receptor (e.g., 'A').
        binder_chain (str): Chain ID of the binder (e.g., 'B').
        hotspots (list): List of residue indices (integers) on the target chain.
        distance_threshold (float): Distance threshold in Angstroms.

    Returns:
        list: Sorted list of unique residue indices (integers) on the binder chain.
    """
    # Initialize PyMOL in headless mode (if not already running)
    # This prevents the GUI from popping up if running locally
    pymol.finish_launching(['pymol', '-qc'])
    
    # Ensure a clean state
    cmd.reinitialize()
    
    # Load the PDB structure
    # We give it the object name 'structure'
    cmd.load(pdb_path, "structure")
    
    # Format hotspots for selection (e.g., [25, 26] -> "25+26")
    hotspots_str = "+.join(map(str, hotspots))"
    
    # 1. Select the target hotspot residues
    target_sel_str = f"chain {target_chain} and resi {hotspots_str}"
    
    # 2. Select binder residues that are within the distance threshold of the target hotspots
    # Logic: 
    #   (target_sel_str) around distance_threshold -> Selects all atoms within distance of target
    #   ... and chain binder_chain -> Restricts those atoms to the binder chain
    #   byres (...) -> Expands the selection to the full residues of those atoms
    selection_query = f"byres (chain {binder_chain} and ({target_sel_str} around {distance_threshold}))"
    
    # Create a named selection for iteration
    cmd.select("mutable_candidates", selection_query)
    
    # Set to store unique residue numbers
    mutable_residues = set()
    
    # Define a callback to collect residue numbers
    # 'resi' in PyMOL is a string, so we convert to int
    def collect_residues(resi):
        mutable_residues.add(int(resi))
        
    # Iterate over the selection. 
    # we use space={'collect_residues': collect_residues} to pass the function scope
    cmd.iterate("mutable_candidates", "collect_residues(resi)", space={"collect_residues": collect_residues})
    
    # Clean up PyMOL session
    cmd.delete("all")
    
    return sorted(list(mutable_residues))

def saturated_mutation(pdb_path, chain_id, mutation_sites, output_dir, current_gen_num, parent_id_in_gen):
    """
    Applies saturated mutagenesis on the given chain's specified sites.
    For each site, 20 PDB files are generated (one for each amino acid).
    Filenames are structured as Gen{N}_P{ParentID}_M{OLD}{site}_to_{NEW}.pdb for lineage tracking.
    
    Args:
        pdb_path (str): Path to the input PDB file (parent structure).
        chain_id (str): Chain ID to mutate.
        mutation_sites (list): List of residue indices (integers) to mutate.
        output_dir (str): Directory where generated PDBs will be saved.
        current_gen_num (int): The current generation number (for filename).
        parent_id_in_gen (int): Unique ID of the parent within its generation (for filename).
        
    Returns:
        list: List of tuples (path_to_generated_pdb, mutation_string).
    """
    # Initialize PyMOL (headless)
    pymol.finish_launching(['pymol', '-qc'])
    
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    generated_files_and_mutations = []
    
    # Standard Amino Acids (3-letter codes)
    amino_acids = [
        "ALA", "ARG", "ASN", "ASP", "CYS", 
        "GLN", "GLU", "GLY", "HIS", "ILE", 
        "LEU", "LYS", "MET", "PHE", "PRO", 
        "SER", "THR", "TRP", "TYR", "VAL"
    ]
    
    # The base filename prefix will include generation and parent ID
    base_filename_prefix = f"Gen{current_gen_num}_P{parent_id_in_gen}"
    
    for site in mutation_sites:
        cmd.reinitialize()
        cmd.load(pdb_path, "structure")
        
        stored = {'resn': ''}
        cmd.iterate(f"chain {chain_id} and resi {site} and name CA", "stored['resn'] = resn", space={'stored': stored})
        original_resn = stored['resn']
        
        if not original_resn:
            print(f"Warning: Could not find residue {site} in chain {chain_id} for {os.path.basename(pdb_path)}. Skipping mutation at this site.")
            cmd.delete("all")
            continue
            
        for aa in amino_acids:
            cmd.reinitialize()
            cmd.load(pdb_path, "structure")
            
            selection = f"chain {chain_id} and resi {site}"
            
            cmd.wizard("mutagenesis")
            cmd.do("refresh_wizard")
            
            cmd.get_wizard().set_mode(aa)
            cmd.get_wizard().do_select(selection)
            cmd.get_wizard().apply()
            
            cmd.set_wizard()
            
            # Construct the filename: Gen{N}_P{ParentID}_M{OLD}{site}_to_{NEW}.pdb
            mutation_string = f"{original_resn}{site}_to_{aa}"
            new_filename = f"{base_filename_prefix}_M{mutation_string}.pdb"
            output_path = os.path.join(output_dir, new_filename)
            
            cmd.save(output_path, "structure")
            generated_files_and_mutations.append((output_path, mutation_string))
            
    cmd.delete("all")
    return generated_files_and_mutations


def align_and_get_tmscore(pdb_file1, pdb_file2):
    """
    Aligns two PDB files using USalign and extracts the TM-score.

    Args:
        pdb_file1 (str): Path to the first PDB file.
        pdb_file2 (str): Path to the second PDB file.

    Returns:
        float: The TM-score of the alignment, or None if an error occurs.
    """
    command = ["USalign", pdb_file1, pdb_file2]
    # print(f"Executing command: {' '.join(command)}") # Debug print (removed for production cleanliness)
    
    try:
        # Run USalign command
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        # print(f"USalign stdout:\n{result.stdout}") # Debug print
        output_lines = result.stdout.splitlines()
        
        tm_score = None
        tm_score_pattern = re.compile(r"TM-score=\s*([0-9.]+)")
        for line in output_lines:
            match = tm_score_pattern.search(line)
            if match:
                tm_score = float(match.group(1))
                break
        
        return tm_score
        
    except subprocess.CalledProcessError as e:
        print(f"Error running USalign: {e}")
        print(f"Stderr: {e.stderr}")
        return None
    except FileNotFoundError:
        print("Error: USalign command not found. Make sure it's in your PATH.")
        return None
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        return None

def load_config(path="config.yaml"):
    """Loads the configuration from a YAML file."""
    if not os.path.exists(path):
        print(f"Error: Configuration file '{path}' not found.")
        sys.exit(1)
    with open(path, 'r') as f:
        return yaml.safe_load(f)

def run_afig(input_dir, output_base_dir):
    """
    Submits the AFig_custom.sh script to the cluster using sbatch.
    Uses --wait to block until the job completes.
    
    Args:
        input_dir (str): Directory containing mutant PDBs.
        output_base_dir (str): Base directory for output (script appends 'af2ig/').
    """
    script_path = "AFig_custom.sh"
    if not os.path.exists(script_path):
        print(f"Error: Script '{script_path}' not found.")
        return False

    # Check if sbatch is available
    if shutil.which("sbatch") is None:
        print("Error: 'sbatch' command not found. Cannot submit job to cluster.")
        print("       (Are you running this on the cluster login node?)")
        return False

    # The script usage: -pdbdir $1 -outpdbdir \"$2\"af2ig/
    
    command = ["sbatch", "--wait", script_path, input_dir, output_base_dir]
    
    print(f"Submitting AFig job: {' '.join(command)}")
    try:
        # capture_output=True allows us to see "Submitted batch job ..." if needed
        result = subprocess.run(command, check=True, text=True, capture_output=True)
        print(f"Job finished. Output:\n{result.stdout}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Error executing sbatch: {e}")
        print(f"Stderr: {e.stderr}")
        return False
    except Exception as e:
        print(f"Unexpected error executing sbatch: {e}")
        return False

def parse_score_file(score_file_path):
    """
    Parses the score file generated by AFig.
    Assumes a whitespace-separated format with 'SCORE:' prefix on each data line.
    Expected columns: 'description' (for name) and 'pae_interaction' (for iPae metric).
    """
    if not os.path.exists(score_file_path):
        print(f"Warning: Score file '{score_file_path}' not found.")
        return pd.DataFrame()

    try:
        # Read the file, skipping lines that don't start with "SCORE:"
        # and removing "SCORE:" prefix
        cleaned_lines = []
        with open(score_file_path, 'r') as f:
            for line in f:
                if line.startswith("SCORE:"):
                    cleaned_lines.append(line[len("SCORE:"):].strip())
        
        # Read into DataFrame, using first cleaned line as header
        from io import StringIO
        df = pd.read_csv(StringIO("\n".join(cleaned_lines)), delim_whitespace=True)
            
        # Normalize column names to lowercase
        df.columns = [c.lower() for c in df.columns]
        
        # Ensure required columns exist
        if 'description' not in df.columns or 'pae_interaction' not in df.columns:
            print(f"Error: Expected 'description' and 'pae_interaction' columns not found in score file. Found: {df.columns.tolist()}")
            return pd.DataFrame()

        # Rename for consistency if needed, but for now we use them directly
        # df.rename(columns={'pae_interaction': 'ipae'}, inplace=True)
        
        return df
    except Exception as e:
        print(f"Error parsing score file '{score_file_path}': {e}")
        return pd.DataFrame()
