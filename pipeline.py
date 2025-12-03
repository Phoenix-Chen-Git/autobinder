import os
import sys
import glob
import shutil
import re
from utils import get_mutable_residues, saturated_mutation, align_and_get_tmscore, load_config, run_afig, parse_score_file

def main():
    config = load_config()
    
    # Configuration
    target_chain = config.get('target_chain')
    binder_chain = config.get('binder_chain')
    hotspots = config.get('hotspots', [])
    distance_threshold = config.get('distance_threshold', 5.0)
    tm_score_threshold = config.get('tm_score_threshold', 0.9)
    max_iterations = config.get('max_iterations', 5)
    beam_width = config.get('beam_width', 5)
    keep_checkpoints = config.get('keep_checkpoints', True)

    print("--- AutoBinder Design Pipeline Started ---")
    print(f"Target Chain: {target_chain}, Binder Chain: {binder_chain}")
    print(f"Hotspots: {hotspots}")
    print(f"TM-score Threshold: {tm_score_threshold}, Beam Width: {beam_width}")
    
    # Initial Population from Iteration 0
    current_iter_dir = "iteration0"
    if not os.path.exists(current_iter_dir):
        print(f"Error: Initial directory '{current_iter_dir}' does not exist.")
        return

    initial_pdbs = glob.glob(os.path.join(current_iter_dir, "*.pdb"))
    if not initial_pdbs:
        print(f"Error: No PDB files found in '{current_iter_dir}'.")
        return

    # Population will be a list of dictionaries, each representing a parent structure:
    # {'path': '/path/to/pdb', 'id_in_gen': 0, 'mutation_string': 'WT' (or specific for mutants)}
    # For the initial population, we assign ID P0.
    population = []
    for idx, pdb_path in enumerate(initial_pdbs):
        population.append({
            'path': pdb_path,
            'id_in_gen': idx,
            'mutation_string': 'WT' # Wild Type for initial structure
        })
        
    print(f"Initial Population: {len(population)} structures.")

    for iteration in range(max_iterations):
        print(f"\n=== Iteration {iteration} ===")
        
        # Define directories for the current iteration
        mutants_dir = os.path.join(current_iter_dir, "mutants")
        af2ig_output_dir = os.path.join(current_iter_dir, "af2ig")
        af2ig_score_file = os.path.join(current_iter_dir, "af2igpae")
        
        # Ensure directories exist
        os.makedirs(mutants_dir, exist_ok=True)
        os.makedirs(af2ig_output_dir, exist_ok=True)
        
        # 1. Generate Mutants
        print(f"--- Step 1: Generating Mutants in {mutants_dir} ---")
        mutant_tracking_map = {} # Maps core_mutant_name (e.g., Gen0_P0_MVAL45_to_ALA) to its parent_path and mutation_string
        generated_count = 0
        
        for parent_item in population:
            parent_pdb_path = parent_item['path']
            parent_id_in_gen = parent_item['id_in_gen']
            
            print(f"Processing parent: {os.path.basename(parent_pdb_path)} (ID: P{parent_id_in_gen})")
            
            # Identify mutable residues
            mutable_residues = get_mutable_residues(
                parent_pdb_path, target_chain, binder_chain, hotspots, distance_threshold
            )
            
            if not mutable_residues:
                print(f"  No mutable residues found for {os.path.basename(parent_pdb_path)}.")
                continue
                
            print(f"  Mutable residues identified: {mutable_residues}")
            
            # Saturated Mutagenesis
            # Returns a list of (path_to_mutant_pdb, mutation_string)
            new_mutant_data = saturated_mutation(
                parent_pdb_path, binder_chain, mutable_residues, mutants_dir, iteration, parent_id_in_gen
            )
            
            for mutant_pdb_path, mutation_str in new_mutant_data:
                # The core mutant name (e.g., Gen0_P0_MVAL45_to_ALA)
                core_mutant_name = os.path.basename(mutant_pdb_path).replace(".pdb", "")
                
                mutant_tracking_map[core_mutant_name] = {
                    'parent_pdb_path': parent_pdb_path,
                    'mutation_string': mutation_str
                }
                
            generated_count += len(new_mutant_data)
            
        print(f"Total mutants generated: {generated_count}")
        
        if generated_count == 0:
            print("No mutants generated for this iteration. Stopping pipeline.")
            break
            
        # 2. Run Structure Prediction (AFig)
        print("--- Step 2: Running Structure Prediction (AFig) ---")
        # AFig output PDBs go to current_iter_dir/af2ig/
        # Score file is current_iter_dir/af2igpae
        
        output_base_arg = current_iter_dir + os.sep # Ensure trailing separator for script
        success = run_afig(mutants_dir, output_base_arg)
        
        if not success:
            print("Structure prediction failed. Check AFig logs/configuration.")
            break # Stop pipeline if AFig fails
            
        # 3. Analyze and Filter
        print("--- Step 3: Analyzing and Filtering ---")
        
        scores_df = parse_score_file(af2ig_score_file)
        if scores_df.empty:
            print("No scores found or failed to parse score file. Cannot rank. (Did AFig run and produce output?)")
            break
            
        name_col = 'description'
        ipae_col = 'pae_interaction'
        
        if name_col not in scores_df.columns or ipae_col not in scores_df.columns:
            print(f"Error: Expected '{name_col}' or '{ipae_col}' columns not found in score file. Found: {scores_df.columns.tolist()}")
            break
            
        candidates = []
        
        for index, row in scores_df.iterrows():
            pdb_description_from_afig = str(row[name_col]) # e.g., Gen0_P0_MVAL45_to_ALA_af2pred
            ipae = float(row[ipae_col])
            
            # Core name without '_af2pred' postfix and '.pdb' suffix for map lookup
            core_name_for_lookup = pdb_description_from_afig.replace("_af2pred", "")
            if core_name_for_lookup.endswith('.pdb'):
                core_name_for_lookup = core_name_for_lookup.replace('.pdb', '')

            mutant_details = mutant_tracking_map.get(core_name_for_lookup)
            
            if not mutant_details:
                print(f"Warning: Mutant details not found in tracking map for AFig output '{pdb_description_from_afig}'. Skipping.")
                continue

            parent_pdb_for_tmscore = mutant_details['parent_pdb_path']
            original_mutation_string = mutant_details['mutation_string']
            
            # Full path to the predicted PDB (e.g., iteration0/af2ig/Gen0_P0_MVAL45_to_ALA_af2pred.pdb)
            pred_pdb_path = os.path.join(af2ig_output_dir, pdb_description_from_afig + ".pdb") 
            
            if not os.path.exists(pred_pdb_path):
                print(f"Warning: Predicted PDB '{pred_pdb_path}' not found. Skipping.")
                continue
                
            # TM-score Filter
            tm_score = align_and_get_tmscore(parent_pdb_for_tmscore, pred_pdb_path)
            
            status = "Discard"
            if tm_score is not None:
                if tm_score >= tm_score_threshold:
                    status = "Keep"
                    candidates.append({
                        'path': pred_pdb_path, # Path to AFig's output PDB
                        'ipae': ipae,
                        'tm_score': tm_score,
                        'parent_pdb_path': parent_pdb_for_tmscore,
                        'mutation_string': original_mutation_string
                    })
                else:
                    status = f"Discard (TM-score {tm_score:.2f} < {tm_score_threshold})"
            else:
                status = "Discard (TM-score calculation failed)"
                
            print(f"  Candidate: {pdb_description_from_afig}.pdb, iPae: {ipae:.2f}, TM: {tm_score if tm_score is not None else 'N/A'}, Status: {status}")

        # 4. Selection
        print("--- Step 4: Selection ---")
        if not candidates:
            print("No candidates passed the filters for this iteration. Stopping pipeline.")
            break
            
        # Sort by iPae (Ascending - lower is better)
        candidates.sort(key=lambda x: x['ipae'])
        
        selected = candidates[:beam_width]
        print(f"Selected Top {len(selected)} candidates for next iteration.")
        
        # 5. Prepare Next Iteration
        print("--- Step 5: Preparing Next Iteration ---")
        next_iter_num = iteration + 1
        next_iter_dir = f"iteration{next_iter_num}"
        
        os.makedirs(next_iter_dir, exist_ok=True)
            
        next_population = []
        for idx, cand_item in enumerate(selected):
            new_id_in_gen = idx # New ID for this candidate within the *next* generation
            src_afig_output_path = cand_item['path'] # This is AFig's output PDB
            
            # Construct the new filename for the promoted PDB (Gen{N+1}_P{new_id}_M{mutation}.pdb)
            new_promoted_filename = f"Gen{next_iter_num}_P{new_id_in_gen}_M{cand_item['mutation_string']}.pdb"
            dst_next_iter_path = os.path.join(next_iter_dir, new_promoted_filename)
            
            shutil.copy2(src_afig_output_path, dst_next_iter_path)
            
            next_population.append({
                'path': dst_next_iter_path,
                'id_in_gen': new_id_in_gen,
                'mutation_string': cand_item['mutation_string']
            })
            print(f"  Promoted: {new_promoted_filename} (iPae: {cand_item['ipae']:.2f})")
            
        # Update current state for next loop
        current_iter_dir = next_iter_dir
        population = next_population
        
        # Optional Cleanup (if keep_checkpoints is False)
        if not keep_checkpoints:
            print(f"  Cleaning up intermediate files for iteration {iteration}...")
            shutil.rmtree(mutants_dir, ignore_errors=True)
            # Keep af2ig_output_dir for now, as score file references it.
            # If score file is consolidated, this could be cleaned too.
            # For now, only remove generated mutants.

    print("\n--- AutoBinder Design Pipeline Completed ---")

if __name__ == "__main__":
    main()
