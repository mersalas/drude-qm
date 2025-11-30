# Python script to renumber residues continuously in a PDB file

import re

def renumber_residues(input_pdb, output_pdb):
    """
    Reads a PDB file, assigns continuous residue numbers starting from 1,
    and writes the updated PDB to output_pdb.
    """
    res_map = {}  # Mapping from original (chain, old_resnum) to new_resnum
    next_resnum = 1

    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if line.startswith(('ATOM  ', 'HETATM')):
                # Extract chain ID and original residue number
                chain_id = line[21]
                old_resnum = line[22:26].strip()
                key = (chain_id, old_resnum)
                
                # Assign new residue number if not already mapped
                if key not in res_map:
                    res_map[key] = next_resnum
                    next_resnum += 1
                
                new_resnum = res_map[key]
                
                # Reformat the residue number field (columns 23-26)
                new_res_str = f"{new_resnum:>4}"
                updated_line = line[:22] + new_res_str + line[26:]
                outfile.write(updated_line)
            else:
                # Copy other lines unchanged
                outfile.write(line)

if __name__ == "__main__":
    input_filename = "3_MRH30.pdb"   # Replace with your input filename
    output_filename = "output.pdb" # Desired output filename
    renumber_residues(input_filename, output_filename)
    print(f"Residues renumbered continuously. Output written to {output_filename}.")

