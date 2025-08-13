import pandas as pd
from Bio.PDB import PDBParser, Polypeptide
from Bio.SeqUtils import seq1
import os
import matplotlib.pyplot as plt

# Input Files
variation_file = "./rv_02_pilE_Pt_PubMLST_variation_score.xlsx"
pdb_file = "./pilE_FA1090_Pt_b470a_unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb"
output_dir = "./PyMol_output"

# Ensure Output Directory Exists
os.makedirs(output_dir, exist_ok=True)
output_mapping_file = os.path.join(output_dir, "v09_5_variation_mapping.txt")
output_pymol_script = os.path.join(output_dir, "V09_5_heatmap_script.pml")
legend_output_file = os.path.join(output_dir, "v09_5_heatmap_legend.png")

# Load Variation Index Data
variation_data = pd.read_excel(variation_file)

# Normalize Variation Index (0-100%) for Gradient Mapping
variation_data["Normalized_Index"] = variation_data["Variation_Index"] / variation_data["Variation_Index"].max()

# Ignore stop codon positions
stop_codon_position = 165  # Define the stop codon position
variation_data = variation_data[variation_data["Amino_Acid_Position"] != stop_codon_position]

# Map Normalized Scores to PyMol-Compatible RGB Colors (Green to Red)
def normalized_to_rgb(normalized_value):
    # Green: #006400; Red: #FF0000
    red = normalized_value  # Proportional red intensity
    green = 1 - normalized_value  # Proportional green intensity
    blue = 0  # No blue for this gradient
    return f"[{red:.2f}, {green:.2f}, {blue:.2f}]"

variation_data["Color"] = variation_data["Normalized_Index"].apply(normalized_to_rgb)

# Validate Positions with PDB File
parser = PDBParser(QUIET=True)
structure = parser.get_structure("pilE", pdb_file)

# Extract Residue IDs and Names from the PDB File
pdb_residues = {residue.id[1]: seq1(residue.resname)  # Use seq1 for conversion
                for model in structure
                for chain in model
                for residue in chain
                if Polypeptide.is_aa(residue)}

# Check for Missing Residues
xlsx_positions = set(variation_data["Amino_Acid_Position"])
missing_positions = xlsx_positions - set(pdb_residues)

if missing_positions:
    print(f"Warning: The following positions are in the XLSX file but not in the PDB file: {missing_positions}")

# Merge Residue Names into the Data
variation_data["Amino_Acid_Name"] = variation_data["Amino_Acid_Position"].map(pdb_residues)

# Save Variation Mapping for PyMol and Legend
with open(output_mapping_file, "w") as f:
    f.write("Amino_Acid_Position\tAmino_Acid_Name\tAmino_Acid_Variation_Index\tColor\n")
    for _, row in variation_data.iterrows():
        f.write(f"{row['Amino_Acid_Position']}\t{row['Amino_Acid_Name']}\t{row['Variation_Index']}\t{row['Color']}\n")
print(f"Variation mapping saved to {output_mapping_file}")

# Generate PyMol Script for Heatmap
with open(output_pymol_script, "w") as script:
    script.write(f"load {pdb_file}, pilE\n")
    for _, row in variation_data.iterrows():
        position = row["Amino_Acid_Position"]
        color = row["Color"]
        script.write(f"set_color color{position}, {color}\n")
        script.write(f"color color{position}, resi {position}\n")
    script.write(f"save {os.path.join(output_dir, 'pilE_heatmap.pdb')}\n")
    script.write(f"png {os.path.join(output_dir, 'pilE_heatmap.png')}, dpi=300, ray=1\n")
print(f"PyMol script saved to {output_pymol_script}")

# Generate Heatmap Legend
# Adjustable Parameters for Legend
x_label_size = 12
tick_label_size = 10
legend_height = 1.5  # Reduced height for minimal white space

# Create figure with minimal margins
plt.figure(figsize=(6, legend_height))
plt.subplots_adjust(left=0.03, right=0.97, top=0.97, bottom=0.23)  # Tighter margins  # Adjusted margins

# Generate exact color gradient for legend
gradient = [[(i / 255, (255 - i) / 255, 0) for i in range(256)]]
plt.imshow(gradient, aspect="auto", extent=[0, 100, 0, 1])

# Configure Tick Style and Appearance
plt.tick_params(axis="x", direction="out", colors="black", labelsize=tick_label_size)  # External ticks
plt.yticks([])  # Remove y-axis ticks
plt.gca().spines[:].set_color("white")  # Set frame to white

# Add Labels
plt.xlabel("PilE Amino Acids Variation Index, %", fontsize=x_label_size, color="black")
plt.xticks(fontsize=tick_label_size, color="black")

# Save with reduced DPI and tight layout
plt.savefig(legend_output_file, dpi=300, bbox_inches="tight")
print(f"Legend saved to {legend_output_file}")
