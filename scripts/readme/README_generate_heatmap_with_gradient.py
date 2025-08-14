README_generate_heatmap_with_gradient.py

This Python script generates a 3D protein structure heatmap visualization of amino acid variation indices for the pilE protein in Neisseria gonorrhoeae. The script maps variation scores to a color gradient (green to red) and creates PyMol-compatible scripts for protein structure visualization.

The script requires two input files:
Variation data: rv_02_pilE_Pt_PubMLST_variation_score.xlsx (located in data/raw/ folder)
PDB structure file: pilE_FA1090_Pt_b470a_unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb (located in data/raw/ folder)

Files should contain amino acid positions and corresponding variation indices.

Install required missing Python packages (if any) using: pip install pandas biopython matplotlib

Analysis Steps
Loads amino acid variation index data from the Excel file
Normalizes variation indices to 0-100% scale for color mapping
Excludes stop codon position (position 165 (*)) from analysis
Maps normalized scores to RGB color gradient (green for low variation, red for high variation)
Validates amino acid positions against the PDB structure file
Extracts residue information from protein structure
Generates PyMol-compatible color mapping and visualization scripts

Output Files
Creates PyMol_output/ folder containing:
v09_5_variation_mapping.txt: Amino acid positions with variation indices and color codes
V09_5_heatmap_script.pml: PyMol script for loading structure and applying heatmap colors
v09_5_heatmap_legend.png: Color legend showing variation index scale
pilE_heatmap.pdb: Colored protein structure file
pilE_heatmap.png: High-resolution heatmap image (300 DPI)

  this script from Python environment.
