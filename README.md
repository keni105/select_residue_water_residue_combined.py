# select_residue_water_residue_combined.py

This tool is particularly useful for studying water-mediated interactions in protein complexes and macromolecular structures.

The select_residue_water_residue_combined plugin for PyMOL is designed to:

1. Identify polar contacts and hydrogen bonds between residues in two specified domains via bridging water molecules.
2. Calculate distances between residues and water molecules, as well as the total distance across the interaction.
3. Visualize the interactions in PyMOL with different colors for polar contacts and hydrogen bonds.

# This plugin works only with the incentive version of PyMOL (tested in version 3.1.1) and not with the open-source PyMOL.

# Installation
	1.	Download the file: select_residue_water_residue_combined.py.
	2.	Open PyMOL.
	3.	Go to Plugin > Plugin Manager.
	4.	Click Install and select the file select_residue_water_residue_combined.py.
	5.	Restart PyMOL. The plugin is now available in the plugin menu.

 # Usage Instructions
 1. Preparing Your Protein Structure
    
    Load Your Protein Structure
    
    Define Domain Selections
    
      Use valid PyMOL selection syntax to define domain1 and domain2.

      Examples:
    
    a. Using Chains:
    
        select domain1, chain A
    
        select domain2, chain B
    b. Using Residue Ranges:

        select domain1, resi 1-100 and chain A
        select domain2, resi 101-200 and chain B
    c. Using Specific Residues:
    
        select domain1, chain A and resi 45,50,55
        select domain2, chain B and resi 60,65,70
 3. Running the Plugin
    
        select_residue_water_residue_combined domain1, domain2
    Replace domain1 and domain2 with the names of your domain selections.
    
 5. Output
    
    # Console Output
    
    The plugin prints detailed information about the interactions, including:
    
	  a. Water Molecule Counts:

        Number of water molecules near domainA: 352
        Number of water molecules near domainB: 361
        Number of bridging water molecules: 60
    b. Total Interactions Found:

        Total interactions: 15
    c. Interaction List:

        Residue model1/A`45/LYS interacts with residue model1/B`102/GLU via water id 567 (Hydrogen Bond). Distances: 2.80 Å (resA-water), 3.00 Å (water-resB), Total distance: 5.80 Å
        Residue model1/A`50/SER interacts with residue model1/B`110/THR via water id 570 (Polar Contact). Distances: 3.10 Å (resA-water), 3.50 Å (water-resB), Total distance: 6.60 Å
        ...

    # Visualizations

    Interacting Residues:
    
    Domain 1 Residues:

        Displayed as yellow sticks.

        Stored in the selection: interacting_residues_domain1.
    
	Domain 2 Residues:

	     Displayed as magenta sticks.

	     Stored in the selection interacting_residues_domain2.

	Interacting Water Molecules:

	    Displayed as cyan spheres with a sphere size of 0.25.

	    Stored in the selection interacting_waters.

	Hydrogen Bonds:

	    Displayed as green dashed lines.

	    Dash width set to 2.0.

	    Distance labels are displayed on the dashes.

	Polar Contacts:

	    Displayed as orange dashed lines.

	    Dash width set to 1.5.

	    Distance labels are displayed on the dashes.

Example of Result:
<img width="1512" alt="Screenshot 2024-11-30 at 2 58 14" src="https://github.com/user-attachments/assets/92700ab3-3335-48cb-ac5a-b4d2c8ee15c4">
<img width="868" alt="Screenshot 2024-11-30 at 3 01 15" src="https://github.com/user-attachments/assets/87822a0b-03b0-44a5-a111-4cd5c8946d99">
