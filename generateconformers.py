import sys
from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter

def generate_conformers(input_files, output_files, num_conformers=10, energy_window=5.0):
    for input_sdf, output_sdf in zip(input_files, output_files):
        # Read molecules from the current input SDF file
        suppl = Chem.SDMolSupplier(input_sdf)
        molecules = [mol for mol in suppl if mol is not None]

        # Initialize an SDWriter for the current output file
        writer = SDWriter(output_sdf)

        # Process each molecule in the current input file
        for mol in molecules:
            mol = Chem.AddHs(mol)  # Add hydrogens

            # Generate conformers
            conformer_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=AllChem.ETKDGv3())

            # Optimize each conformer and write to the output file
            for conf_id in conformer_ids:
                AllChem.UFFOptimizeMolecule(mol, confId=conf_id)
                writer.write(mol, confId=conf_id)
        
        # Close the writer for the current output file
        writer.close()
        print(f"Conformers for {input_sdf} written to {output_sdf}.")

if __name__ == "__main__":
    # Ensure correct number of arguments are passed
    if len(sys.argv) < 3:
        print("Usage: python generate_conformers.py <input_files> <output_files>")
        sys.exit(1)

    # Get input and output file paths from command line arguments
    input_files = sys.argv[1].split(",")  # e.g., "input1.sdf,input2.sdf"
    output_files = sys.argv[2].split(",")  # e.g., "output1.sdf,output2.sdf"

    # Call the function to generate conformers
    generate_conformers(input_files, output_files)
