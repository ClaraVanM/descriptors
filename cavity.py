import pymolPy3


def get_cavity_atoms(protein_file, fpocket_out):
    """
    Selects all atoms/AA around the predicted cavity and puts them in new file
    :param protein_file: full protein.pdb file
    :param fpocket_out: output folder of fpocket after cavity prediction
    :return: None
    """
    pm = pymolPy3.pymolPy3(0)
    pm(f"load {protein_file}")
    pm(f"load {fpocket_out}/pockets/pocket1_atm.pdb")
    pm(f"select cavity, pocket1_atm")
    protein_name = protein_file.split("/")[-1].split(".")[0]
    pm(f"select prot, {protein_name}")
    pm(f"select neighborhood, prot near_to 10 of cavity")
    pm(f"save {protein_name}_neighbor.pdb, neighborhood")
    #remove HETATM = water + ligand
    with open(f"{protein_name}_neighbor.pdb", 'r') as input:
        with open(f"{protein_name}_cavity.pdb", "w") as output:
            for line in input:
                if line.startswith("HETATM"):
                    continue
                else:
                    output.write(line)
        output.close()
    input.close()
