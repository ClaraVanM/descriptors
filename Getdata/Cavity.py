from Getdata.Structures import Structures
import pymolPy3

class Cavity(Structures):
    def __init__(self,pdbfile, fpocket_out, pocket):
        cavity_file = self.get_cavity(pdbfile, fpocket_out, pocket)
        super().__init__(cavity_file)
        self.cavity = self.extract_cavity()
        self.ligand = self.extract_ligand()


    @staticmethod
    def get_neighborhood(pdbfile, fpocket_out, pocket):
        """
        Selects all atoms/AA around the predicted cavity and puts them in new file
        :param protein_file: full protein.pdb file
        :param fpocket_out: output folder of fpocket after cavity prediction
        :return: None
        """
        pm = pymolPy3.pymolPy3(0)
        pm(f"load {pdbfile}")
        pm(f"load {fpocket_out}/pockets/{pocket}")
        pm(f"select cavity, {pocket.split('.')[0]}")
        protein_name = pdbfile.split("/")[-1].split(".")[0]
        pm(f"select prot, {protein_name}")
        pm(f"select neighborhood, prot near_to 10 of cavity")
        pm(f"save {protein_name}_neighbor.pdb, neighborhood")
        return protein_name + '_neighbor.pdb'

    @staticmethod
    def get_cavity(pdbfile, fpocket_out, pocket):
        file = Cavity.get_neighborhood(pdbfile, fpocket_out, pocket)
        return file

    def extract_cavity(self):
        return self.structures[self.structures['type'] == 'ATOM']

    def extract_ligand(self):
        return self.structures[self.structures['type'] == 'HETATM']

