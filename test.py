'''test if cavity in object change in Distance print self.cavity before and after calling getDescriptors()'''
from sequence.Sequence import Sequence
from Getdata.Cavity import Cavity
from shape.Shape import Shape
from Distance.Distance import Distance

protein = "C:/Users/32496/Desktop/not_3.2.1/structures/1A82.pdb"
fpocket = "C:/Users/32496/Desktop/not_3.2.1/fpocket/1A82_out"
pocket = 'pocket1_atm.pdb'
c = Cavity(protein, fpocket, pocket)
d = Distance(c.cavity, c.ligand)
