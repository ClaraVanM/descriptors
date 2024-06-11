import sequence.amPseAAC as ampseAAC
import sequence.autocorrelation as autocorrelation
import sequence.composition as composition
import sequence.CTD as CTD
import sequence.pseAAC as pseAAC
import sequence.sequence_order as sequence_order


class Sequence:
    AA_symb = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
               'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
               'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
               'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    def __init__(self, seq_or_cavity):
        if isinstance(seq_or_cavity,str):
            self.seq = seq_or_cavity
        else:
            self.seq = self.getSequence(seq_or_cavity)


    @staticmethod
    def getSequence(cavity):
        """
        Calculate the amino acid sequence from structure

        :param cavity: pandas dataframe with protein structure data from pdb file
        :return: string with sequence of structure
        """
        cavity['AA_number'].astype(float)
        cavity = cavity.sort_values(by=['AA_number'], ascending=True)
        cavity = cavity.drop_duplicates(subset=["AA_number"])
        sequence = ''.join(cavity["AA"].map(Sequence.AA_symb))
        return sequence

    def sequence_descriptors(self):
        """

        :return: sequence-based descriptors
        """
        seq_descr = dict()
        seq_descr.update(composition.aa_composition(self.seq))
        seq_descr.update(composition.dipeptide_composition(self.seq))
        seq_descr.update(composition.tripeptide_composition(self.seq))
        seq_descr.update(composition.conjoint_triad(self.seq))
        CTD_comp = CTD.ctd_composition(self.seq)
        for values in CTD_comp.values():
            seq_descr.update(values)
        CTD_trans = CTD.ctd_transition(self.seq)
        for values in CTD_trans.values():
            seq_descr.update(values)
        CTD_distr = CTD.ctd_distribution(self.seq)
        for values in CTD_distr.values():
            seq_descr.update(values)
        seq_descr.update(pseAAC.pseaac(self.seq))
        seq_descr.update(ampseAAC.amp_pse_AAC(self.seq))
        moreau_broto = autocorrelation.autocorrelation(self.seq, autocorrelation.moreaubroto_ac)
        for values in moreau_broto.values():
            seq_descr.update(values)
        moran = autocorrelation.autocorrelation(self.seq, autocorrelation.moran_ac)
        for values in moran.values():
            seq_descr.update(values)
        geary = autocorrelation.autocorrelation(self.seq, autocorrelation.geary_ac)
        for values in geary.values():
            seq_descr.update(values)
        seq_descr.update(sequence_order.tau_qsoc(self.seq))
        return seq_descr

    def partial_descriptors(self):
        '''

        :return: partial set of sequence-based descriptors
        '''
        descriptors = dict()
        pseaac = pseAAC.pseaac(self.seq)
        moran = autocorrelation.autocorrelation(self.seq, autocorrelation.moran_ac)
        qsoc = sequence_order.tau_qsoc(self.seq)
        descriptors.update(pseaac)
        descriptors.update(qsoc)
        for values in moran.values():
            descriptors.update(values)
        return descriptors