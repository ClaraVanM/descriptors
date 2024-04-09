import pandas as pd


AA_symb = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


def AA_per_buriedness(df):
    df = df.drop_duplicates(subset=["AA_number"])
    total_count = {}
    for i in df["buriedness"].unique():
        df_temp = df[df["buriedness"]==i]
        aa_count = df_temp['AA'].value_counts()
        total_count[i] = {aa:aa_count.get(aa,0) for aa in AA_symb.keys()}
    return total_count


#take only nonzero narrowness depths into account or 1 0 to have the floor of the pocket
def exposed_aa(df, narrow_list):
    if not narrow_list.count(0.0) == 0:
        index = narrow_list.count(0.0)-1
    else: index = 0
    depth = list(range(index,5))
    df = df[df['buriedness'].isin(depth)]
    #choose closest point of aa to axis as representative for aa
    subset_index = df.groupby('AA')['dist_from_axis'].idxmin()
    df = df.loc[subset_index]
    # take closest 5 of aa to axes
    df = df.sort_values(by='dist_from_axis').head(int(len(df)*0.3))
    return df
