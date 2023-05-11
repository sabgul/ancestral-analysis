# - nacitat viacnasobne zarovnanie z msa.fasta
# - nacitat strukturu fylogenetickeho stromu (Phylo lib z BioPython)
# - nacitat data z ancestrals.csv
# - na zaklade posteriornych pravdepodobnosti v ancestrals odhadnut najviac pravdepodobne
#   ancestralne sekvencie, vratane doplnenia ancestralnych medzier


# objekt na subor msa {array objektov na zarovnanie, metoda na parsovanie suboru}
# objekt na zaorvnanie -> name, sequence

# objekt na strom (obsahuje vsetky nody, root)
# objekt na node (meno, potomkovia, vzdialenost ku kazdemu potomkovi)
import argparse
import pandas as pd
from typing import List


def parse_args() -> argparse.Namespace:
    args = argparse.ArgumentParser()
    args.add_argument('--ancestrals_path', type=str, default='data/ancestrals.csv',
                      help='Path to ancestrals file of input data.')
    args.add_argument('--msa_path', type=str, default='data/msa.fasta',
                      help='Path to multiple alignments file of input.')
    args.add_argument('--phylo_tree_path', type=str, default='data/tree.tre',
                      help='Path to phylogenetic tree.')

    return args.parse_args()


class PhylogeneticTree:
    def __init__(self, input_path: str) -> None:
        self.input_path = input_path
        # self.tree_data =
        # TODO read phylo tree into the structure
        pass


class ProbabilityTable:
    def __init__(self, input_path: str) -> None:
        self.input_path = input_path
        self.probabilities = pd.read_csv(self.input_path)

    def replace_missing_probabilities(self) -> None:
        self.probabilities = self.probabilities.replace('-', 0.0)


class MultSeqAlignment:
    def __init__(self, name: str) -> None:
        self.name = name
        self.alignment = ''
        pass

    def get_alignment(self, line_one: str, line_two: str):
        # Strip both leading > and newline from the name of alignment
        self.name = self.name.rstrip('\n')
        self.name = self.name[1:]

        # Concatenate the sequence from subsequences split by fasta requirements
        self.alignment = line_one.rstrip('\n') + line_two.rstrip('\n')


class AncestryAnalyzer:
    def __init__(self,
                 phylo_tree: PhylogeneticTree,
                 probabilities: ProbabilityTable,
                 alignments: List[MultSeqAlignment]) -> None:
        pass


if __name__ == '__main__':
    args = parse_args()

    ''' Parse input data from data/files '''
    phylo_tree = PhylogeneticTree(args.phylo_tree_path)
    prob_table = ProbabilityTable(args.ancestrals_path)
    prob_table.replace_missing_probabilities()

    ''' Process the MSA from input fasta file by associating alignment with its name'''
    msa = []
    msa_file = open('data/msa.fasta', 'r')
    lines = msa_file.readlines()
    for i in range(0, len(lines), 3):
        temp_alignment = MultSeqAlignment(lines[i])
        temp_alignment.get_alignment(lines[i+1], lines[i+2])
        msa.append(temp_alignment)

    ''' Perform ancestry analysis '''
    analyzer = AncestryAnalyzer(phylo_tree, prob_table, msa)

    #     mame sekvencie na listoch, budeme vypocitavat sekvencie na vnutornych uzloch
    #       teda si najskor pre cely strom vytiahneme uzol, na zaklade tabulky vyextrahujem
    #        najpravdepodobnejsi aminokyselinu a skonkatenujem na danu poziciu

    # ked budem mat poskladany cely strom, budem sa pre kazdy uzol pozerat na jeho potomkov
    # ulozim si skore: space_replacement_weight = 0 -> budem prechadzat potomkami a vzidalenostami k nim
    # ak bude u potomka na danom mieste medzera -> space_replacement_weight + vzdialenost k nemu
    # ak u potomka na danom meiste medzera nebude -> -vzdialenost k nemu

    # ak bude space_replacement_weight > 0 -> vlozim medzeru, inak necham tak

    # vytvorim subor > node_id, na dalsi riadok hodim sekvenciu, pripadne pridam newline po xtom symbole
    # set_most_probable
    # append_spaces
    # print_outputs

    # x najskor si v ancestrals nahradis - za 0
    # najskor vyberies najpravdepodobnejsiu aminokyselinu, vytvoris sekvencie bez medzier
    # na zaklade potomkov sa pre kazdy stlpec pozries ci u potomkov prevazuju v danom stlpci medzery,
    #   ak ano, das tam medzeru, ak nie, nechas

    # vypises do osobitnych suborov


    # file1 = open('data/msa.fasta', 'r')
    # count = 0
    # Lines = file1.readlines()
    #
    # for line in Lines:
    #     # if(count % 3 == 1):
    #     count += 1
    #     print("Line: ", line, count)