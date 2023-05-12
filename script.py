
import argparse
import pandas as pd
from typing import List
from Bio import Phylo


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
        self.tree = next(Phylo.parse(input_path, "newick"))
        # TODO read phylo tree into the structure

    def display_tree(self) -> None:
        Phylo.draw_ascii(self.tree)


class ProbabilityTable:
    def __init__(self, input_path: str) -> None:
        self.input_path = input_path
        self.data = pd.read_csv(self.input_path)

    def replace_missing_probabilities(self) -> None:
        self.data = self.data.replace('-', 0)


class MultSeqAlignment:
    def __init__(self, name: str) -> None:
        self.name = name
        self.alignment = ''
        pass

    def get_alignment(self, line_one: str, line_two: str) -> None:
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
        # iterate through the tree
        self.phylo_tree = phylo_tree
        self.leaf_nodes = self.phylo_tree.tree.get_terminals()

        self.probabilities = probabilities.data
        self.alignments = alignments

        self.nodes_values = dict()

    def get_sequences_without_spaces(self):
        for node in self.phylo_tree.tree.find_clades(terminal=False):
            # TODO check if name is at confidence or name attribute
            node_id = node.confidence
            if node_id in self.nodes_values:
                continue
            else:
                alignment = self.obtain_ml_sequence(node_id)
                self.nodes_values[node_id] = alignment
        # for leaf_node in self.leaf_nodes:
        #     # current_node = leaf_node
        #     current_node = leaf_node.ancestor
        #     ancestors = current_node.ancestors()
        #     for ancestor in ancestors:
        #
        #     while current_node is not None:
        #         # Perform operations on the current node
        #         # skip duplicates
        #         node_id = current_node.confidence
        #         if node_id in self.nodes_values:
        #             continue
        #         else:
        #             alignment = self.obtain_ml_sequence(node_id)
        #             self.nodes_values[node_id] = alignment
        #         print(current_node)
        #
        #         current_node = current_node.parent

    def obtain_ml_sequence(self, node_id: str):
        result = ''
        # df.loc[df['column_name'] == some_value]
        extracted_values = self.probabilities.loc[self.probabilities['node'] == node_id]
        extracted_values = extracted_values.drop(['node'], axis=1)
        extracted_values = extracted_values.sort_values(by=['position'])
        extracted_values = extracted_values.drop(['position'], axis=1)

        for index, row in extracted_values.iterrows():
            row = pd.to_numeric(row, errors='coerce')
            max_column = row.idxmax()
            result = result + max_column

        # TODO calculate alignment, based on the probabilities df
        # extract rows whose first column == node_id
        # df[df['Column1'].isin(['B', 'D'])]
        # related_rows = self.probabilities[self.probabilities['node'].isin[node_id]]
        # related_rows = self.probabilities['nodes'].isin[node_id]
        # vyfiltrovat hodnoty kt maju spravnu hodnotu nodes, vymrdat tento stlpec prec
        # sort by 'position', vymrdat tento stlpec prec
        # potom pojdeme riadok po riadku, zistime na ktorom key je najvyssia hodnota,
        #   tento key appendujeme do result
        return result

    def add_spaces_to_sequences(self):
        # prechod stromom od deti vyssie
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
    analyzer.get_sequences_without_spaces()
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