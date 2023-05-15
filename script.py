"""
    Project for BIF/2023L course: Ancestry analysis

    :author :Sabina Gulcikova
            :(xgulci00@stud.fit.vutbr.cz)
"""
import os
import argparse
import pandas as pd
from Bio import Phylo
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
        self.tree = next(Phylo.parse(input_path, "newick"))

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
        self.phylo_tree = phylo_tree.tree
        self.leaf_nodes = self.phylo_tree.get_terminals()

        self.probabilities = probabilities.data
        self.alignments = alignments

        self.nodes_values = dict()
        self.nodes_ids = []
        for node in self.phylo_tree.get_nonterminals():
            self.nodes_ids.append(node.confidence)
        for node in self.phylo_tree.get_terminals():
            self.nodes_ids.append(node.name)

    def get_sequences_without_spaces(self):
        for node in self.phylo_tree.find_clades(terminal=False):
            node_id = node.confidence
            if node_id in self.nodes_values:
                continue
            else:
                self.nodes_values[node_id] = self.obtain_ml_sequence(node_id)

    def obtain_ml_sequence(self, node_id: str):
        result = ''
        extracted_values = self.probabilities.loc[self.probabilities['node'] == node_id] \
            .drop(['node'], axis=1) \
            .sort_values(by=['position']) \
            .drop(['position'], axis=1)

        for index, row in extracted_values.iterrows():
            row = pd.to_numeric(row, errors='coerce')
            max_column = row.idxmax()
            result = result + max_column

        return result

    def add_spaces_to_sequences(self):
        # Initialize mask for calculating weights which will determine if
        # corresponding aminoacid should be replaced by space
        weighted_space_mask = {node_id: [0] * 96 for node_id in self.nodes_ids}

        terminals = self.phylo_tree.get_terminals()

        # Parse the phylogenetic tree from a file
        elements = self.phylo_tree.get_nonterminals(order='level')
        elements.reverse()

        for node in elements:
            child_nodes = node.clades
            alignment_id = node.confidence
            curr_alignment = self.nodes_values[alignment_id]
            for id_, symb in enumerate(curr_alignment):
                weight = 0
                for child in child_nodes:
                    if child in terminals:
                        for instance in self.alignments:
                            if instance.name == child.name:
                                child_align = instance.alignment
                        weight += weighted_space_mask[child.name][id_]
                    else:
                        # Add the child's weight of necessity do add sequence to
                        # the weight of the current node, so that we do not have to traverse
                        # the child nodes again in future iterations
                        child_align = self.nodes_values[child.confidence]
                        weight += weighted_space_mask[child.confidence][id_]

                    if child_align[id_] == '-':
                        weight = weight + child.branch_length
                    else:
                        weight = weight - child.branch_length
                weighted_space_mask[node.confidence][id_] += weight
                if weight > 0:
                    curr_alignment = curr_alignment[:id_] + '-' + curr_alignment[id_ + 1:]
            self.nodes_values[alignment_id] = curr_alignment

    def print_outputs(self):
        os.makedirs("out", exist_ok=True)
        for id_, value in self.nodes_values.items():
            filename = os.path.join("out", "node_" + str(id_) + ".fas")
            value = '>' + str(id_) + '\n' + value[:60] + '\n' + value[60:]
            with open(filename, "w") as file:
                file.write(value)


if __name__ == '__main__':
    args = parse_args()

    ''' Parse input data from data/files '''
    phylo_tree = PhylogeneticTree(args.phylo_tree_path)
    prob_table = ProbabilityTable(args.ancestrals_path)
    prob_table.replace_missing_probabilities()

    ''' Process the MSAs from input fasta file by associating alignment with its name'''
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
    analyzer.add_spaces_to_sequences()

    ''' Display outputs '''
    analyzer.print_outputs()
