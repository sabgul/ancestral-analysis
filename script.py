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


def parse_args() -> argparse.Namespace:
    args = argparse.ArgumentParser()
    args.add_argument('--ancestrals-csv-path', type=str, default='ancestrals.csv',
                      help='Path to ancestrals file of input data.')
    args.add_argument('--msa-fasta-path', type=str, default='msa.fasta',
                      help='Path to multiple alignments file of input.')
    args.add_argument('--tree-path', type=str, default='tree.tre',
                      help='Path to phylogenetic tree.')

    return args.parse_args()

class PhylogeneticTree:
    def __init__(self, input_path: str) -> None:
        self.input_path = input_path
        pass


if __name__ == '__main__':
    args = parse_args()

    phylo_tree = PhylogeneticTree(args.tree-path)

    # file1 = open('msa.fasta', 'r')
    # Lines = file1.readlines()
    # for line in Lines:
    #     count += 1
    #     print("Line{}: {}".format(count, line.strip()))