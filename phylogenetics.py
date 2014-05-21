"""
References: http://biopython.org/DIST/docs/tutorial/Tutorial.html and http://biopython.org/wiki/Phylo
Fasta sequences: http://pycogent.org/cookbook/building_a_tree_of_life.html
"""
import os

from Bio.Align.Applications import ClustalwCommandline
from Bio import Phylo

class PhyloTree(object):
    """
    Is a really simple object responsible for performing a multiple
    sequence alignment with sequences in a fasta file, loading a 
    phylogenetic tree from the results and displaying it.
    """
    def __init__(self, path=None):
        """
        Saves the fasta file path and the filename with extension removed.
        """
        self.path = path
        self.filename = os.path.splitext(path)[0]
        self.tree = None

    def align(self):
        """
        Runs the clustalw alg. on fasta sequences, which writes a dnd file.
        """
        msa_func = ClustalwCommandline('clustalw2', infile=self.path)
        msa_func()

    def display(self, isascii=False):
        """
        Loads the tree from a file. And displays it.
        """
        self.tree = Phylo.read('{}.dnd'.format(self.filename), 'newick')
        
        if isascii:
            Phylo.draw_ascii(self.tree)

        try:
            import pylab
            Phylo.draw_graphviz(self.tree)
            pylab.show()
        except:
            print('Warning: failed to display using graphviz')
            Phylo.draw_ascii(self.tree)


if __name__ == '__main__':
    ptree = PhyloTree('silva_sequences.fasta')
    ptree.align()
    ptree.display()

