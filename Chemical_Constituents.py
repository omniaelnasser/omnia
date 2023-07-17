def chemical_constituents(file):
    """
        this function used to determine the chemical Constituents of a different proteins sets to be in a file
        file : is the desired file to calculate the chemical constituent of the proteins
        the output file is a file contain all 20 amino acids percentages in tab delimited text format
        :return:

        creator : Ahmed Mohamed Aziza
        Date : 21 Jun 2023
    """
    try:                                                                    # Try is used to work only on files
        from Bio.SeqUtils.ProtParam import ProteinAnalysis                  # Import Protein analysis functions to do the analysis
        from Bio import SeqIO                                               # Import SeqIO to parse the sequences from the desired file
        amino_ = {'A': ('A', 'ALA', 'alanine'),                             # Dictionary with all the abbreviations of amino acids
                  'R': ('R', 'ARG', 'arginine'),
                  'N': ('N', 'ASN', 'asparagine'),
                  'D': ('D', 'ASP', 'aspartic acid'),
                  'C': ('C', 'CYS', 'cysteine'),
                  'Q': ('Q', 'GLN', 'glutamine'),
                  'E': ('E', 'GLU', 'glutamic acid'),
                  'G': ('G', 'GLY', 'glycine'),
                  'H': ('H', 'HIS', 'histidine'),
                  'I': ('I', 'ILE', 'isoleucine'),
                  'L': ('L', 'LEU', 'leucine'),
                  'K': ('K', 'LYS', 'lysine'),
                  'M': ('M', 'MET', 'methionine'),
                  'F': ('F', 'PHE', 'phenylalanine'),
                  'P': ('P', 'PRO', 'proline'),
                  'S': ('S', 'SER', 'serine'),
                  'T': ('T', 'THR', 'threonine'),
                  'W': ('W', 'TRP', 'tryptophan'),
                  'Y': ('Y', 'TYR', 'tyrosine'),
                  'V': ('V', 'VAL', 'valine')}
        amino_acids = list(amino_.keys())                                  # Make a list of all the Keys
        result_path = './Chemical_Constituents.txt'                        # Specify The affinity File output Location by relative path
        result = open(result_path, 'w')                                    # Open the output file in write mode
        result.write(
            "Sequence_ID\tA\tR\tN\tD\tC\tQ\tE\tG\tH\tI\tL\tK\tM\tF\tP\tS\tT\tW\tY\tV\n")   # The header line in file
        for seq_record in SeqIO.parse(file, 'fasta'):                      # Parsing on all the sequences in the file
            x = ProteinAnalysis(str(seq_record))                           # Select the sequence record
            result.write('%s' % seq_record.id)                             # Write the Sequence ID in the result file
            for i in amino_acids:                                          # Iterating on the amino acids
                y = x.count_amino_acids()[i]                               # Calculate the Percentage of amino acid
                result.write("\t%s" % y)                                   # Write the result in the text file
            result.write("\n The Amino Acids: \n")                         # Add a header for the Dictionary to abbreviation
        for i in amino_acids:
            q = amino_[i
            result.write('%s\n' % str(q))
        result.close()
        print('The Process is successful\nThe out put file name is Chemical_Constituents.txt\n')
    except FileNotFoundError:                                              # Terminate the program if the file doesn't exist
        print("File not found\nPlease provide the desired file containing sequences with their IDs")
        exit()
