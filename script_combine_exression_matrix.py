'''Combines rows in two matricies based on having the same gene.

USAGE
script_combine_expression_matrix.py MATRIX_A MATRIX_B

MATRIX_A and MATRIX_B need to be tab-delimited matricies. The first line in each
   matrix is assumed to be a header, and all subsequent rows are data.
   The first column is assumed to contain the name of the genes in the row.

OUTPUT
A combined matrix named as MATRIX_A-MATRIX_B
'''

import sys

def main(matrixA, matrixB):
    # Load in the first matrix
    matrixA_dict = {}
    matrixA_file = open(matrixA, 'r')
    matrixA_header = matrixA_file.readline().strip()
    for line in matrixA_file:
        tab = line.strip().split('\t')
        gene = tab[0].upper()
        data = '\t'.join(tab[1:])
        matrixA_dict[gene] = data
    matrixA_file.close()
    # Go through matrixB and concatonate when you can
    
    output = open('%s-%s' % (matrixA, matrixB), 'w')
    matrixB_file = open(matrixB, 'r')
    matrixB_header = matrixB_file.readline().strip()
    matrixB_header = '\t'.join(matrixB_header.split('\t')[1:])
    output.write('%s\t%s\n' % (matrixA_header, matrixB_header))
    for line in matrixB_file:
        tab = line.strip().split('\t')
        gene = tab[0].upper()
        data = '\t'.join(tab[1:])
        print (gene)
        try:
            matrixA_data = matrixA_dict[gene]
        except KeyError:
            print ('Could not find gene in first matrix: ', gene)
        else:
            outline = '%s\t%s\t%s\n' % (gene, matrixA_data, data)
            output.write(outline)
    output.close()
    matrixB_file.close()

if __name__ == "__main__":
    if len(sys.argv) <= 2:
        print (__doc__)
        sys.exit()
    else:
        main(sys.argv[1], sys.argv[2])
