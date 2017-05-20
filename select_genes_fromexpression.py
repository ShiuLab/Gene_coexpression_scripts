import sys

gene_list = sys.argv[1]
expression_file = sys.argv[2]
a = sys.argv[3]
expression = open(expression_file, 'r')


line = expression.readline()
#print line
output = open('selected_gene_expression_%s' % a, 'w')
output.write('%s' % line)

#n = 5208

while line:
    info = line.strip().split()
    list = open(gene_list, 'r')
    list_line = list.readline()
    while list_line:
        gene_info = list_line.strip('\n').strip('\r').strip().split()
        gene = gene_info[0]
        go = info[0]
        go = go.replace('"', '')
        #print go
        if gene == go:
            output.write('%s' % line)
        list_line = list.readline()
    list.close()
    line = expression.readline()
    