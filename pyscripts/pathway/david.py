"""
Processes DAVID outputs into a format that is parsable into R iGraph packages.
"""
import os
import re


def split_cluster(file):
    '''
    split_cluster(file): splits the native DAVID output file into separate
                         files for processing each cluster by itself.
    
    @param file: the input file that can be downloaded from DAVID output
    
    '''
    with open(file,'r') as infile:
        for line in infile:
            if line.startswith("Annotation Cluster"):
                F = open(line[:20]+".txt", 'w')
            else:
                F.write(line)


def adj_matrix():
    '''
    adj_matrix(): processes each file created by split_cluster() and creates 
              adjacency matrices to be processed by the R packages iGraph
              or RedeR.
    '''
    clusters = []
    files = []
    
    for file in os.listdir(os.getcwd()):
        if re.match("Annotation Cluster \d.txt",file):
            print("Processing {0}".format(file))
            with open(file, 'r') as infile:
                targets = 0 # number of targets included in each vertex (pathway)
                for line in infile:
                    if line[0] != "\n":
                        line = line.split("\t")
                        if line[5] == ("Genes"): # if not the header
                            outfile = open(file.replace(".txt",".adj"), 'w')
                            files.append(file.replace(".txt",".adj"))
                        else:
                            targets = targets + 1
                            outfile.write(line[5].replace(', ','\t')+"\n")
                clusters.append(targets)
    outfile.close()
    
    for z in range(0,len(clusters)):
        count = [] # n x n adjacency matrix that counts the weights of the shared targets between two nodes
        for b in range(0, clusters[z]): # initialize the count[]
            new = []
            for j in range(0, clusters[z]):
                new.append(0)
            count.append(new)
        i = 0
        with open(files[z],'r') as infile:
            vertices = [] # one line for each vertex (pathway)
            for line in infile:
                line = line.replace('\n','').split('\t')
                vertices.append(line) # reads and inputs each target to its corresponding pathway
            for i in range(0,len(vertices)): # for each pathway
                for j in range(0,len(vertices[i])): # for each target
                    term = vertices[i][j] # each target in one pathway
                    for k in range(0, len(vertices)): # check if it exists within other pathways. The more matches, the higher the weight.
                        if term in vertices[k]:
                            count[i][k] = count[i][k] + 1
        infile.close()
        with open(files[z],'w') as outfile:
            for i in range(0, clusters[z]):
                for j in range(0,len(count)):
                    outfile.write("{0}\t".format(count[i][j]))
                outfile.write("\n")


def main():
    split_cluster("files/david.txt")
    adj_matrix()
    
if __name__ == '__main__':
    main()