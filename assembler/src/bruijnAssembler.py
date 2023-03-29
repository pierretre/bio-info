import time

def readSequencesFromFastaFile(filepath):

    reads = []
    currentSeq = ""

    f = open(filepath, "r")
    for line in f:

        if(line[0] == '>'):
            if(currentSeq != ""):
                reads.append(currentSeq)
                currentSeq = ""
        else:
            currentSeq += line.strip("\r\n")
            
    if(currentSeq != ""):
        reads.append(currentSeq)

    f.close()
    return reads


def assemble(reads, k):
    nodes = {}
    
    k = round(k * len(reads[0])) if 0.5 < k < round((len(reads[0]) - 1) / len(reads[0])) else len(reads[0]) - 1

    # utilisation d'un set pour éviter les redondances de kmers
    kmers = set()
    for r in reads:

        while len(r) >= k:
            kmers.add(r[:k])
            r = r[1:]

    # Add a kmer to the graph 
    # nodes are represented with 2 attributes:
    # - next nodes (name)
    # - previous nodes (name)
    for kmer in kmers:
        key1 = kmer[:k - 1]
        key2 = kmer[1:]

        #First node :
        node1 = nodes.get(key1)
        if not node1:
            node1 = {'next': {}, 'previous': {}}
        
        if not node1['next'].get(key2) :
            node1['next'][key2] = 0
        node1['next'][key2] += 1

        # Add second node :
        node2 = nodes.get(key2)
        if not node2:
            node2 = {'next': {}, 'previous': {}}
                  
        if not node2['previous'].get(key1) :
            node2['previous'][key1] = 0
        node2['previous'][key1] += 1

        nodes[key1] = node1
        nodes[key2] = node2


    # assemble the graph :
    base_key = next(iter(nodes))
    node = nodes[base_key]
    sequence = base_key
    
    # run through nodes in both left and right directions
    while len(node['previous']) > 0:
        key = next(iter(node['previous']))
        node = nodes.pop(key)
        sequence = key[:1] + sequence

    if len(node['previous']) > 1:
        raise Exception("Erreur, au moins une séquence ne contient pas de chevauchement (de longueur > 0) dans les chaînes proposées")

    node = nodes.pop(base_key)
    while len(node['next']) > 0:
        key = next(iter(node['next']))
        node = nodes.pop(key)
        sequence += key[-1:]

    if len(node['next']) > 1:
        raise Exception("Erreur, au moins une séquence ne contient pas de chevauchement (de longueur > 0) dans les chaînes proposées")
    
    if len(nodes) > 0:
        raise Exception("Erreur, au moins une séquence ne contient pas de chevauchement (de longueur > 0) dans les chaînes proposées")

    return sequence


# Lecture des reads fasta du Coronavirus :
PATH_TO_FASTA = './fasta/'

coroReads100 = readSequencesFromFastaFile(PATH_TO_FASTA+'MN908947.3_fraction0000100_READS_MIXED.fasta')
coroReads500 = readSequencesFromFastaFile(PATH_TO_FASTA+'MN908947.3_fraction0000500_READS_MIXED.fasta')
coroReads1000 = readSequencesFromFastaFile(PATH_TO_FASTA+'MN908947.3_fraction0001000_READS_MIXED.fasta')
coroReads5000 = readSequencesFromFastaFile(PATH_TO_FASTA+'MN908947.3_fraction0005000_READS_MIXED.fasta')
coroReads10000 = readSequencesFromFastaFile(PATH_TO_FASTA+'MN908947.3_fraction0010000_READS_MIXED.fasta')
coroReads = readSequencesFromFastaFile(PATH_TO_FASTA+'MN908947.3_READS_MIXED.fasta')

# Tests temps de calcul :

# Test pour 100 reads :
start_time = time.time()
assemble(coroReads100, 0)
print("--- %s seconds ---" % (time.time() - start_time)) 

# Test pour 500 reads :
start_time = time.time()
assemble(coroReads500, 0)
print("--- %s seconds ---" % (time.time() - start_time)) 

# Test pour 1000 reads :
start_time = time.time()
assemble(coroReads1000, 0)
print("--- %s seconds ---" % (time.time() - start_time)) 

# Test pour 5000 reads :
start_time = time.time()
assemble(coroReads5000, 0)
print("--- %s seconds ---" % (time.time() - start_time)) 

# Test pour 10000 reads :
start_time = time.time()
assemble(coroReads10000, 0)
print("--- %s seconds ---" % (time.time() - start_time)) 

# Test pour tous les reads :
start_time = time.time()
assemble(coroReads, 0)
print("--- %s seconds ---" % (time.time() - start_time)) 



