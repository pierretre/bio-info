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


def assemble(reads):
    main_seq = reads[0]
    while(len(reads) > 1): 
        
        max_overlap = 0
        r = -1
        reads2 = [main_seq]
        
        for i in range(1, len(reads)):

            readi = reads[i]
            min_length = min(len(main_seq), len(readi))

            # On compare le début de la séquence principale avec la fin de la séquence reads[i] :
            #              MAINSEQ --> 
            #          <-- READI
            readi_slice = readi[len(readi) - min_length + 1:]
            main_seq_slice = main_seq[:min_length - 1]

            while(len(readi_slice) > 0 and readi_slice != main_seq_slice):
                readi_slice = readi_slice[1:]
                main_seq_slice = main_seq_slice[:-1]
                

            if(len(readi_slice) > abs(max_overlap)):
                max_overlap = -len(readi_slice)
                if(r != -1): 
                    reads2.append(reads[r])
                r = i
            
            if(abs(max_overlap) == min_length - 1):
                reads2.extend(reads[i+1:])
                break

            # On compare la fin de la séquence principale avec le début de la séquence reads[i] :
            #           <-- MAINSEQ
            #                 READI --> 
            readi_slice = readi[:min_length - 1]
            main_seq_slice = main_seq[len(main_seq) - min_length + 1:]

            while(len(readi_slice) > 0 and readi_slice != main_seq_slice):
                readi_slice = readi_slice[:-1]
                main_seq_slice = main_seq_slice[1:]

            if(len(readi_slice) > abs(max_overlap)):
                if(r != i and r!=-1): 
                    reads2.append(reads[r])
                max_overlap = len(readi_slice)
                r = i
            
            if(abs(max_overlap) == min_length - 1):
                reads2.extend(reads[i+1:])
                break

            if(len(readi_slice) < abs(max_overlap) and r != i):
                reads2.append(readi)

        if(max_overlap != 0):
            r_seq = reads[r]
            if(max_overlap < 0):
                main_seq = r_seq[:len(r_seq) + max_overlap] + main_seq
            else:
                main_seq = main_seq + r_seq[max_overlap:]
        else:
            # print(main_seq)
            # print(reads)
            # raise Exception("Erreur, au moins une séquence ne contient pas de chevauchement (de longueur > 0) dans les chaînes proposées")
            print("Erreur : au moins une séquence ne contient pas de chevauchement (de longueur > 0) dans les chaînes proposées")
            return main_seq

        reads = reads2
    return main_seq


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
assemble(coroReads100)
print("--- %s seconds ---" % (time.time() - start_time)) 

# Test pour 500 reads :
start_time = time.time()
assemble(coroReads500)
print("--- %s seconds ---" % (time.time() - start_time)) 

# Test pour 1000 reads :
start_time = time.time()
assemble(coroReads1000)
print("--- %s seconds ---" % (time.time() - start_time)) 

# Test pour 5000 reads :
start_time = time.time()
assemble(coroReads5000)
print("--- %s seconds ---" % (time.time() - start_time)) 

# Test pour 10000 reads :
start_time = time.time()
assemble(coroReads10000)
print("--- %s seconds ---" % (time.time() - start_time)) 

# Test pour tous les reads :
start_time = time.time()
assemble(coroReads)
print("--- %s seconds ---" % (time.time() - start_time)) 


