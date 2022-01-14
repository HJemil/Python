 #A function to count all nucleotides
def get_counts(seq, base=3):
    counts = {}
    for i in range(0,len(seq),base):
        slice = seq[i:i+base]
        if slice in counts:
            counts[slice ] += 1
        else:
            counts[slice ] = 1
    return counts

# A function to translate codons into aminoacids
def translate_dna(dna_seq): 
    bases = ["T","C","A","G"]
    codons = [a+b+c for a in bases for b in bases for c in bases]
    amino_acids = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
    codon_table = dict(zip(codons, amino_acids))
    protein_seq = []
    
    for i in range(0,len(dna_seq)-3,3):
        codon = dna_seq[i:i+3]
        if codon in codon_table:
            protein_seq +=codon_table[codon]
        else:
            protein_seq += "X"
    return ''.join(protein_seq)

# Function to get CG content 
def get_cg_content (gene_seq):
    at = gene_seq.count("A") + gene_seq.count("T")
    cg = gene_seq.count("C") + gene_seq.count("G")
    if ((at+cg) > 0):
        gc_content = cg/(at+cg)
    else:
        gc_content = 0
    return gc_content

# Function to get start and stop codon
def get_start_stop_codon (use_genome):
    import matplotlib.pyplot as plt
    start_codon = {}
    stop_codon = {}
    for seq_key in use_genome:
        tmp_seq = use_genome[seq_key]["sequence"]
        tmp_start = tmp_seq[0:3]
        if(tmp_start not in start_codon.keys()):
            start_codon[tmp_start] = 1
        else:
            start_codon[tmp_start] += 1

        tmp_stop = tmp_seq[(len(tmp_seq)-3):(len(tmp_seq)+1)]
        if(tmp_stop not in stop_codon.keys()):
            stop_codon[tmp_stop] = 1
        else:
            stop_codon[tmp_stop] += 1
            
    plt.subplots(ncols=2,dpi=150)
    plt.bar(start_codon.keys(),start_codon.values())
    plt.bar(stop_codon.keys(),stop_codon.values())
    plt.xticks(rotation="vertical",fontsize=6) ;
    plt.show()
    
  # function to get genome features

def show_genome_features(use_genome,
                        show_no_of_genes = True,
                        show_average_gene_length = True,
                        show_average_coding_length = True,
                        plot_gene_length_distri = True,
                        plot_coding_length_distri = True,
                        plot_coding_vs_gene_length = True,
                        plot_gc_content_distri = True,
                        plot_gc_vs_coding_len = True):
    
    import matplotlib.pyplot as plt
    
    # extract features
        # gene lengths
    gene_len = []
    for seq_key in use_genome:
        tmp = [use_genome[seq_key]["stop"]-use_genome[seq_key]["start"]]
        gene_len += tmp
        # coding lengths
    coding_len = []
    for seq_key in use_genome:
        coding_len += [len(use_genome[seq_key]["sequence"])]    
        # gc content
    gc_per_gene = []    
    for seq_key in use_genome:
        gc_per_gene += [get_cg_content(use_genome[seq_key]["sequence"])]
    
    
    # show no of genes
    if(show_no_of_genes):
        print("This genome consists of " + str(len(use_genome))  + " genes")
    
    # show average gene length
    if(show_average_gene_length):
        print("The average gene length is " + str(sum(gene_len)/len(gene_len)) + " basepairs")   
        
    # show average coding length    
    if(show_average_coding_length):
        print("The average coding region is " + str(sum(coding_len)/len(coding_len)) + " basepairs")   
    
    # plot gene length distribution
    if(plot_gene_length_distri):
        plt.hist(gene_len,bins=1000)
        plt.title("Gene Length")
        plt.xlim(0,20000)
        plt.show()
     
    # plot coding length distribution
    if(plot_coding_length_distri):
        plt.hist(coding_len,bins=1000)
        plt.title("Coding Length")
        plt.xlim(0,20000)        
        plt.show()
        
    # plot coding vs gene length
    if(plot_coding_vs_gene_length):
        plt.plot(gene_len,coding_len,linestyle="",marker=".")
        plt.xlim(0,20000)
        plt.ylim(0,20000)
        plt.xlabel("Gene Length")
        plt.ylabel("Coding Length")
        plt.show()
    
    # plot gc_content_distri
    if(plot_gc_content_distri):
        plt.hist(gc_per_gene,bins=200)
        plt.xlim(0.1,0.9)
        plt.ylim(0,1000)
        plt.title("GC content")
        plt.show()
     
    # plot gc vs coding length
    if(plot_gc_vs_coding_len):
        plt.plot(coding_len,gc_per_gene,linestyle="",marker=".")
        plt.xlabel("Gene Length")
        plt.ylabel("GC")
        plt.xlim(0,20000)
        plt.show()
        
    get_start_stop_codon(use_genome)
    
    # get structured genome data from a file
def get_genome_data(use_file):
    
    genome_file = open(use_file,"r")

    genome_genes = {}
    for line in genome_file:
        line = line.rstrip()
        if ( len(line) > 0):
            data = line.split("|")
            if (data[0][0] == ">"):
                gene_info = {"start":int(data[1]),"stop":int(data[2]),"contig":data[3],"sequence":""}
                use_key = data[0]
                genome_genes[use_key] = gene_info
            else:
                genome_genes[use_key]["sequence"] += data[0]

    genome_file.close()
    return genome_genes

    
