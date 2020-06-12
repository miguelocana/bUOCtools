from time import sleep
def translate(seq):        
    table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    } 
    protein ="" 
    if len(seq)%3 == 0: 
        for i in range(0, len(seq), 3): 
            cod = seq[i:i + 3] 
            protein += table[cod] 
    return protein 
def posiciones(codon):
    codon = codon.upper()
    nucleotidos = ['A','G','T','C']
    split_codon = [i for i in codon]

    # Primer check, tamaño (3)
    if len(codon)!=3:
        print('Introduce una secuencia múltiplo de 3 nucleótidos.')
    
    # Segundo check, los codones
    for i in split_codon:
        if i in nucleotidos:
            pass
        else:
            print('Introduce nucleótidos válidos (carácter erróneo: {0})'.format(i))
    
    # Vamos a modificar las posiciones de los nucleótidos
    else: 
        translated_codon = translate(codon)
        
        codones = []
        aas = []
        for i in range(0,len(split_codon)):
            split_codon = [i for i in codon]
            nucleotidos = ['A','G','T','C']
            nucleotidos.remove(split_codon[i])
            columna_codones = []
            columna_aas = []
            for j in range(0,len(codon)):
                split_codon[i] = nucleotidos[j]
                columna_codones.append(''.join(split_codon))
                columna_aas.append(translate(''.join(split_codon)))
            codones.append(columna_codones)
            aas.append(columna_aas)
    
    # Comparamos los aminoácidos con el aminoácido principal para obtener s y n
    codon_aa = translate(codon)
    sn = []
    for i in aas:
        s, n = 0, 0
        for j in i:
            if j == codon_aa:
                s += 1
            else:
                n += 1
        try:
            sn.append([s/len(codon),n/len(codon)])
        except ZeroDivisionError:
            return None
    
    # Sumamos todo
    s = 0
    n = 0
    for i in sn:
        s += i[0]
        n += i[1]
    
    n_posiciones = [s,n]
    
    for i in '.'*5:
        print(i)
        sleep(0.3)
        
    
    print('Codón: ', codon)
    print('Aminoácido: ', codon_aa)
    sleep(0.5)
    print('-'*20)
    print('Posiciones cambiadas por codón: ')
    print(codon[0],': ',codones[0])
    print(codon[1],': ',codones[1])
    print(codon[2],': ',codones[2])
    sleep(0.5)
    print('-'*20)
    print('Aminoácidos: ')
    print(aas[0])
    print(aas[1])
    print(aas[2])
    sleep(0.5)
    print('-'*20)
    print('Resultados por nucleótido: ')
    for i in range(0,len(codon)):
        print('s{0} = '.format(codon[i]),sn[i][0])
        print('n{0} = '.format(codon[i]),sn[i][1])
        print('-'*5)
    sleep(0.5)
    print('-'*20)
    print('Resultado codón')
    print('s{0} = '.format(codon),n_posiciones[0])
    print('n{0} = '.format(codon),n_posiciones[1])
    
    return n_posiciones

def posiciones_seq(seq):
    seq = seq.upper()
    nucleotidos = ['A','G','T','C']
    seq_lista = [i for i in seq]
    if len(seq)%3 != 0:
        print('Error, introduzca una secuencia múltiplo de 3.')
    else:
        count_codones = len(seq) / 3
        codones = []
        x_count = -3
        y_count = 0
        
        for i in range(0,int(count_codones)):
            codon = []
            x_count += 3
            y_count += 3
            for j in range(x_count,y_count):
                codon.append(seq_lista[j])
            codones.append(''.join(codon))
        
        for i in '.'*5:
            print(i)
            sleep(0.3)
        
        print('Codones detectados: ', len(codones))
        for i in codones:
            print(i)
        
        # Se calculan las posiciones sinónimas y no sinónimas
        resultado = []
        for i in codones:
            resultado.append(posiciones(i))
            
        s_total = 0
        n_total = 0
        for i in resultado:    
            s_total += i[0]
            n_total += i[1]
    
    for i in '.'*5:
        print(i)
        sleep(0.3)
        
    print('Resultado de la secuencia: ')
    print('s{0} = '.format(seq),s_total)
    print('n{0} = '.format(seq),n_total)
    sleep(0.5)
    return [s_total, n_total]

def posiciones_par_seq(seq1,seq2):
    seqs = [seq1,seq2]
    s_par_total = 0
    n_par_total = 0
    for i in seqs:
        s,n = posiciones_seq(i)
        s_par_total += s
        n_par_total += n  
    for i in '.'*5:
        print(i)
        sleep(0.3)
    print('Resultado de las posiciones del par de secuencias: ')
    print('s: ',s_par_total/2)
    print('n: ',n_par_total/2)
    return [s_par_total/2,n_par_total/2]

print('Bienvenido a la calculadora de número de posiciones sinónimas y no sinónimas de un par de secuencias.')
seq_input1 = input('Introduce la secuencia 1: ')
seq_input2 = input('Introduce la secuencia 2: ')
posiciones_par_seq(seq_input1,seq_input2)