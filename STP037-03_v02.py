#coding:utf-8


import csv
import re

def IL_replace(list):

    x = ''
    list_replaced = []

    for i in list:
        if 'I' in i:
            x = i.replace('I', 'L')
            list_replaced.append(x)
        else:
            list_replaced.append(i)

    return list_replaced



def IL_replace2(str):

    if 'I' in str:
        str_replaced = str.replace('I', 'L')
            
    else:
        str_replaced = str

    return str_replaced
    
    
def interval_calculate(fragment_N, fragment_C, replaced_proteinSeq):
    
##N側fragment→C側fragmentの順でprotein内に存在した時の最短距離を計算
    
    #左から検索、最初にヒットした位置を返す
    fragment_N_end=[me.end() for me in re.finditer(fragment_N, replaced_proteinSeq)] #複数ヒットした場合にすべてのfragment終了位置をリストで取得
    fragment_C_start=[os.start() for os in re.finditer(fragment_C, replaced_proteinSeq)] #複数ヒットした場合にすべてのfragment開始位置をリストで取得
    
    #距離候補計算
    interval_Candidate=[]
    for e in fragment_N_end:
        for s in fragment_C_start:
            interval_Candidate.append(s-e)
    #正の数字のみ
    interval_Candidate_plus=[]     
    for p in interval_Candidate:
        if p > 0:
            interval_Candidate_plus.append(p)
    
    if len(interval_Candidate_plus) > 0: #正の数字（Nfragment→Cfragmentの順でprotein内に存在）がひとつでもあれば
        #正順での最小値
        interval_plus = min(interval_Candidate_plus)
    else: 
        interval_plus=0

##C側fragment→N側fragmentの順でprotein内に存在した時の最短距離を計算
    fragment_N_start=[ms.start() for ms in re.finditer(fragment_N, replaced_proteinSeq)] 
    fragment_C_end=[oe.end() for oe in re.finditer(fragment_C, replaced_proteinSeq)] 
    
    #距離候補計算
    interval_Candidate2=[]
    for s2 in fragment_N_start:
        for e2 in fragment_C_end:
            interval_Candidate2.append(s2-e2)
    #正の数字のみ
    interval_Candidate_plus2=[]     
    for p2 in interval_Candidate2:
        if p2 > 0:
            interval_Candidate_plus2.append(p2)
    
    if len(interval_Candidate_plus2) > 0: #正の数字（Cfragment→Nfragmentの順でprotein内に存在）がひとつでもあれば
        #正順での最小値
        interval_plus2 = min(interval_Candidate_plus2)
    else: 
        interval_plus2=0
    
##2つの最小値が得られた場合はより小さい方を採用
    if interval_plus==0: #もしN→Cには候補距離なければC→Nの最小距離で決定
        if interval_plus2==0: #C→Nも候補距離がない(=どちらも負の数字しかなかった）のは、fragmentのタンパク内位置が重複しているため。
            interval='NA'
            splice_type='NA'
        else:
            interval=interval_plus2 #C→Nに最小距離候補があれば決定
            splice_type='Reverse'
    else: #N→Cに候補距離あって
        if interval_plus2==0: #C→Nには候補距離なければN→Cの最小距離で決定
            interval=interval_plus
            splice_type='Standard'
        else: #どちらにも候補距離があれば、より小さい方で決定
            if interval_plus <= interval_plus2:
                interval=interval_plus
                splice_type='Standard'
            else:
                interval=interval_plus2
                splice_type='Reverse'
    
    return interval,splice_type

def gene_get(annotation): 
    if annotation[0:4] == '>sp|': #swissの場合
        
        if 'GN=' in annotation:
            if 'PE=' in annotation:
                gene=annotation.split(' PE=')[0].split('GN=')[1]
            else:
                gene=annotation.split('GN=')[1].split('\n')[0]
        else:
            gene='NA'
            


    return gene

def expression_get(gene):
    f=open(samplename+'/'+'expression.csv')
    EX = f.readlines()
    f.close
    
    flag = 0
    for t in range(1, len(EX)):
        if EX[t].split(',')[0] == gene:
            flag = 1
            expression=EX[t].split(',')[1][:-1]
    
    if flag == 0:
        expression='NA'
            
    return expression

#main======================================================

print ('sample name=')
samplename = raw_input()

print ('Reference=Swiss/Other?')
ref = raw_input()

if ref == 'Swiss':
    print ('Taxonomy=Human/Mouse?')
    taxonomy = raw_input()
    
    if taxonomy == 'Human':
        f = open ('../../CommonDatabases/patho1_SwissProt_human_20180314.delLF.fasta')
        FASTA = f.readlines()
        f.close
        
    elif taxonomy == 'Mouse':
        f = open ('../../CommonDatabases/patho1_SwissProt_mouse_20190510_delLF.fasta')
        FASTA = f.readlines()
        f.close

elif ref == 'Other': 
    print ('Reference FASTA = XX.fasta')
    reference = raw_input()
       
    f = open (samplename+'/'+ reference + '.fasta')
    FASTA = f.readlines()
    f.close

print ('Expression=TPM/FPKM/RPKM/ND?')
exp = raw_input()


f = open (samplename+'/'+samplename + '.csv')
peptides = f.readlines()
f.close

peptide = []
length = []
Modification =[]
Affinity = []
Scan =[]
ALC =[]
ALC_rank = []
mz = []
z_ = []
RT = []
Mass = []
ppm = []
localconfidence = []


for i in range (1, len(peptides)):
    peptide.append(peptides[i].split(',')[0])
    length.append(int(peptides[i].split(',')[1]))
    Modification.append(peptides[i].split(',')[2])
    Affinity.append(peptides[i].split(',')[3])
    Scan.append(peptides[i].split(',')[5])
    ALC.append(peptides[i].split(',')[6])
    ALC_rank.append(peptides[i].split(',')[7])
    mz.append(peptides[i].split(',')[8])
    z_.append(peptides[i].split(',')[9])
    RT.append(peptides[i].split(',')[10])
    Mass.append(peptides[i].split(',')[11])
    ppm.append(peptides[i].split(',')[12])
    localconfidence.append(peptides[i].split(',')[13][:-1])
    
        
#peptideのIL置換
peptide_replaced = IL_replace(peptide)


#FASTAのIL置換
FASTA_replaced = []
for f in range (0, len(FASTA)-1, 2):
    FASTA_replaced.append(FASTA[f][:-1])
    
    FASTA_seq_replaced = IL_replace2(FASTA[f+1][:-1])
    FASTA_replaced.append(FASTA_seq_replaced)

#print FASTA_replaced

f = open (samplename+'/'+samplename + '_splicedv02.csv', 'a')
f.write(\
'Peptide' + ','+\
'Peptide_length' + ','+\
'Fragment1' + ','+\
'Fragment2' + ','+\
'Fragment1_length' + ','+\
'Fragment2_length' + ','+\
'Fragment距離(AA)' + ','+\
'Splice_type' + ','+\
'Accession' + ','+\
'Gene' + ','+\
'Expression(' + exp + ')' + ','+\
'Modification' + ','+\
'%Rank(LI注意)' + ','+\
'Scan' + ','+\
'ALC(%)' + ','+\
'ALC_rank' + ','+\
'm/z' + ','+\
'z' + ','+\
'RT' + ','+\
'Mass' + ','+\
'ppm' + ','+\
'localconfidence' + ','+\
'\n')
f.close


for m in range (0, len(peptide_replaced)):

    print peptide_replaced[m]
    
    #N末側からのみfragmentを作成して検索を行う  
    for z in range (2, length[m]-1):
        #N末から2fragment以上(ABCDEFGHI:9merだったらAB～ABCDEFG)
        fragment1_Nver = peptide_replaced[m][0:z]
        #print fragment1_Nver
        
        for n in range (1, len(FASTA_replaced), 2):
            #print FASTA_replaced[n]
            
            if fragment1_Nver in FASTA_replaced[n]:
                       
                #同じエントリ内にもう片方の断片が含まれているか
                fragment2_Nver = peptide_replaced[m][z:]
                if fragment2_Nver in FASTA_replaced[n]:
                
                    #print (FASTA_replaced[n-1].split('|')[0] + '_' + FASTA_replaced[n-1].split('|')[1])
                    interval,splice_type=interval_calculate(fragment1_Nver, fragment2_Nver, FASTA_replaced[n])
                    gene=gene_get(FASTA_replaced[n-1])
                    
                    if exp == 'ND':
                        expression='ND' 
                    else:
                        expression=expression_get(gene)
                    
                    if interval != 'NA': #NAなのは2つのfragmentのタンパク内位置が重複している時
                    
                        f = open (samplename+'/'+samplename + '_splicedv02.csv', 'a')
                        #peptide, length書き出し
                        f.write(peptide[m]+',' + str(length[m]) + ',')
			#fragment情報書き出し
                        f.write(fragment1_Nver+',')
                        f.write(fragment2_Nver+',')
                        f.write(str(len(fragment1_Nver)) + ',' + str(len(fragment2_Nver)) + ',')
                        f.write(str(interval) + ',')
                        f.write(splice_type+',')
                        #accession書き出し
                        f.write(FASTA_replaced[n-1].split('|')[0] + '_' + FASTA_replaced[n-1].split('|')[1] + ',')                        
                        #Gene情報書き出し
                        f.write(gene + ',' + str(expression)+',')
                        #その他denovo情報書き出し
                        f.write(\
                        Modification[m] + ','+\
                        Affinity[m] + ','+\
                        Scan[m] + ','+\
                        ALC[m] + ','+\
                        ALC_rank[m] + ','+\
                        mz[m] + ','+\
                        z_[m] + ','+\
                        RT[m] + ','+\
                        Mass[m] + ','+\
                        ppm[m] + ','+\
                        localconfidence[m] + ','+\
                        '\n')
                        f.close
                        

                    
        