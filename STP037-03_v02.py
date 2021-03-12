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
    
##N��fragment��C��fragment�̏���protein���ɑ��݂������̍ŒZ�������v�Z
    
    #�����猟���A�ŏ��Ƀq�b�g�����ʒu��Ԃ�
    fragment_N_end=[me.end() for me in re.finditer(fragment_N, replaced_proteinSeq)] #�����q�b�g�����ꍇ�ɂ��ׂĂ�fragment�I���ʒu�����X�g�Ŏ擾
    fragment_C_start=[os.start() for os in re.finditer(fragment_C, replaced_proteinSeq)] #�����q�b�g�����ꍇ�ɂ��ׂĂ�fragment�J�n�ʒu�����X�g�Ŏ擾
    
    #�������v�Z
    interval_Candidate=[]
    for e in fragment_N_end:
        for s in fragment_C_start:
            interval_Candidate.append(s-e)
    #���̐����̂�
    interval_Candidate_plus=[]     
    for p in interval_Candidate:
        if p > 0:
            interval_Candidate_plus.append(p)
    
    if len(interval_Candidate_plus) > 0: #���̐����iNfragment��Cfragment�̏���protein���ɑ��݁j���ЂƂł������
        #�����ł̍ŏ��l
        interval_plus = min(interval_Candidate_plus)
    else: 
        interval_plus=0

##C��fragment��N��fragment�̏���protein���ɑ��݂������̍ŒZ�������v�Z
    fragment_N_start=[ms.start() for ms in re.finditer(fragment_N, replaced_proteinSeq)] 
    fragment_C_end=[oe.end() for oe in re.finditer(fragment_C, replaced_proteinSeq)] 
    
    #�������v�Z
    interval_Candidate2=[]
    for s2 in fragment_N_start:
        for e2 in fragment_C_end:
            interval_Candidate2.append(s2-e2)
    #���̐����̂�
    interval_Candidate_plus2=[]     
    for p2 in interval_Candidate2:
        if p2 > 0:
            interval_Candidate_plus2.append(p2)
    
    if len(interval_Candidate_plus2) > 0: #���̐����iCfragment��Nfragment�̏���protein���ɑ��݁j���ЂƂł������
        #�����ł̍ŏ��l
        interval_plus2 = min(interval_Candidate_plus2)
    else: 
        interval_plus2=0
    
##2�̍ŏ��l������ꂽ�ꍇ�͂�菬���������̗p
    if interval_plus==0: #����N��C�ɂ͌�⋗���Ȃ����C��N�̍ŏ������Ō���
        if interval_plus2==0: #C��N����⋗�����Ȃ�(=�ǂ�������̐��������Ȃ������j�̂́Afragment�̃^���p�N���ʒu���d�����Ă��邽�߁B
            interval='NA'
            splice_type='NA'
        else:
            interval=interval_plus2 #C��N�ɍŏ�������₪����Ό���
            splice_type='Reverse'
    else: #N��C�Ɍ�⋗��������
        if interval_plus2==0: #C��N�ɂ͌�⋗���Ȃ����N��C�̍ŏ������Ō���
            interval=interval_plus
            splice_type='Standard'
        else: #�ǂ���ɂ���⋗��������΁A��菬�������Ō���
            if interval_plus <= interval_plus2:
                interval=interval_plus
                splice_type='Standard'
            else:
                interval=interval_plus2
                splice_type='Reverse'
    
    return interval,splice_type

def gene_get(annotation): 
    if annotation[0:4] == '>sp|': #swiss�̏ꍇ
        
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
    
        
#peptide��IL�u��
peptide_replaced = IL_replace(peptide)


#FASTA��IL�u��
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
'Fragment����(AA)' + ','+\
'Splice_type' + ','+\
'Accession' + ','+\
'Gene' + ','+\
'Expression(' + exp + ')' + ','+\
'Modification' + ','+\
'%Rank(LI����)' + ','+\
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
    
    #N��������̂�fragment���쐬���Č������s��  
    for z in range (2, length[m]-1):
        #N������2fragment�ȏ�(ABCDEFGHI:9mer��������AB�`ABCDEFG)
        fragment1_Nver = peptide_replaced[m][0:z]
        #print fragment1_Nver
        
        for n in range (1, len(FASTA_replaced), 2):
            #print FASTA_replaced[n]
            
            if fragment1_Nver in FASTA_replaced[n]:
                       
                #�����G���g�����ɂ����Е��̒f�Ђ��܂܂�Ă��邩
                fragment2_Nver = peptide_replaced[m][z:]
                if fragment2_Nver in FASTA_replaced[n]:
                
                    #print (FASTA_replaced[n-1].split('|')[0] + '_' + FASTA_replaced[n-1].split('|')[1])
                    interval,splice_type=interval_calculate(fragment1_Nver, fragment2_Nver, FASTA_replaced[n])
                    gene=gene_get(FASTA_replaced[n-1])
                    
                    if exp == 'ND':
                        expression='ND' 
                    else:
                        expression=expression_get(gene)
                    
                    if interval != 'NA': #NA�Ȃ̂�2��fragment�̃^���p�N���ʒu���d�����Ă��鎞
                    
                        f = open (samplename+'/'+samplename + '_splicedv02.csv', 'a')
                        #peptide, length�����o��
                        f.write(peptide[m]+',' + str(length[m]) + ',')
			#fragment��񏑂��o��
                        f.write(fragment1_Nver+',')
                        f.write(fragment2_Nver+',')
                        f.write(str(len(fragment1_Nver)) + ',' + str(len(fragment2_Nver)) + ',')
                        f.write(str(interval) + ',')
                        f.write(splice_type+',')
                        #accession�����o��
                        f.write(FASTA_replaced[n-1].split('|')[0] + '_' + FASTA_replaced[n-1].split('|')[1] + ',')                        
                        #Gene��񏑂��o��
                        f.write(gene + ',' + str(expression)+',')
                        #���̑�denovo��񏑂��o��
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
                        

                    
        