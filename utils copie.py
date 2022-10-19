#from webbrowser import get
import pandas as pd
import math

def load_data(file_name):

    f = open(file_name,'r')
    lines = f.readlines()
    useful_data = []
    count = 0
    for line in lines:
        if(count%3 == 1):
            useful_data.append([line.strip()])
        if(count%3 == 2):
            useful_data[(count-2)//3].append(line.strip())
        count+=1

    return pd.DataFrame(useful_data, columns = ['sequence','clivage_place'])

data = load_data('EUKSIG_13.txt')
data_procaryotes = load_data('GRAM+SIG_13.txt')


p = 13
q = 2

g = {}
f = {}

for i in range(26):
    g[chr(ord('A')+i)] = 0

for i in range(26):
    for j in range(-p,q):
        f[chr(ord('A')+i),j] = 0

def lign_to_be_deleted():
    return

#PREPORCESSING
# pas utile ici mais peut être pour autre jeu de donnéees

def pre_processing(data):
    for i in range(len(data)):         
        sequence = data.at[i,'sequence']
        clivage_place = data.at[i,'clivage_place']
        length = len(sequence)
        beginning_mature_prot_index = clivage_place.find('C')
        if beginning_mature_prot_index - p < 0 or beginning_mature_prot_index + q - 1 >= length:
            data.drop(i,inplace = True)

pre_processing(data)
data = data.reset_index(drop = True)
N_train = len(data)-100
print(N_train)
N_test = 100

def get_information(data,sequence_number):
    sequence = data.at[sequence_number,'sequence']
    clivage_place = data.at[sequence_number,'clivage_place']
    length = len(sequence)
    beginning_mature_prot_index = clivage_place.find('C')
    return length, beginning_mature_prot_index

number_of_AA = 0
for i in range(N_train):
    number_of_AA += get_information(data,i)[0]
number_of_AA -= N_train #on enlève les acides aminés méthionine qui sont toujours au début d'une protéine

#TRAITEMENT STATISTIQUE DU TRAINING SET

def fill_in_freq(data, sequence_number):
    sequence = data.at[sequence_number,'sequence']
    clivage_place = data.at[sequence_number,'clivage_place']
    length = len(sequence)
    beginning_mature_prot_index = clivage_place.find('C')

    for i in range(1,length):
        g[sequence[i]] += 1.0/number_of_AA

    if beginning_mature_prot_index - p < 0 or beginning_mature_prot_index + q - 1 >= length:
        print('spotted')
        return
    else:
        for j in range(beginning_mature_prot_index-p,beginning_mature_prot_index+q-1):
            f[(sequence[j],j-beginning_mature_prot_index)] += 1/N_train

for k in range(N_train):
    fill_in_freq(data, k)

#print(sum(g.values())) # cela doit vaaloir 1

#TESTING SET

def score(sequence):
    score = 0
    for j in range(-p,q):
        score += math.log(f[(sequence[p+j],j)]+0.0000000001) - math.log(g[sequence[p+j]]+0.000000001)
    return score

def predict_clivage(data, sequence_number):

    sequence = data.at[sequence_number,'sequence']
    clivage_place = data.at[sequence_number,'clivage_place']
    length = len(sequence)
    beginning_mature_prot_index = clivage_place.find('C')
    predicted_beginning_mature_prot_index = p
    best_score = score(sequence[0:p+q])

    for potential_beginning_mature_prot_index in range(p,length-q+1):
        temp_score = score(sequence[potential_beginning_mature_prot_index-p:potential_beginning_mature_prot_index+q])
        if temp_score>best_score:
            predicted_beginning_mature_prot_index = potential_beginning_mature_prot_index
            best_score = temp_score

    return predicted_beginning_mature_prot_index == beginning_mature_prot_index


success_rate = 0
for i in range(len(data)-100, len(data)):
    success_rate += predict_clivage(data,i)/N_test

print('success_rate_eucaryotes = ',success_rate, '(on a entrainé sur un training set procaryotes)')

# SUR PROCARYOTES

def predict_clivage(data,_sequence,_clivage_place):
    sequence = _sequence
    clivage_place = _clivage_place
    length = len(sequence)
    beginning_mature_prot_index = clivage_place.find('C')
    predicted_beginning_mature_prot_index = p
    best_score = score(sequence[0:p+q])

    for potential_beginning_mature_prot_index in range(p,length-q+1):
        temp_score = score(sequence[potential_beginning_mature_prot_index-p:potential_beginning_mature_prot_index+q])
        if temp_score>best_score:
            predicted_beginning_mature_prot_index = potential_beginning_mature_prot_index
            best_score = temp_score

    return predicted_beginning_mature_prot_index == beginning_mature_prot_index


success_rate_procaryote = 0
for i in range(len(data_procaryotes)):
    sequence = data_procaryotes.at[i,'sequence']
    clivage_place = data_procaryotes.at[i,'clivage_place']
    success_rate_procaryote += predict_clivage(data_procaryotes, sequence,clivage_place)/len(data_procaryotes)

print('success_rate_procaryote = ', success_rate_procaryote, '(les clivages procaryotes ne suivent pas les mêmes lois mais bien que le hasard : ancêtre commun ?)')

