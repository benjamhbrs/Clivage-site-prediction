import pandas as pd
import math

#FONCTIONS

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


def initiate_dict(p,q):
    f = {}
    g = {}
    for i in range(26):
        g[chr(ord('A')+i)] = 0

    for i in range(26):
        for j in range(-p,q):
            f[chr(ord('A')+i),j] = 0
    return f,g


#PREPORCESSING
# pas utile ici mais peut être pour autre jeu de donnéees

def pre_processing(data,p,q):
    for i in range(len(data)):         
        sequence = data.at[i,'sequence']
        clivage_place = data.at[i,'clivage_place']
        length = len(sequence)
        beginning_mature_prot_index = clivage_place.find('C')
        if (beginning_mature_prot_index - p < 0) or (beginning_mature_prot_index + q - 1 >= length):
            data.drop(i,inplace = True)
    return data.reset_index(drop = True)


def get_information(data,sequence_number):
    sequence = data.at[sequence_number,'sequence']
    clivage_place = data.at[sequence_number,'clivage_place']
    length = len(sequence)
    beginning_mature_prot_index = clivage_place.find('C')
    return length, beginning_mature_prot_index

def get_number_of_AA(train_data):
    number_of_AA = 0
    for i in range(len(train_data)):
        number_of_AA += get_information(train_data,i)[0]
    number_of_AA -= len(train_data) #on enlève les acides aminés méthionine qui sont toujours au début d'une protéine
    return number_of_AA

#TRAITEMENT STATISTIQUE DU TRAINING SET

def fill_in_freq(train_data, sequence_number, length_train_data, number_of_AA):
    sequence = train_data.at[sequence_number,'sequence']
    clivage_place = train_data.at[sequence_number,'clivage_place']
    length = len(sequence)
    beginning_mature_prot_index = clivage_place.find('C')
    for i in range(1,length):
        g[sequence[i]] += 1.0/number_of_AA
    if (beginning_mature_prot_index - p < 0) or (beginning_mature_prot_index + q - 1 >= length):
        print('abnormal_sequence_spotted')
    else:
        for j in range(beginning_mature_prot_index-p,beginning_mature_prot_index+q-1):
            f[(sequence[j],j-beginning_mature_prot_index)] += 1/length_train_data


def train(train_data):
    length = len(train_data)
    number_of_AA = get_number_of_AA(train_data)
    for k in range(length):
            fill_in_freq(train_data,k,length,number_of_AA)


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


def test(test_data):
    length = len(test_data)
    success_rate = 0
    for i in range(length):
        success_rate += predict_clivage(test_data,i)/length
    return success_rate


## UTLILISATION
if __name__=="__main__":

    p = 13
    q = 2 #on remarque que dans EUKSIG pour toute séquence, length-beginning_mature_prot_index vaut 30 donc on peut prendre q dans [1,30] 
            #en gardant toutes les protéines

    data_eucaryotes = load_data('EUKSIG_13.txt')
    data_procaryotes = load_data('GRAM+SIG_13.txt')
    f,g = initiate_dict(p,q)

    data_eucaryotes = pre_processing(data_eucaryotes,p,q) #on vire les données telles que le site ne soit pas entourée d'une séquence de longueur p+q
    data_procaryotes = pre_processing(data_procaryotes,p,q) #en fait ne fait rien sur les données fournies, et si jamais il faaut dégager une ligne penser à changer 
                                                #la fonction pour réindexer la dataframe, sinon erreur : FAIT
    #print(data_eucaryotes)
    train_data_euc = data_eucaryotes.iloc[0:905,:]
    test_data_euc = data_eucaryotes.iloc[905:1005,:].reset_index(drop = True)
    test_data_proc = data_procaryotes

    train(train_data_euc)
    print(f'success rate for training set euc. and testing set euc. is {test(test_data_euc)}')
    print(f'success rate for training set euc. and testing set proc. is {test(test_data_proc)}')
