'''

Phenotype prediction model using topological summaries 
generated from persistence images

'''

import dcor
import math
import joblib
import numpy as np
import pandas as pd
from tqdm import tqdm
from gudhi import RipsComplex
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.neural_network import MLPClassifier
from gudhi.representations import PersistenceImage
from sklearn.model_selection import train_test_split

class Preprocess:

    def __init__(self, patient_data):
        
        self.patient_data = patient_data

    def preprocess(self, ref):
        """
        Read pickle file and pre-select cancer signalling genes from Genes.csv file curated from KEGG
        :return: a preselected gene expression dataset 
        """
        self.df = pd.read_pickle(self.patient_data)
        self.genes = pd.read_csv("/Users/lebohangmashatola/downloads/DATA/Genes.csv", sep=';')
        self.genes = list(self.genes["Genes"])
        self.data = self.df[self.df.columns.intersection(self.genes)]        
        return self.data

class Prediction:

    def __init__(self, resistance, sensitive, signedTOM, output_dir_MLP):
        
        self.resistance = resistance
        self.sensitive = sensitive
        self.signedTOM = signedTOM
        self.output_dir_MLP = output_dir_MLP

    def split(self):
        """"
        Split gene expression dataset into train and test subsets
        :return: train expression & labels; test expression & labels
        """
        self.all_patients = pd.concat([self.cancer, self.healthy])
        self.labelencoding = preprocessing.LabelEncoder()
        self.labels = self.labelencoding.fit_transform(self.all_patients.index.to_list())

        self.train, self.test = train_test_split(self.all_patients, test_size=0.3, shuffle=True)
        self.train_labs = self.labelencoding.fit_transform(self.train.index.to_list())
        self.test_labs = self.labelencoding.fit_transform(self.test.index.to_list())

        self.train, self.test = np.array(self.train), np.array(self.test)

        return self.train, self.test

    def intergene_correlation_measure(self):
        """"
        Calculate the distance correalation measures across training patients
        :return: distance correlation matrix across training patients
        """
        self.num_genes = np.array(self.cancer.shape[1])
        self.dist = np.zeros((self.num_genes, self.num_genes))
        
        for i in tqdm(range(self.num_genes)):
            for j in range(i + 1, self.num_genes):
                self.dist[i, j] = dcor.distance_correlation(np.array(self.cancer)[:, i], np.array(self.cancer)[:, j])
                
        self.dist = self.dist + self.dist.T + np.eye(self.num_genes)
        self.dist = 1 - self.dist
        return self.dist

    def intergene_signed_TOM(self):
        """"
        Read a calculated signed-TOM across training patients
        :return: signed-TOM matrix across training patients
        """
        self.dist = pd.read_csv(self.signedTOM, index_col=0)
        return self.dist
    
    return np.array(df)

    def plot_corr(self):
        """"
        Plot signed-TOM/Distance correlation matrix
        :return: Visualised correlation matrix
        """
        
        f = plt.figure(figsize=(19, 15))
        
        plt.matshow(self.dist, fignum=f.number)
        plt.xticks(range(self.dist.shape[1]), fontsize=12, rotation=45)
        plt.yticks(range(self.dist.shape[1]), fontsize=12)
        
        cb = plt.colorbar()
        cb.ax.tick_params(labelsize=12)
        
        plt.title('Distance Correlation Matrix', fontsize=16)
        plt.show()


    def patient_correlation_measure(self, F):
        """
        Calculate the adjusted distance correlation measure for individual patients  
        :return: adjusted distance correlaion measure
        """
        F = F.T
        self.patient_dist = np.zeros((self.num_genes, self.num_genes))
        
        for i in range(self.num_genes):
            for j in range(i+1, self.num_genes):
                self.patient_dist[i,j] = self.dist[i,j] + math.sqrt(F[i] + F[j])/10 
        
        self.patient_dist = self.patient_dist + self.patient_dist.T + np.eye(self.num_genes)

        return self.patient_dist
    

    def filter_persistence(self, persistence_intervals, top_percentile):]
        """
        Filters based on the lifespan and select the top percentile persisting features
        :return: filtered persistence betti-numbers based on lifespan
        """
        dict_birth_death = dict(enumerate(persistence_intervals))
        dict_lifespan = dict(enumerate([(value[1] - value[0]) for index, value in dict_birth_death.items()])) 
        dict_sorted = dict(sorted(dict_lifespan.items(), key=lambda x: x[1])) 

        subset = math.floor(len(persistence_intervals) * int(top_percentile) / 100) 
        subset_dict = {k: dict_sorted[k] for k in list(dict_sorted)[:subset]}
        keys = [i for i in subset_dict.keys()]

        persistence_pairs = [persistence_intervals[i] for i in keys] 
        
        return (np.array(persistence_pairs))


    def PersistentHomology(self, patients, top_percentile=50): #default take the 50% most persistent Betti numbers
        """Computes the persistence of the simplicial complex. Saves as persistent betti pairs
        """
        Persistent_diagrams1, Persistent_diagrams2 = [], [] #only consider > Betti-0
        
        for patient in tqdm(range(len(patients))):
            
            self.patient_exprs = patients[patient]
            self.distance_matrix = self.patient_correlation_measure(self.patient_exprs) #Weights used include per-patient gene expressions
            self.rips_complex = RipsComplex(self.distance_matrix, max_edge_length = float('Inf'), sparse=0.5).create_simplex_tree(max_dimension=1) 
        
            self.rips_complex.collapse_edges()
            self.rips_complex.expansion(3)
            self.rips_complex.persistence()

            self.Persistent_diagrams1.append(self.filter_persistence(self.rips_complex.persistence_intervals_in_dimension(1), top_percentile))
            self.Persistent_diagrams2.append(self.filter_persistence(self.rips_complex.persistence_intervals_in_dimension(2), top_percentile))

        remove_infinity = lambda barcode : np.array([bars for bars in barcode if bars[1]!= np.inf]) #remove infinity Betti-numbers 
        
        self.Persistent_diagrams1 = list(map(remove_infinity, self.Persistent_diagrams1))
        self.Persistent_diagrams2 = list(map(remove_infinity, self.Persistent_diagrams2))
        
        self.image = PersistenceImage(resolution=[20, 20])
        self.samplelandscape1img = self.image.fit_transform(self.Persistent_diagrams1)
        self.samplelandscape2img = self.image.fit_transform(self.Persistent_diagrams2)
    
        return np.column_stack((self.samplelandscape1img, self.samplelandscape2img))


    def MLP_Model(self, image_train, image_test):
         """Training a MLP neural network classifier
            :return:saved MLP neural weights 
        """
        bst = MLPClassifier(random_state=1, activation='logistic', learning_rate=0.001, max_iter=300, validation_fraction=0.1)
        bst.fit(image_train, self.train_labs)
        joblib.dump(bst, self.output_dir_XGBoost) #save model

    def test_MLP_Model(self, image_test, patient_number):
        """Testing a MLP neural network classifier
            :return:probability score for a particular patient
        """
        mlp_model_latest = MLPClassifier() 
        mlp_model_latest = joblib.load(self.output_dir_MLP)
        
        prediction = mlp_model_latest.predict(image_test)
        prediction_proba = mlp_model_latest.predict_proba(image_test)
        pred = prediction_proba[patient_number]
        
        if prediction[patient_number] == 0:
            print("Label: FOLFOX Resistance")
            print("Confidence Score: " + str(round(pred[0],2) * 100) + "%")
            
        elif prediction[patient_number] == 1:
            print("Label: FOLFOX Sensitive")
            print("Confidence Score: " + str(round(pred[1],2) * 100) + "%")


if __name__ == '__main__':

    resistance = Preprocess("/dir/Resistance.pkl")
    resistance_patients = resistance.preprocess("R") 

    sensitive = Preprocess("/dir/documents/Sensitive.pkl")
    sensitive_patients = sensitive.preprocess("S") 

    prediction = Prediction(cancer_patients, normal_patients, "dir/MLP_Model.json")
    X_train, X_test = prediction.split() 
    
    prediction.intergene_correlation_measure() 
    prediction.intergenes_signed_TOM() 
    
    prediction.plot_corr() 
    
    Train_samples = prediction.PersistentHomology(X_train, 50) #Change Percentile persisting features
    Test_samples = prediction.PersistentHomology(X_test, 50) #Change Percentile persisting features

    prediction.MLP_Model(Train_samples, Test_samples) 
    prediction.test_MLP_Model(Test_samples, 10) #Select patient number you want to evaluate