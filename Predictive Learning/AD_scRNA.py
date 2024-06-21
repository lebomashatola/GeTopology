
import time
import numpy as np
import pandas as pd
import scanpy as sc
from tqdm import tqdm
import anndata as ad
import seaborn as sns
from sklearn.svm import SVC
import matplotlib.pyplot as plt
from imblearn.pipeline import Pipeline
from qiskit.primitives import Sampler
from qiskit_algorithms.optimizers import COBYLA
from sklearn.metrics import classification_report
from sklearn.metrics import f1_score, make_scorer
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from imblearn.under_sampling import RandomUnderSampler
from qiskit_machine_learning.algorithms.classifiers import VQC
from qiskit.circuit.library import ZZFeatureMap, RealAmplitudes
from sklearn.model_selection import ShuffleSplit, cross_val_score
from qiskit_machine_learning.circuit.library import RawFeatureVector


class VQCWrapper:
    def __init__(self, feature_dimension: int, reps=3, maxiter=10):
        self.feature_dimension = feature_dimension
        self.reps = reps
        self.maxiter = maxiter
        self.vqc = self.prepare_vqc()

    def prepare_vqc(self):
        """
        Iniialise the parameters for the VQC
        """
        sampler = Sampler()
        feature_map = RawFeatureVector(feature_dimension=self.feature_dimension)
        ansatz = RealAmplitudes(feature_map.num_qubits, reps=self.reps)
        
        vqc = VQC(
            sampler=sampler,
            feature_map=feature_map,
            ansatz=ansatz,
            optimizer=COBYLA(maxiter=150),
            callback = None
        )
        return vqc

    def fit(self, X, y):
        self.vqc.fit(X, y)
        return self

    def predict(self, X):
        return self.vqc.predict(X)

    def get_params(self, deep=True):
        return {
            "feature_dimension": self.feature_dimension,
            "reps": self.reps,
            "maxiter": self.maxiter,
        }

    def set_params(self, **params):
        for key, value in params.items():
            setattr(self, key, value)
        self.vqc = self.prepare_vqc()  # Re-prepare VQC with updated parameters
        return self
    

class scGNN:
    def __init__(self, file_in : str, cell_labels : str, num_genes: int, num_cells: int):
        """
        Initialise the number of genes and cells
        """
        self.file_in = file_in
        self.cell_labels = cell_labels
        self.num_genes = num_genes
        self.num_cells = num_cells

    def define_data(self):
        
        # Read data and labels 
        counts = pd.read_csv('./GSE138852_counts.csv', index_col=0)
        cell_types = pd.read_csv('./GSE138852_covariates.csv', index_col=0)
        counts = counts.transpose()
        
        # Create a new column
        counts['keys'] = counts.index
        cell_types['keys'] = cell_types.index
        
        # Create key to match UID and cell type
        labels = counts.merge(cell_types, on='keys', how='inner')[['keys', 'oupSample.cellType']]
        counts = counts.drop(['keys'],axis =1)

        select_labels = ['astro', 'OPC', 'neuron']
        filtered_labels = labels[labels['oupSample.cellType'].isin(select_labels)]
        filter_cells = filtered_labels['keys'].to_list()

        mask = counts.index.isin(filter_cells)
        filtered_df_counts = counts[mask]

        # Create AnnData Object
        adata_ah = ad.AnnData(X=filtered_df_counts.values.T) # Rows correspond to cells and columns to genes
        adata_ah.var_names = filtered_df_counts.index
        adata_ah.obs_names = filtered_df_counts.columns
        
        # Filter and normalise scRNA-Seq selecting min gene counts
        sc.pp.filter_genes(adata_ah, min_counts=1) 
        sc.pp.filter_cells(adata_ah, min_counts=1)
        sc.pp.normalize_total(adata_ah, target_sum=1)
        sc.pp.scale(adata_ah)

        adata_ah = adata_ah.transpose()

        # Obtain highly variable genes
        sc.pp.highly_variable_genes(adata_ah, n_top_genes=64, flavor='cell_ranger')
        adata = adata_ah[:, adata_ah.var.highly_variable]
        selected_indices = np.random.choice(len(adata), 500, replace=False)
        filtered_adata_output = adata.X[selected_indices, :]
        selected_labels = filtered_labels.iloc[selected_indices]

        # Print the distribution of labels after selection
        label_counts = selected_labels['oupSample.cellType'].value_counts()
        print("Distribution of cell types after selection:")
        print(label_counts)
        
        self.y = selected_labels['oupSample.cellType']
        self.X = np.array(filtered_adata_output)

        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(self.X, self.y, test_size=0.33, random_state=0)
        
        return self.X

    def perform_one_vs_all_classification(self):

        self.X = np.array(self.X)
        unique_classes = np.unique(self.y)

        self.results_df = pd.DataFrame(
            columns=[
                "Class",
                "Classifier",
                "Sampling_Technique",
                "F1-Score",
                "Mean F1-Score CV",
                "Time",
            ]
        )

        classifiers = {
            "RandomForest "
            + str(self.num_genes): RandomForestClassifier(
                n_estimators=100, random_state=42
            ),
            "SVM " + str(self.num_genes): SVC(kernel="linear", random_state=42),
            "VQC " + str(self.num_genes): VQCWrapper(feature_dimension=self.X.shape[1]),
        }

        # Define the F1 score as the evaluation metric
        f1_scorer = make_scorer(f1_score, average="weighted")

        for class1 in tqdm(unique_classes, desc="One-vs-All Classifications"):
            y_binary = np.where(self.y == class1, 1, 0)
            self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(
                self.X, y_binary, test_size=0.33, random_state=0
            )

            pipeline = Pipeline(
                steps=[("resample", RandomUnderSampler(sampling_strategy="majority"))]
            )

            X_train_res, y_train_res = pipeline.fit_resample(self.X_train, self.y_train)

            for clf_name, clf in classifiers.items():

                start_time = time.time()
                # Perform cross-validation on the resampled training set and calculate the mean F1 score
                scores = cross_val_score(
                    clf,
                    X_train_res,
                    y_train_res,
                    cv=ShuffleSplit(n_splits=5, test_size=0.33, random_state=42),
                    scoring=f1_scorer,
                )
                mean_f1_score_cv = np.mean(scores)

                clf.fit(X_train_res, y_train_res)
                y_pred = clf.predict(self.X_test)
                end_time = time.time()

                time_taken = end_time - start_time
                report = classification_report(self.y_test, y_pred, output_dict=True)
                f1_score_test = report["weighted avg"]["f1-score"]

                new_row_df = pd.DataFrame(
                    [
                        {
                            "Class": class1,
                            "Classifier": clf_name,
                            "Sampling_Technique": "RandomUnderSampler",
                            "F1-Score": f1_score_test,
                            "Mean F1-Score CV": mean_f1_score_cv,
                            "Time": time_taken,
                        }
                    ]
                )

                self.results_df = pd.concat(
                    [self.results_df, new_row_df], ignore_index=True
                )

        self.results_df = self.results_df.sort_values(by=["Classifier", "Class"])
        self.results_df.to_csv(
            "./F1-Scores-" + str(self.num_genes) + "_" + str(self.num_cells), index=None, header=True
        )

        return self.results_df

    def plot_f1_matrix(self):

        # Pivot the results to create a matrix for F1-Scores
        f1_pivot = self.results_df.pivot_table(
            index="Class", columns=["Classifier"], values="F1-Score"
        )

        # Plot the heatmap of F1-Scores
        plt.figure(figsize=(14, 7))

        # Define the heatmap
        sns.heatmap(f1_pivot, annot=False, cmap="viridis", fmt=".2f")
        plt.title("F1-Score Matrix for One-vs-All Classification")
        plt.ylabel("Class")
        plt.xlabel("Classifier and Sampling Technique")
        plt.xticks(rotation=35, ha="right")

        # Modify figure
        plt.rcParams["axes.grid"] = False
        plt.rcParams["figure.dpi"] = 150

        # Save figure
        plt.savefig("./heatmap_" + str(self.num_genes) + "_" + str(self.num_cells) + ".png", bbox_inches="tight")
        #plt.show()
        
        
if __name__ == "__main__":

    # 8 Dimensions
    experiment1 = scGNN("./GSE138852_counts.csv", "./GSE138852_covariates.csv", 8, 500)
    experiment1.define_data()
    experiment1.perform_one_vs_all_classification()
    experiment1.plot_f1_matrix()
    
    # 16 Dimensions
    experiment1 = scGNN("./GSE138852_counts.csv", "./GSE138852_covariates.csv", 16, 500)
    experiment1.define_data()
    experiment1.perform_one_vs_all_classification()
    experiment1.plot_f1_matrix()
    
    # 32 Dimensions
    experiment1 = scGNN("./GSE138852_counts.csv", "./GSE138852_covariates.csv", 32, 500)
    experiment1.define_data()
    experiment1.perform_one_vs_all_classification()
    experiment1.plot_f1_matrix()
    
    # 64 Dimensions
    experiment1 = scGNN("./GSE138852_counts.csv", "./GSE138852_covariates.csv", 64, 500)
    experiment1.define_data()
    experiment1.perform_one_vs_all_classification()
    experiment1.plot_f1_matrix()
    
    # 128 Dimensions
    experiment1 = scGNN("./GSE138852_counts.csv", "./GSE138852_covariates.csv", 128, 500)
    experiment1.define_data()
    experiment1.perform_one_vs_all_classification()
    experiment1.plot_f1_matrix()
    
    # 256 Dimensions
    experiment1 = scGNN("./GSE138852_counts.csv", "./GSE138852_covariates.csv", 256, 500)
    experiment1.define_data()
    experiment1.perform_one_vs_all_classification()
    experiment1.plot_f1_matrix()
    
    # 512 Dimensions
    experiment1 = scGNN("./GSE138852_counts.csv", "./GSE138852_covariates.csv", 512, 1000)
    experiment1.define_data()
    experiment1.perform_one_vs_all_classification()
    experiment1.plot_f1_matrix()