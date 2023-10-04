'''
Generate Topological Summaries for Each Patient
Represented by persistence diagrams 

'''

import os
import dcor
import math
import errno
import warnings
import cv2 as cv
import numpy as np
from os import path
import pandas as pd
from tqdm import tqdm
from gudhi import RipsComplex
import matplotlib.pyplot as plt
from functools import lru_cache
from sklearn import preprocessing
from sklearn.model_selection import train_test_split

from gudhi.reader_utils import read_persistence_intervals_in_dimension
from gudhi.reader_utils import read_persistence_intervals_grouped_by_dimension

__author__ = "Vincent Rouvreau, Bertrand Michel, Theo Lacombe"
__copyright__ = "Copyright (C) 2016 Inria"
__license__ = "MIT"

_gudhi_matplotlib_use_tex = True

def _array_handler(a):
    """
    :param a: if array, assumes it is a (n x 2) np.array and returns a
                persistence-compatible list (padding with 0), so that the
                plot can be performed seamlessly.
    :returns: * List[dimension, [birth, death]] Persistence, compatible with plot functions, list.
              * boolean Modification status (True if output is different from input)
    """
    if isinstance(a[0][1], (np.floating, float)):
        return [[0, x] for x in a], True
    else:
        return a, False


def _limit_to_max_intervals(persistence, max_intervals, key):
    """This function returns truncated persistence if length is bigger than max_intervals.
    :param persistence: Persistence intervals values list. Can be grouped by dimension or not.
    :type persistence: an array of (dimension, (birth, death)) or an array of (birth, death).
    :param max_intervals: maximal number of intervals to display.
        Selected intervals are those with the longest life time. Set it
        to 0 to see all. Default value is 1000.
    :type max_intervals: int.
    :param key: key function for sort algorithm.
    :type key: function or lambda.
    """
    if max_intervals > 0 and max_intervals < len(persistence):
        warnings.warn(
            "There are %s intervals given as input, whereas max_intervals is set to %s."
            % (len(persistence), max_intervals)
        )
        # Sort by life time, then takes only the max_intervals elements
        return sorted(persistence, key=key, reverse=True)[:max_intervals]
    else:
        return persistence


@lru_cache(maxsize=1)
def _matplotlib_can_use_tex():
    """This function returns True if matplotlib can deal with LaTeX, False otherwise.
    The returned value is cached.
    """
    try:
        from matplotlib import checkdep_usetex

        return checkdep_usetex(True)
    except ImportError as import_error:
        warnings.warn(f"This function is not available.\nModuleNotFoundError: No module named '{import_error.name}'.")


import gudhi as gd
import matplotlib.pyplot as plt
import os
import numpy as np
import errno
from sys import path
import warnings

def plot_persistence_diagram(
    persistence=[],
    persistence_file="",
    alpha=1.0,
    band=0.0,
    max_intervals=1000000,
    inf_delta=0.1,
    colormap=None,
    axes=None,
    fontsize=16,
    legend = None,
    greyblock=False,
):
    r"""This function plots the persistence diagram from persistence values
    list, a np.array of shape (N x 2) representing a diagram in a single
    homology dimension, or from a `persistence diagram <fileformats.html#persistence-diagram>`_ file`.
    :param persistence: Persistence intervals values list. Can be grouped by dimension or not.
    :type persistence: an array of (dimension, (birth, death)) or an array of (birth, death)
    :param persistence_file: A `persistence diagram <fileformats.html#persistence-diagram>`_ file style name
        (reset persistence if both are set).
    :type persistence_file: string
    :param alpha: plot transparency value (0.0 transparent through 1.0
        opaque - default is 0.6).
    :type alpha: float
    :param band: band (not displayed if :math:`\leq` 0. - default is 0.)
    :type band: float
    :param max_intervals: maximal number of intervals to display.
        Selected intervals are those with the longest life time. Set it
        to 0 to see all. Default value is 1000000.
    :type max_intervals: int
    :param inf_delta: Infinity is placed at :code:`((max_death - min_birth) x
        inf_delta)` above :code:`max_death` value. A reasonable value is
        between 0.05 and 0.5 - default is 0.1.
    :type inf_delta: float
    :param legend: Display the dimension color legend. Default is None, meaning the legend is displayed if dimension
        is specified in the persistence argument, and not displayed if dimension is not specified.
    :type legend: boolean or None
    :param colormap: A matplotlib-like qualitative colormaps. Default is None
        which means :code:`matplotlib.cm.Set1.colors`.
    :type colormap: tuple of colors (3-tuple of float between 0. and 1.)
    :param axes: A matplotlib-like subplot axes. If None, the plot is drawn on
        a new set of axes.
    :type axes: `matplotlib.axes.Axes`
    :param fontsize: Fontsize to use in axis.
    :type fontsize: int
    :param greyblock: if we want to plot a grey patch on the lower half plane for nicer rendering. Default True.
    :type greyblock: boolean
    :returns: (`matplotlib.axes.Axes`): The axes on which the plot was drawn.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
        from matplotlib import rc

        if _gudhi_matplotlib_use_tex and _matplotlib_can_use_tex():
            plt.rc("text", usetex=True)
            plt.rc("font", family="serif")
        else:
            plt.rc("text", usetex=False)
            plt.rc("font", family="DejaVu Sans")

        # By default, let's say the persistence is not an array of shape (N x 2) - Can be from a persistence file
        nx2_array = False
        if persistence_file != "":
            if path.isfile(persistence_file):
                # Reset persistence
                persistence = []
                diag = read_persistence_intervals_grouped_by_dimension(persistence_file=persistence_file)
                for key in diag.keys():
                    for persistence_interval in diag[key]:
                        persistence.append((key, persistence_interval))
            else:
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), persistence_file)

        try:
            persistence, nx2_array = _array_handler(persistence)
            persistence = _limit_to_max_intervals(
                persistence, max_intervals, key=lambda life_time: life_time[1][1] - life_time[1][0]
            )
            min_birth, max_death = 0.2, 1.4
        except IndexError:
            min_birth, max_death = 0.2, 1.4
            pass

        delta = (max_death - min_birth) * inf_delta
        # Replace infinity values with max_death + delta for diagram to be more
        # readable
        
        infinity = max_death + delta
        axis_end = max_death + delta / 2
        axis_start = min_birth - delta

        if axes is None:
            _, axes = plt.subplots(1, 1)
        if colormap is None:
            colormap = plt.cm.Set2.colors
        # bootstrap band
        if band > 0.0:
            x = np.linspace(axis_start, infinity, 1000)
            axes.fill_between(x, x, x + band, alpha=alpha, facecolor="red")
        # lower diag patch
        if greyblock:
            axes.add_patch(
                mpatches.Polygon(
                    [[axis_start, axis_start], [axis_end, axis_start], [axis_end, axis_end]],
                    fill=True,
                    color="lightgrey",
                )
            )
        # line display of equation : birth = death
        #axes.plot([axis_start, axis_end], [axis_start, axis_end], linewidth=1.0, color="k")

        x = [birth for (dim, (birth, death)) in persistence]
        y = [death if death != float("inf") else infinity for (dim, (birth, death)) in persistence]
        c = [colormap[dim] for (dim, (birth, death)) in persistence]

        axes.scatter(x, y, alpha=alpha, color=c)
        if float("inf") in (death for (dim, (birth, death)) in persistence):
            # infinity line and text
            #axes.plot([axis_start, axis_end], [infinity, infinity], linewidth=1.0, color="k", alpha=alpha)
            # Infinity label
            yt = axes.get_yticks()
            yt = yt[np.where(yt < axis_end)]  # to avoid plotting ticklabel higher than infinity
            yt = np.append(yt, infinity)
            ytl = ["%.3f" % e for e in yt]  # to avoid float precision error
            ytl[-1] = r"$+\infty$"
            
            #axes.set_yticks(yt)
            #axes.set_yticklabels(ytl)

        if legend is None and not nx2_array:
            # By default, if persistence is an array of (dimension, (birth, death)), display the legend
            legend = False

        if legend:
            dimensions = list({item[0] for item in persistence})
            axes.legend(
                handles=[mpatches.Patch(color=colormap[dim], label=str(dim)) for dim in dimensions],
                loc="lower right",
            )

        axes.set_xlabel("Birth", fontsize=fontsize)
        axes.set_ylabel("Death", fontsize=fontsize)
        axes.set_title("Persistence diagram", fontsize=fontsize)
        # Ends plot on infinity value and starts a little bit before min_birth

        axes.axis([axis_start, axis_end, axis_start, infinity + delta / 2])
        return axes

    except ImportError as import_error:
        warnings.warn(f"This function is not available.\nModuleNotFoundError: No module named '{import_error.name}'.")


class Preprocess:

    def __init__(self, patient_data):
        
        self.patient_data = patient_data

    def preprocess(self, ref):

        self.df = pd.read_pickle(self.patient_data)
        self.genes = pd.read_csv("/dir/Genes.csv", sep=';')
        self.genes = list(self.genes["Genes"])
        self.data = self.df[self.df.columns.intersection(self.genes)]        
        
        self.data.index = [ref] * self.data.shape[0]
        print(self.data.shape)

        return self.data

class Prediction:

    def __init__(self, resistance, sensitive):
        
        self.resistance = resistance
        self.sensitive = sensitive

    def split(self):
        """"
        Split gene expression dataset into train and test subsets
        :return: train expression & labels; test expression & labels
        """
        self.all_patients = pd.concat([self.resistance, self.sensitive])
        self.labelencoding = preprocessing.LabelEncoder()
        self.labels = self.labelencoding.fit_transform(self.all_patients.index.to_list())

        self.train, self.test = train_test_split(self.all_patients, test_size=0.2, shuffle=True)
        self.train_labs = self.labelencoding.fit_transform(self.train.index.to_list())
        self.test_labs = self.labelencoding.fit_transform(self.test.index.to_list())
        self.train, self.test = np.array(self.train), np.array(self.test)
        return self.train, self.test

    def intergene_correlation_measure(self):
        """"
        Calculate the distance correalation measures across training patients
        :return: distance correlation matrix across training patients
        """
        
        self.num_genes = self.train.shape[1]
        self.dist = np.zeros((self.num_genes, self.num_genes))
        
        for i in tqdm(range(self.num_genes)):
            for j in range(i + 1, self.num_genes):
                self.dist[i, j] = dcor.distance_correlation(self.train[:, i], self.train[:, j])#Distance Correlations 
                
        self.dist = self.dist + self.dist.T + np.eye(self.num_genes)
        self.dist = 1 - self.dist

    def intergene_signed_TOM(self):
        """"
        Read a calculated signed-TOM across training patients
        :return: signed-TOM matrix across training patients
        """
        self.dist = pd.read_csv(self.signedTOM, index_col=0)
        return self.dist

    def patient_correlation_measure(self, F):
         """
        Calculate the adjusted distance correlation measure for individual patients  
        :return: adjusted distance correlaion measure
        """
        F = F.T
        self.patient_dist = np.zeros((self.num_genes, self.num_genes))
        
        for i in range(self.num_genes):
            for j in range(i+1, self.num_genes):
                self.patient_dist[i,j] = self.dist[i,j] + math.sqrt(F[i] + F[j])/10 #smaller values are shrunk and larger values are inflated 
        
        self.patient_dist = self.patient_dist + self.patient_dist.T + np.eye(self.num_genes)

        return self.patient_dist
    

    def filter_persistence(self, persistence_intervals, top_percentile=50):
        """
        Filters based on the lifespan and select the top percentile persisting features
        :return: filtered persistence betti-numbers based on lifespan
        """
        dict_birth_death = list(enumerate(persistence_intervals))
        filter = [i[1][1][1] - i[1][1][0] for i in dict_birth_death ]
        
        dict_lifespan = dict(enumerate(filter))
        dict_sorted = dict(sorted(dict_lifespan.items(), key=lambda x: x[1])) 

        subset = math.floor(len(dict_birth_death) * int(top_percentile) / 100) #select the number of topological features to select
        subset_dict = {k: dict_sorted[k] for k in list(dict_sorted)[:subset]}
        keys = [i for i in subset_dict.keys()]
        
        persistence_pairs = [persistence_intervals[i] for i in keys] 

        return (persistence_pairs)

    def PeristentHomology(self, patients, lbl, top_percentile=30): #default take the 50% most persistent Betti numbers
        """Computes the persistence of the simplicial complex. Saves as persistent betti pairs
        ":return:Persistence diagram for each patient  
        """
        count = 0
        
        for patient in tqdm(range(len(patients))):
            patient_exprs = patients[patient]

            if lbl == 'train':
                patient_label = self.train_labs[patient]
            elif lbl == 'test':
                patient_label = self.test_labs[patient]

            self.distance_matrix = self.patient_correlation_measure(patient_exprs) #Weights used include per-patient gene expressions
            
            self.rips_complex = gd.RipsComplex(distance_matrix=self.distance_matrix, max_edge_length=np.inf, sparse=1)
            self.rips_complex = self.rips_complex.create_simplex_tree(max_dimension=1)
            self.rips_complex.collapse_edges()
            self.rips_complex.expansion(3)
            
            diag = self.rips_complex.persistence()
            diag = [i for i in diag if i[0] > 0] #remove Betti-0

            try:
                diag_filter = self.filter_persistence(diag)
                diag_plot = plot_persistence_diagram(diag_filter, greyblock=False, legend=None,
                                                     colormap=None)

                image = '/dir/Images/' + str(patient_label) + '_' + str(count) + '.png'
                diag_plot.figure.savefig(image, dpi=600)

                count += 1
            
            except:
                continue

if __name__ == '__main__':

    resistance = Preprocess("/dir/Resistance.pkl") #Format input gene expression path
    resistance_patients = resistance.preprocess("R") 

    sensitive = Preprocess("/dir/Sensitive.pkl") #Format input gene expression path
    sensitive_patients = sensitive.preprocess("S") 

    prediction = Prediction(resistance_patients, sensitive_patients)
    X_train, X_test = prediction.split() 

    prediction.intergene_correlation_measure() 
    prediction.intergene_signed_TOM() 
    
    Train_samples = prediction.PHomology(X_train, 'train') 
    Test_samples = prediction.PHomology(X_test, 'test') 

    os.system("mkdir /dir/Images") #Format output folder path
    process_images('/dir/Images') 