
from tensorflow.keras.preprocessing.image import ImageDataGenerator
from tensorflow.keras.metrics import categorical_crossentropy
from tensorflow.keras.models import Sequential
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import Model
from tensorflow import keras
import tensorflow as tf
import pandas as pd
import numpy as np
import itertools
import glob
import os

class Integrate:

    def __init__(self, images, topology):
        """
        Initialise input data for topological summaries and images
        """
        self.images = images
        self.topology = topology

    def MobileNet(self):
        """
        Load pre-trained Mobile Net mobile and remove last 3 layers for concatenation
        """
        self.mobileNet = tf.keras.applications.mobilenet.MobileNet()
        self.x = mobileNet.layers[-3].output
        model.add(Dense(units=64, activation='softmax'))(self.x)

    
    def MLP(X_train):
        """
        Build a sequential neural network with 3 layers
        """
        model = Sequential()
        model.add(Dense(128, activation='relu', input_shape=(self.X_train.shape[1], )))
        model.add(Dense(64, activation='relu'))
        model.add(Dense(16, activation='relu'))
        
        return(model)

    def image_data(self):
        """
        Load training and testing images 
        """
        self.input_shape = (224, 224)
        self.batch_size = 8

        self.train_ds = image_dataset_from_directory(
            self.images + "/training", label_mode="binary", shuffle=True, subset=None,
            image_size=input_shape, batch_size=batch_size)

        self.test_ds = image_dataset_from_directory(
            self.images + "/testing", label_mode="binary", shuffle=True,
            subset=None, image_size=input_shape, batch_size=batch_size)


    def process(self):
        """
        Load training and testing topological summaries 
        """
        self.train = pd.read_csv(self.topology + '/train')
        self.train_labs = pd.read_csv(self.topology + '/train_labs')
        self.test = pd.read_csv(self.topology + '/test')
        self.test_labs = pd.read_csv(self.topology + '/test_labs')
        
    def input_models(self):
        """
        Form a hybrid model using image and topological summaries
        """
        self.mlp_model = MLP(self.train)
        self.mobileNet_model = MobileNet()
        self.combinedInput = concatenate([self.mlp_model.output, self.MobileNet.output])
        
        self.x_layer = Dense(9, activation="relu")(self.combinedInput)
        self.x_layer = Dense(1, activation="sigmoid")(self.x_layer)
        self.model = Model(inputs=[self.mlp_model.input, self.cnn_model.input], outputs=self.x_layer)
        plot_model(self.model, to_file='demo.pdf', show_shapes=True)

        self.model.compile(optimizer=keras.optimizers.Adam(1e-3),loss="binary_crossentropy", metrics=["accuracy", "loss"])
        
        self.model.fit([self.train, self.train_ds], self.train_labs, epochs=1000, batch_size=8, verbose=True)
        self.unseen = self.model.evaluate([self.test, self.test_ds], self.test_labs)

if __name__ == '__main__':
    
    while True:
        
        print(" GeTopology (v0.1) \n Welcome to Phenotype Prediction using Bimodal Learning!")
        print("###################################")

        dir_images = input("Enter file path for training and testing samples for cancer stage image data \n : ")
        dir_topology = input("Enter file path for training and testing samples for topological summaries data \n : ")
        
        initialise = Integrate(dir_images, dir_topology)
        initialise.MobileNet()
        initialise.MLP()
        initialise.image_data()
        initialise.process()
        intialise.input_models()

        print('Process Complete !')
        break


        