
"""
Transfer Learning for Image Classification for Cancer Stages
"""

from tensorflow.keras.layers import Dense, Flatten, BatchNormalization, Conv2D, MaxPool2D
from tensorflow.keras.preprocessing.image import ImageDataGenerator
from tensorflow.keras.metrics import categorical_crossentropy
from sklearn.metrics import plot_confusion_matrix
from tensorflow.keras.models import Sequential
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import Model
from sklearn.metrics import confusion_matrix
from IPython.display import Image
from keras import backend as K
from tensorflow import keras
import tensorflow as tf
import numpy as np
import itertools
import glob
import os

"""
Run code for Nvidia GPU application
"""
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
config = tf.compat.v1.ConfigProto()
gpu_options = tf.compat.v1.GPUOptions(per_process_gpu_memory_fraction = 0.6)
sess = tf.compat.v1.Session(config = tf.compat.v1.ConfigProto(gpu_options = gpu_options))
tf.compat.v1.keras.backend.set_session(sess)


class Prediction_Histology:

    def __init__(self, image_gen, shape, model_call, rm, label):

        self.image_gen = image_gen
        self.shape = shape
        self.model_call = model_call
        self.rm = rm
        self.label = label

    def preprocess(self):
        
        def plot_confusion_matrix(cm, classes, normalize=False,
                                  title='Confusion matrix', cmap=plt.cm.Blues):
            """
            Outputs model performance metrics applicable to classification task
            :return:TF,TN,FP,FN
            """
            plt.imshow(cm, interpolation='nearest', cmap=cmap)
            plt.title(title)
            plt.colorbar()
            tick_marks = np.arange(len(classes))
            plt.xticks(tick_marks, classes, rotation=45)
            plt.yticks(tick_marks, classes)

            if normalize:
                cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
                print("Normalized confusion matrix")
            else:
                print('Confusion matrix, without normalization')

            print(cm)

            thresh = cm.max() / 2.
            for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
                plt.text(j, i, cm[i, j],
                         horizontalalignment="center",
                         color="white" if cm[i, j] > thresh else "black")

            plt.tight_layout()
            plt.ylabel('True label')
            plt.xlabel('Predicted label')
            plt.show()

        """
        Read training, validation and testing images 
        """

        os.chdir('/dir/images')

        train_path = '/images/train'
        valid_path = '/images/validate'
        test_path = '/images/test'
        
        train_batch = ImageDataGenerator(preprocessing_function=self.image_gen) \
        .flow_from_directory(directory=train_path, target_size=self.shape,
        classes=['partial', 'responsive', 'unresponsive'], batch_size=15)

        test_batch = ImageDataGenerator(preprocessing_function=self.image_gen) \
        .flow_from_directory(directory=test_path, target_size=self.shape,
        classes=['partial', 'responsive', 'unresponsive'], batch_size=15)

        validate_batch = ImageDataGenerator(preprocessing_function=self.image_gen) \
        .flow_from_directory(directory=valid_path, target_size=self.shape,
        classes=['partial', 'responsive', 'unresponsive'], batch_size=15)

        imgs, labels = next(train_batch)
        x = self.model_call.layers[self.rm].output

        config = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=12,
                                inter_op_parallelism_threads=6,
                                allow_soft_placement=True,
                                device_count = {'CPU': 12})

        session = tf.compat.v1.Session(config=config)

        os.environ["6"] = "6"
        os.environ["10"] = "30"
        os.environ["KMP_SETTINGS"] = "1"
        os.environ["KMP_AFFINITY"]= "granularity=fine,verbose,compact,1,0"

        output = Dense(units=3, activation='softmax')(x)
        model = Model(inputs=self.model_call.input, outputs=output)

        model.compile(optimizer=Adam(learning_rate=0.0001),
                            loss='categorical_crossentropy', metrics=['accuracy'])

        history = model.fit(x=train_batch, validation_data=validate_batch, epochs=2, verbose=2)
        predictions = model.predict(x=test_batch, verbose=0)

        plt.plot(history.history['accuracy'], label='Training')
        plt.plot(history.history['val_accuracy'], label='Validation')

        plt.xlabel('Number of Epochs')
        plt.ylabel('Accuracy(%)')
        plt.tight_layout()

        plt.legend()
        plt.show()

        plt.plot(history.history['loss'], label='Training')
        plt.plot(history.history['val_loss'], label='Validation')

        plt.xlabel('Number of Epochs')
        plt.ylabel('Loss)')
        plt.tight_layout()

        plt.legend()
        plt.show()

        cm = confusion_matrix(y_true=test_batch.classes, y_pred=np.argmax(predictions, axis=1))
        cm_plot_labels = ['Responsive', 'Partial', 'Unresponsive']

        plot_confusion_matrix(cm, cm_plot_labels, title='Confusion Matrix')


if __name__ == '__main__':

    resent_preprocess = tf.keras.applications.resnet_v2.preprocess_input
    nasnet_preprocess = tf.keras.applications.nasnet.preprocess_input
    mobile_preprocess = tf.keras.applications.mobilenet.preprocess_input
    incept_preprocess = tf.keras.applications.inception_resnet_v2.preprocess_input

    resnet = tf.keras.applications.ResNet152V2()
    nesnet = tf.keras.applications.NASNetLarge()
    mobileNet = tf.keras.applications.mobilenet.MobileNet()
    incept = tf.keras.applications.InceptionResNetV2()

    hist_mobile = Prediction_Histology(mobile_preprocess, (224,224), mobileNet, -6, 'Mobile')
    hist_resnet = Prediction_Histology(resent_preprocess, (224,224), resnet, -2, 'ResNet')
    hist_nesnet = Prediction_Histology(nasnet_preprocess, (331,331), nesnet, -2, 'NesNet')
    hist_incept = Prediction_Histology(incept_preprocess, (299,299), incept, -2, 'InceptNet')

    while True:
        
        print(" GeTopology (v0.1) \n Welcome to Phenotype Prediction using Transfer Learning on WSI!")
        print("###################################")
        opt = input('Select option: \n 1. ResNet \n 2. NasNet \n 3. MobileNet \n 4. Incept \n 5. Exit \n : ')

        if opt == 1:
            folds = input('Enter number of cross-validation folds to compute: ')
            for j in range(folds):
                hist_resnet.process()
                print('fold number ' + j + ' complete!')

        if opt == 2:
            folds = input('Enter number of cross-validation folds to compute: ')
            for j in range(folds):
                hist_nesnet.process()
                print('fold number ' + j + ' complete!')

        if opt == 3:
            folds = input('Enter number of cross-validation folds to compute: ')
            for j in range(folds):
                hist_mobile.process()
                print('fold number ' + j + ' complete!')

        if opt == 4:
            folds = input('Enter number of cross-validation folds to compute: ')
            for j in range(folds):
                hist_incept.process()
                print('fold number ' + j + ' complete!')

        if opt == 5:
            break