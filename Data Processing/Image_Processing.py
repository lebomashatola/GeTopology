'''
Image segmentation and image tiling of WSI obtained from https://portal.gdc.cancer.gov
'''

import os 
import openslide
from shutil import copyfile

class Image:

    def __init__(self, inputImages, OutputImages, DestinationImages):
        """ 
        Initialise input, intermediate and output directory of WSI 
        """
        self.inputImages = inputImages
        self.OutputImages = OutputImages
        self.DestinationImages = DestinationImages

    def segmentation_tiling(self):
        """ 
        Perform image segmentation and tiling 
        :return:For each sample/folder, tiled WSI's   
        """
        for file in os.listdir(self.inputImages):
            self.filename = os.fsdecode(file)
            self.read_data = self.inputImages + '/' + self.filename
            
            try:
                self.files = os.listdir(self.read_data)
            except:
                continue
                    
            for i in self.files:
                if i.endswith(".svs"):
                    self.image_in = self.inputImages + '/' + file  + '/' + i
                    os.system('python3 pyhist.py --method "graph" --patch-size 64 --output-downsample 16 ' +
                                '--content-threshold 0.8 --k-const 1024 --save-patches --save-tilecrossed-image --output ' 
                                + self.OutputImages + '/' + ' ' + self.image_in)
 
    def move_images(self):
         """ 
        Move separated tiled images into a single folder of same phenotype 
        :return:Folder with all tiles images generated from each patient  
        """
        for file in os.listdir(self.OutputImages):
            self.filename = os.fsdecode(file)
            self.read_data = self.OutputImages + '/' + self.filename
            
            for image in os.listdir(self.read_data):
                if ("_tiles" in image):
                    self.soc = self.read_data + '/' + image
                    imgs = os.listdir(self.soc)
                    
                    for i in imgs:
                        os.system('mv ' + self.soc + '/' + i + ' ' + self.destination)

if __name__ == '__main__':

    while True:
        
        print(" GeTopology (v0.1) \n Welcome to WSI Image Processing!")
        print("###################################")
        opt = int(input('Select option: \n 1. Process Stage 1 \n 2. Process Stage 4 \n 3. Process both stages \n 4. Exit \n : '))

        if opt == 1:
            Read_Imgs_Stg1 = Image("/dir/Images_Stage_1", "/dir/stage_1", "/dir/stg_1")
            Read_Imgs_Stg1.segmentation_tiling()
            Read_Imgs_Stg1.move_images()
            print('Process Complete! \n')

        if opt == 2:
            Read_Imgs_Stg4 = Image("/dir/Images_Stage_4", "/dir/stage_4", "/dir/stg_4")
            Read_Imgs_Stg4.segmentation_tiling()
            Read_Imgs_Stg4.move_images()
            print('Process Complete! \n')
    
        if opt == 3:
            Read_Imgs_Stg1 = Image("/dir/Images_Stage_1", "/dir/stage_1", "/dir/stg_1")
            Read_Imgs_Stg1.segmentation_tiling()
            Read_Imgs_Stg1.move_images()
            print('Complete Stage 1 \n')

            Read_Imgs_Stg4 = Image("/dir/Images_Stage_4", "/dir/stage_4", "/dir/stg_4")
            Read_Imgs_Stg4.segmentation_tiling()
            Read_Imgs_Stg4.move_images()
            print('Complete Stage 4 \n')

        if opt == 4:
            break
