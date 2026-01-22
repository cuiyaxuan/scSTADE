# scSTADE : An Integrated and Robust Deep Learning Framework for Denoising and Analyzing Single-Cell Spatial Transcriptomics
![image](https://github.com/cuiyaxuan/DenoiseST/blob/master/Image/%E5%B9%BB%E7%81%AF%E7%89%871.png)
## Tip: To facilitate researchers' usage, examples of our project can be run in the demo example folder's IPython notebooks (after configuring the environment dependencies as described in the README).  Additionally, we have set up a web tutorial to make it easier for users to access. [Website Link](https://denoisest-tutorial1.readthedocs.io/en/latest/). <br>


##### Due to the protocol issues of various space technology platforms, the data format is very different, and various platforms do not provide morphological images. For the convenience of users, we have changed the way of reading data to make it easier to use.<br>
##### For specific details, please refer to the file "1_Example_lowresolve_test.ipynb".Instructions for running other types of spatial transcriptomics data can be found in the IPython notebooks under the "full" directory.<br>
##### If your current system environment does not support both R and Python, you can refer to the User Manual for more detailed instructions. Alternatively, you can create two separate environments to run the program. <br>


# Identifying spatial domains
##### Using python virtual environment with conda. Please create a Pytorch environment, install Pytorch and some other packages, such as "numpy","pandas", "scikit-learn" and "scanpy". See the requirements.txt file for an overview of the packages in the environment we used to produce our results. Alternatively, you can install the environment dependencies in the following sequence to minimize environment conflicts. <br>

```R
conda create -n pipeline
source activate pipeline

conda search r-base
conda install r-base=4.2.0
conda install python=3.8

conda install conda-forge::gmp
conda install conda-forge::r-seurat==4.4.0
conda install conda-forge::r-hdf5r
conda install bioconda::bioconductor-sc3

conda install conda-forge::pot
conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
conda search -c conda-forge r-mclust
conda install -n pipeline -c conda-forge r-mclust=5.4.10

pip install scanpy
pip install anndata==0.8.0
pip install pandas==1.4.2

conda install -c conda-forge rpy2=3.5.1
pip install scikit-learn==1.1.1
pip install scipy==1.8.1
pip install tqdm==4.64.0
pip install scikit-misc

```

##### to install some R packages. <br>
```R
import rpy2.robjects as robjects
robjects.r('''
install.packages('foreach')
install.packages('doParallel')
           ''')
```



# FVG:identifying functionally variable genes

##### Then, we execute the FVG model in the R environment <br>
##### First, cd /home/.../DenoiseST-main/FVG <br>

```R
conda create -n r4
source activate r4

conda search r-base
conda install r-base=4.2.0
conda install conda-forge::r-seurat==4.4.0
conda install conda-forge::r-hdf5r
onda install r-devtools
```

##### Using R virtual environment with conda <br>
```R

library("devtools")
devtools::install_github("shaoqiangzhang/DEGman")

install.packages("hdf5r")
install.packages('philentropy')
install.packages('dplyr')
install.packages('foreach')
install.packages('parallel')
install.packages('doParallel')


```




