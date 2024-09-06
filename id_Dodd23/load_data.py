import numpy as np
import os
data_folder =  str(os.path.dirname(os.path.dirname(__file__))) +"/Data/"
test_data_file = data_folder + "test_data.csv"
# Load Fits and Groups
named_Groups = np.load(data_folder+"original_Dodd23_Groups.npy")
group_mean = np.load(data_folder+"original_Dodd23_Group_Mean_ELzLp.npy")
group_covar = np.load(data_folder+"original_Dodd23_Group_Covar_ELzLp.npy")
