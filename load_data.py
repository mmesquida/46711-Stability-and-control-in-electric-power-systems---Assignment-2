import scipy.io as sio
import os


# This function loads the provided data for the assignement
def load_data(question):
    
    current_path = os.getcwd()
    additional_path = f'/Assignment_data/system_{question}.mat'
    file_path = current_path + additional_path
    
    data = sio.loadmat(file_path)
    return data

