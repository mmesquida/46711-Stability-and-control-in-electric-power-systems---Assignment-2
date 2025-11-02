from tabulate import tabulate
import numpy as np

def print_eigenvalues(lambda_A, f_d, zeta):
    
    column_headers = ['No.', 'Re{λ}', 'Im{λ}', 'f_d [Hz]', 'ζ']
    
    data = []
    for i in range(len(lambda_A)):
        row = [i+1, 
               f"{np.real(lambda_A[i]):.4f}",
               f"{np.imag(lambda_A[i]):.4f}" if np.imag(lambda_A[i]) != 0 else '−',
               f"{f_d[i]:.4f}" if np.imag(lambda_A[i]) != 0 else '−',
               f"{zeta[i]:.4f}" if np.imag(lambda_A[i]) != 0 else '−']
        data.append(row)
    
    table = tabulate(data, headers = column_headers, tablefmt = "fancy_grid")
    
    print('Eigenvalues, Frequency and Damping Ratio:')
    print(table)
 
    
def print_P_matrix(P, row_headers):

    column_headers = []
    for i in range(P.shape[1]):
        column_headers.append(f'λ{i+1}')
        
    data = np.empty(P.shape)
    for i in range(P.shape[0]):
        for j in range(P.shape[1]):
            data[i,j] = np.abs(P[i,j])
    
    table = tabulate(data, headers = column_headers, showindex = row_headers, tablefmt="fancy_grid", floatfmt = '.2f')
    
    print('Participation Matrix (Absolute Values):')
    print(table)
    