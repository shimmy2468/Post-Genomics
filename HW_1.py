# Arbitrary matrices
A = [[1,0,1],
    [1,1,1],
    [2,2,1]]

B = [[2,1,0],
    [0,1,2],
    [2,0,0]]

#Function creation
def matrix_mult(Matrix_1, Matrix_2):
    result_matrix = [[0 for j in range(len(Matrix_2[0]))] for i in range(len(Matrix_1))]
    if len(Matrix_1[0]) != len(Matrix_2):
        return "Matrices not compatible"
    for i in range(len(Matrix_1)):
        for j in range(len(Matrix_2[0])):
            for k in range(len(Matrix_2)):
                result_matrix[i][j] += Matrix_1[i][k] * Matrix_2[k][j]
    
    return result_matrix

# Creating an incompatible matrix to test functionality
C =  [[1,2,3]]            
print(matrix_mult(A, C))


# Part 2.2 Seeing which is faster
import numpy as np
import time

start_time_my_function = time.perf_counter()
print(matrix_mult(A, B))
end_time_my_function = time.perf_counter()
elapsed_time_my_function = end_time_my_function - start_time_my_function
print(f"Time elapsed for my function: {elapsed_time_my_function}")

print("\n")

start_time = time.perf_counter()
print(np.dot(A, B))
end_time = time.perf_counter()
elapsed_time = end_time - start_time
print(f"Time elapsed for the numpy dot function: {elapsed_time}")

# My function is faster, but by a very small margin