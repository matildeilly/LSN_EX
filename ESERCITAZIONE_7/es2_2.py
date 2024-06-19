import numpy as np

fase = ["gas","liquid","solid"]

def vector_variance(V):
    dim = len(V)
    if dim == 1:
        return 0
    else:
        sum_val = np.sum(V)
        sum2_val = np.sum(np.square(V))
        return np.sqrt(1.0 / dim * (1.0 / dim * sum2_val - np.power(1.0 / dim * sum_val, 2)))

for name in fase:
    data = np.loadtxt('instant_en_'+ name + '_out.txt')
    file_out = open("Error_en_block_"+name+".txt", "w")

    U = data[:,1]
    

    for L in range (10,5000,10):
        mean = []
        N = len(U)/L

        for i in range(int(N)):
            mean.append(np.sum(U[i*L:(i+1)*L])/L)     

        line = str(L) + "\t" + str(vector_variance(mean)) + "\n"

        file_out.write(line)        

    file_out.close()    