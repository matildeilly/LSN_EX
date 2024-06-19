import numpy as np

fase = ["gas","liquid","solid"]

for name in fase:
    data = np.loadtxt('instant_en_'+ name + '_out.txt')
    file_out = open("Autocorrelation_"+name+".txt", "w")

    m = data[:,1]

    t_max = len(data)

    m_sum = np.sum(m)
    m_sum2 = np.sum(m * m)

    den = 1.0/t_max * m_sum2 - (1.0/t_max * m_sum)**2

    for t in range(t_max):
        a = np.sum(m[:t_max-t] * m[t:])
        b = np.sum(m[:t_max-t])
        c = np.sum(m[t:])

        num = 1.0/(t_max-t) * (a - 1.0/(t_max-t) * b * c)

        line = str(t) + "\t" + str(num/den) + "\n"

        file_out.write(line)        

    file_out.close()