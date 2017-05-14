import pandas as pd
import numpy as np
from math import sqrt
from scipy.stats import norm
from sklearn import preprocessing



input_file = "gwas_data.txt"
# delimiter=' ':    space separated
# index_col = 0:    use the first column as row names
# usecols=list(range(0, 100000)): Return a subset of the columns, specified by column number [0, 1, 2, ... 100000]
# low_memory=False: use low memory
# a = read_csv(input_file, delimiter=' ', index_col=0, usecols=list(range(0, 9511)), low_memory=False,  nrows=10)
df = pd.read_csv(input_file, delimiter=' ', index_col=0, usecols=list(range(0, 9511)), low_memory=False)
(m,n) = df.shape

print (m,n)
# m is number of individuals
# n-1 is number of SNPs
df = df.astype(float)
matrix = df.as_matrix()
X = np.array(matrix)
X_scaled = preprocessing.scale(X)

df_scaled = df
phenotype = df['Phenotype']

# print (phenotype)
for i in range(0, n-1):
    if i %400 == 0:
        print(i)
    for j in range(0, m):
        df_scaled.iloc[j, i] = X_scaled[j][i]
SNPs = df_scaled.iloc[:,0:n-1]
# print (SNPs)
# print(df['Phenotype'].value_counts())

mu = phenotype / m

SNPs_transpose = SNPs.transpose()
beta =  SNPs_transpose.dot(phenotype)
e = phenotype - mu - SNPs.dot(beta)
et = e.transpose()
sigma2 = (et.dot(e)) / (m-2)
sigma = sqrt(sigma2)
S = beta * sqrt(m) / sigma

# print(df_scaled)
xi = SNPs.iloc[0]
xj = SNPs.iloc[1]
# calculate norm of each individual
a = np.sqrt(np.square(SNPs).sum(axis=1))
kin_matrix = pd.DataFrame(np.nan, index= list(range(0, m)), columns=list(range(0, m)))

# kin_matrix = pd.DataFrame(np.nan, index=list(range(0, m)), columns=list(range(0,m))))

i_norm = np.sqrt(np.square(SNPs).sum(axis=1))
for i in range (0,m):
    for j in range(0, m):
        dot_product = SNPs.iloc[i].dot(SNPs.iloc[j])
        kin_matrix[i][j] = dot_product / (m * i_norm[i] * i_norm[j])
# print(kin_matrix)
a = kin_matrix.mean()
# print(a)
avg_k = a.mean()
# print(avg_k)


alpha = 0.05
alpha_s = alpha / (n-1)
threshold = -norm.ppf(alpha_s / 2)
print(threshold)
p_values_significant = []

for i in range(0,n-1):
    if (S[i] > threshold or S[i] < threshold * (-1) )  :
        p_values_significant.append([S.index[i], S[i]])



# print(p_values_significant)



with open('./output.txt', 'w') as output_file:
    output_file.write('<A>\n')
    for i in range(0, n-1):
        output_file.write(str(S.index[i])+ ':' + str(S[i]) + '\n')
    output_file.write('</A>' + '\n')
    output_file.write('<B>'+ '\n')
    for x in p_values_significant:
        output_file.write(str(x[0]) + ':' + str(x[1]) + '\n')
    output_file.write('</B>' + '\n')
    output_file.write('<C>' + '\n')
    output_file.write('AVG_K:' + str(avg_k) + '\n')
    output_file.write('</C>' + '\n')


