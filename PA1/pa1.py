import pandas as pd
from scipy.stats import norm, chisquare
from math import sqrt


# def compute_association_statistics(n, p_pos, p_neg):
#     p = (p_pos + p_neg) / 2
#     s = (p_pos - p_neg) / (sqrt(2 / n) * sqrt(p * (1 - p)))
#     return s


def compute_frequency(df, n):
    if 1 in df.index and 2 in df.index:
        # print('b[1]', b[1])
        # print('b[2]', b[2])
        p = (df[1] + 2 * df[2]) / n
    elif 1 in df.index:
        # print('b[1]', b[1])
        p = (df[1]) / n

    elif 2 in df.index:
        # print('b[2]', b[2])
        p = ( 2 * df[2]) / n
    else:
        # print('none')
        p = 1

    return p



def compute_association_statistics(n, case, control):

    p_pos = compute_frequency(case, n)
    # print('p_pos = ', p_pos)
    p_neg = compute_frequency(control, n)
    # print('p_neg = ', p_neg)

    if (p_pos != 1) and (p_neg != 1):
        p = (p_pos + p_neg) / 2
        s = ((p_pos - p_neg) / (sqrt(2 / n) * sqrt(p * (1 - p) ) ) )

    else :
        s = 0


    return s





input_file = "SNP_status.txt"
# delimiter=' ':    space separated
# index_col = 0:    use the first column as row names
# usecols=list(range(0, 100000)): Return a subset of the columns, specified by column number [0, 1, 2, ... 100000]
# low_memory=False: use low memory
# a = read_csv(input_file, delimiter=' ', index_col=0, usecols=list(range(0, 100002)), low_memory=False,  nrows=10)
df = pd.read_csv(input_file, delimiter=' ', index_col=0, usecols=list(range(0, 100002)), low_memory=False)

alpha = 0.05
N = 2000
alpha_s = alpha/100000
threshold = -norm.ppf(alpha_s / 2)

# split the dataframne by Case and Control
case = df[df['Status'] == 'Case']
control = df[df['Status'] == 'Control']

# print(case)
# print(control)
p_values = []
p_values_significant = []
chi_square = []
for i in range(0,100000):
    # print(i)
    if i < 10:

        SNPname = 'SNP0000{}'.format(i)
    elif i < 100:
        SNPname = 'SNP000{}'.format(i)
    elif i < 1000:
        SNPname = 'SNP00{}'.format(i)
    elif i < 10000:
        SNPname = 'SNP0{}'.format(i)
    else:
        if i == 10000:
            print(10000)
        elif i == 20000:
            print(20000)
        elif i == 30000:
            print(30000)
        elif i == 40000:
            print(40000)
        elif i == 50000:
            print(50000)
        elif i == 60000:
            print(60000)
        elif i == 70000:
            print(70000)
        elif i == 80000:
            print(80000)
        elif i == 90000:
            print(90000)
        SNPname = 'SNP{}'.format(i)





    a = case[SNPname].value_counts()
    b = control[SNPname].value_counts()
    # print(a)
    # print(b)


    s = compute_association_statistics(N, a, b)

    # print (s, '\n\n\n')
    p_values.append([SNPname, s])
    c2 = s*s
    chi_square.append(c2)
    if (s > threshold):
        p_values_significant.append([SNPname, s])


chi_square.sort()
median = (chi_square[49999] + chi_square[50000]) / 2

for i in range(0,100000):
    chi_square.append(p_values[i][1] * p_values[i][1])
median2 = 49/81

gc = median / median2



with open('./output.txt', 'w') as output_file:
    output_file.write('<A>\n')
    for x in p_values:
        output_file.write(str(x[0]) + ':' + str(x[1]) + '\n')
    output_file.write('</A>'+ '\n')
    output_file.write('<B>'+ '\n')
    for x in p_values_significant:
        output_file.write(str(x[0]) + ':' + str(x[1]) + '\n')
    output_file.write('</B>'+ '\n')
    output_file.write('<C>'+ '\n')
    output_file.write('Lambda_gc:' + str(gc) + '\n')
    output_file.write('</C>' + '\n')


