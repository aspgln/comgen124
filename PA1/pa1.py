import pandas as pd
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


    elif 1 in b.index:
        # print('b[1]', b[1])
        p = (df[1]) / n

    elif 2 in b.index:
        # print('b[2]', b[2])
        p = ( 2 * df[2]) / n
    else:
        # print('none')
        p = 0

    return p



def compute_association_statistics(n, case, control):

    p_pos = compute_frequency(case, n)
    print('p_pos = ', p_pos)
    p_neg = compute_frequency(control, n)
    print('p_neg = ', p_neg)

    if (p_pos !=1) and (p_neg != 1):
        p = (p_pos + p_neg) / 2
        s = (p_pos - p_neg) / (sqrt(2 / n) * sqrt(p * (1 - p)))

    else :
        s = 0

    return s


N = 10
alpha = 0.05



input_file = "SNP_status.txt"
# delimiter=' ':    space separated
# index_col = 0:    use the first column as row names
# usecols=list(range(0, 100000)): Return a subset of the columns, specified by column number [0, 1, 2, ... 100000]
# low_memory=False: use low memory
# a = read_csv(input_file, delimiter=' ', index_col=0, usecols=list(range(0, 100002)), low_memory=False,  nrows=10)
df = pd.read_csv(input_file, delimiter=' ', index_col=0, usecols=list(range(0, 100002)), low_memory=False, nrows=2000)

N = 2000

# split the dataframne by Case and Control
case = df[df['Status'] == 'Case']
control = df[df['Status'] == 'Control']

print(case)
print(control)

a = case['SNP00000'].value_counts()
b = control['SNP00000'].value_counts()
print(a)
print(b)


s = compute_association_statistics(N, a, b )
print (s)

