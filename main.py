import matplotlib.pyplot as plt
from math import *

#Дано:
M = [78, 80, 82, 83, 84, 86]
cf = [0.35, 2.27, 11.56, 11.52, 56.9, 17.40]
for i in range(len(cf)):
    cf[i] = cf[i] / 100
print("SUM CF =", sum(cf))
print("cf =", *cf)
q_0 = 1.1
print("q_0 =", q_0)
N = 50
print("N =", N)
f = 20
print("f =", f)
n_component = 3
print("n_component =", n_component+1)
M_star = (M[3] + M[4]) / 2
print("M* =", M_star)
# Расчет параметров каскада
g = [pow(q_0, M_star - M[i]) for i in range(len(M))] #Расчет g
print("g =", *g)
sum_cp = 0
sum_cw = 0
iter = len(cf)
for i in range(iter):
        sum_cw += cf[i] * (pow(g[i], N + 1 - f) - 1) / (pow(g[i], N + 1) - 1)
        sum_cp += cf[i] * (1 - pow(g[i], -f)) / (1 - pow(g[i], - N - 1))

#Расчет концентраций cf и cp
cp = [cf[i] * (1 - pow(g[i], -f)) / ((1 - pow(g[i], - N - 1)) * sum_cp) for i in range(iter)]
cw = [cf[i] * (pow(g[i], N + 1 - f) - 1) / (pow(g[i], N + 1) - 1) / sum_cw for i in range(iter)]
"""
cw = []
cp = []
for i in range(iter):
        cp.append(cf[i] * (1 - pow(g[i], -f)) / (1 - pow(g[i], -N - 1)) / sum_cp)
        cw.append(cf[i] * (pow(g[i], N + 1 - f) - 1) / sum_cw)
"""
print("SUM CP", sum(cp))
print("SUM CW", sum(cw))
F_P = 1 / sum_cp
print("F/P =", F_P)
P_F = sum_cp
print("P_F =", P_F)
W_P = sum_cw / sum_cp
print("W/P =", W_P)
W_F = sum_cw
print("W/F =", W_F)
print("cp =", *cp)
print("cw =", *cw)

L = [0]*iter
L_is = [0]*iter
for i in range(iter):
    L[i] = []
    for s in range(1, N+1):
        if(s < f):
            L[i].append(cw[i] * W_P * g[i] * (pow(g[i], s) - 1) / (g[i] - 1))
        else:
            L[i].append(cp[i] * g[i] * (1 - pow(g[i], s - N - 1)) / (g[i] - 1))

for i in range(iter):
    L_is[i] = []
    for s in range(N):
        L_is[i].append(L[i][s] * (g[i] + 1) / g[i])

print("L =", *L_is)
print("L' =", *L)

L_smm = 0
L_sel = []
for s in range(N):
    summ = 0
    for i in range(iter):
        summ += L[i][s] * (g[i] + 1) / g[i]
    L_sel.append(summ)
L_smm = sum(L_sel)
print(' summ( L_s/P ) =', L_smm)
print("L_sel =", *L_sel)
Rp_nk = cp[n_component] / cp[n_component + 1]
Rw_nk = cw[n_component] / cw[n_component + 1]
Rf_nk = cf[n_component] / cf[n_component + 1]
L_sum_test = []
for i in range(iter):
    temp_1 = cp[i] * log(Rp_nk) + W_P * cw[i] * log(Rw_nk) - F_P * cf[i] * log(Rf_nk)
    L_sum_test.append(temp_1 * (g[i] + 1) / log(g[n_component]) / (g[i] - 1))
print("L_sum_test", sum(L_sum_test))


C = [0] * iter
for i in range(iter):
    C[i] = []
    for s in range(N):
        C[i].append(L_is[i][s] / L_sel[s])

temp = 0
for i in range(iter):
    temp += C[i][5]
    print("C ", i + 1, C[i][24])
print(temp)
print("Степень извлечния целевого компонента =", cp[n_component] * sum_cp / cf[n_component])

n_array = [i for i in range(N)]
fig = plt.figure()
plt.title('L_s/P')
plt.plot(n_array, L_sel)
plt.scatter(n_array, L_sel)
plt.xlabel("s")
plt.ylabel('L(s)/P')
plt.grid()
plt.savefig('L.png', format='png', dpi=300)
plt.show()

Simbols = ['o','*','x','s','.', '^']
fig = plt.figure()
plt.title('$C_i$')
for i in range(iter):
    plt.plot(range(N), C[i], linewidth = 0.8)
    plt.scatter(range(N), C[i],marker=Simbols[i], s=15, label='Компонент '+str(i+1))
plt.xlabel("$s$")
plt.ylabel('$C_{i,s}$')
plt.legend()
plt.grid()
plt.savefig('C_i.png', format='png', dpi=300)
plt.show()

