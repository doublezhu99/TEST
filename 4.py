import numpy as np

# 设置需求、起始库存、最小库存和最大库存
D = np.loadtxt('C:/Users/jingcon/Documents/MATLAB/data_test.csv')
I0 = 20
smin = 0
smax = 40

# 初始化
f = float('inf') # 设置目标函数初始值为正无穷
S = 0 # 设置初始再订货点
s = I0 # 设置初始库存
Q = 50 # 设置初始订货量
Qmin = 0 # 设置订货量Q的下界
Qmax = 1000 # 设置订货量Q的上界
rhoQ = 1 # 设置订货量Q的变化方向
rhoS = 1 # 设置再订货点S的变化方向
B = [] # 设置订单数量的集合
P = set() # 设置库存策略及其相应目标函数值的集合

# 步骤2：进行模拟
for t in range(len(D)):
    # 更新期末库存
    It = s + Q - D[t]
    # 计算缺货量
    shortage = max(-It, 0) if D[t] > 0 else 0
    # 更新目标函数值
    f_prime = 0.5 * (s + It) + 2 * shortage + 0.5 * max(0, It)
    # 如果目标函数值更优，则更新最优解
    if f_prime < f:
        s_best = s
        S_best = S
        f = f_prime
    # 如果目标函数值相同，则选择再订货点更小的解
    elif f_prime == f and S < S_best:
        s_best = s
        S_best = S
    # 如果当前再订货点和订货量没有被搜索过，则加入集合B
    if (s,Q) not in B:
        B.add((s,Q))
    # 将当前库存策略及其相应目标函数值加入集合P
    P.add((s, S, f_prime))
    # 如果当前订货量为0或达到上下限，则改变搜索方向
    if Q == 0:
        rhoQ = 1
    elif Q == smax - s or s - Q == smin:
        rhoQ = -1
    # 根据搜索方向改变订货量Q
    Q += rhoQ * 10
    # 如果订货量超出上下限，则将搜索方向反转
    if Q > smax - s:
        Q = smax - s
        rhoQ = -1
    elif Q < 0:
        Q = 0
        rhoQ = 1
    # 如果当前再订货点为0或达到上下限，则改变搜索方向
    if S == smin:
        rhoS = 1
    elif S == smax:
        rhoS = -1
    # 根据搜索方向改变再订货点S
    S += rhoS
    
    # 根据订货量Q的变化方向和Q的最小值更新Q
    Qrho = min(I0 - s, smax - s)
    if rhoQ == 1:
        Q = Q + Qrho
    else:
        Q = Q - Qrho

    # 检查(Q-, Q+)区间内是否有Q
    if Q < Qmin or Q > Qmax:
        rhoQ = -rhoQ
        Qrho = min(I0 - s, smax - s)
        if rhoQ == 1:
            Q = Q + Qrho
        else:
            Q = Q - Qrho
        if Q < Qmin or Q > Qmax:
            print("Optimal solution found: (s, S) = ({}, {}), Q = {}".format(s, S, Q))
            break

    # 检查B中是否存在Q，若不存在则加入B中
    if Q not in B:
        B.append(Q)

    # 移动到新的(s', S')并设置ρs
    if rhoQ == 1:
        s_new = s - Qrho
        S_new = S
        rhoS = 1
    else:
        s_new = s
        S_new = S + Qrho
        rhoS = -1

    # 执行模拟
    I = np.zeros(N+1) # 初始化库存为0
    for t in range(1, N+1):
        I[t] = I[t-1] + D[t] - Q
        if I[t] < 0:
            shortage += abs(I[t])
            I[t] = 0
        if t in reorder_points:
            Q = Q + rhoS * Qrho
            Q = max(Qmin, min(Qmax, Q))
