import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from astropy.utils.masked.function_helpers import zeros_like


# 更新图示的函数
def update_population_plot(X, fitness, gbest, title):
    plt.clf()  # 清除旧图形
    plt.scatter(X[:, 0], X[:, 1], c=fitness, cmap='coolwarm', s=50, alpha=0.6, vmin=min(fitness), vmax=max(fitness))
    plt.colorbar(label="Fitness")
    plt.scatter(gbest[0], gbest[1], c='green', s=200, marker='X', label='Best Individual')
    plt.title(title)
    plt.xlabel('X1')
    plt.ylabel('X2')
    plt.legend()
    plt.draw()  # 强制刷新图形
    plt.pause(0.01)  # 强制刷新图形并暂停

def tent_map(popsize, dim, r=1.99, lb=-1, ub=1):
    """
    使用 Tent Mapping 生成种群位置
    r: Tent映射的常数，默认1.99
    lb: 下界
    ub: 上界
    """
    # 初始化种群矩阵
    X = np.zeros((popsize, dim))

    # 随机初始化一个种群位置
    x = np.random.rand(popsize, dim)  # x 取值在 [0, 1] 之间

    for i in range(popsize):
        for j in range(dim):
            # 使用 Tent Map 更新每个维度的位置
            X[i, j] = lb + (ub - lb) * (r * abs(x[i, j] - 0.5))  # 映射到 [lb, ub]

    return X

# 目标函数，f(x, y) = x^4 + y^4
def example_fhd(x):
    return np.sum(x ** 2)  # x^4 + y^4

def TBESO(fhd, pos,popsize, maxgen, lb, ub, dim,*args):
    # Initial settings
    flag=",initial population";
    # Initialize the population (using the provided initial positions)
    X = np.copy(pos)
    fitness = np.array([fhd(x) for x in X])
    vec_flag = [1, -1]
    Threshold = 0.25
    Threshold2 = 0.6
    C1 = 0.5
    C2 = 0.05
    C3 = 2

    # Best fitness and position initialization
    gbest = np.argmin(fitness) # 最佳个体的index
    bestf=np.min(fitness) #最佳fitness
    bestX=X[gbest, :] #最佳个体的位置，变量x，y值

    # Divide the swarm into two equal groups: males and females
    Nm = round(popsize / 2)
    Nf = popsize - Nm
    Xm = X[:Nm, :]
    Xf = X[Nm:, :]
    fitness_m = fitness[:Nm]
    fitness_f = fitness[Nm:]

    # Find the best male and female
    fitnessBest_m = np.min(fitness_m)
    gbest1 = np.argmin(fitness_m)
    Xbest_m = Xm[gbest1, :]
    fitnessBest_f = np.min(fitness_f)
    gbest2 = np.argmin(fitness_f)
    Xbest_f = Xf[gbest2, :]

    # 主迭代过程
    curve = []
    draw_plot(0,X,fitness,bestX,f"TBESO:Generation {0}"+flag)
    save_csv(X, fitness, 0, "./data/Initial_population.csv")

    for t in range(maxgen):

        # Calculate Temperature and Food Quality
        Temp = np.exp(-t / maxgen)  # Higher t lower Temp
        Q = C1 * np.exp((t - maxgen) / maxgen)  # Higher t higher Q
        if Q>1: Q=1

        ## the variable for BDS
        gmworst = np.argmax(fitness_m)
        Xworst_m = Xm[gmworst, :]
        gfworst = np.argmax(fitness_f)
        Xworst_f = Xf[gfworst, :]

        Xnew = np.zeros_like(X)
        Xnewm = Xnew[:Nm, :]
        Xnewf = Xnew[Nm:, :]
        fitness_new=np.zeros_like(fitness)

        if Q < Threshold:
            flag = ",Exploration";
            # Exploration Phase (no Food)
            Xnewmfit=np.zeros_like(Xm)
            Xnewffit=np.zeros_like(Xf)
            save_csv(X, fitness, 0, f"./data/BDS/before/before_BDS{t}.csv")
            for i in range(Nm):
                rand_leader_index = np.random.randint(0, Nm)
                X_randm = Xm[rand_leader_index, :]
                Flag = np.random.choice([-1, 1])
                Am = np.exp(-fitness_m[rand_leader_index] / (fitness_m[i] + np.finfo(float).eps))
                for j in range(dim):
                    Xnewm1 = X_randm[j] + Flag * C2 * Am * ((ub[j] - lb[j]) * np.random.rand() + lb[j])
                    Xnewm2= Xm[i, j] + np.random.rand() * (Xbest_m[j] - Xm[i, j]) - np.random.rand() * (
                            Xworst_m[j] - Xm[i, j])
                    fit1=fhd(Xnewm1)
                    fit2=fhd(Xnewm2)
                    if fit1<fit2:
                        Xnewmfit[i,j]=Xnewm1
                    else:
                        Xnewmfit[i,j]=Xnewm2
                # 选择产生的最佳个体
                Xnewm[i, :] = Xnewmfit[i,:]
            for i in range(Nf):
                rand_leader_index = np.random.randint(0, Nf)
                X_randf = Xf[rand_leader_index, :]
                Flag = np.random.choice([-1, 1])
                Af = np.exp(-fitness_f[rand_leader_index] / (fitness_f[i] + np.finfo(float).eps))
                for j in range(dim):
                    Xnewf1 = X_randf[j] + Flag * C2 * Af * ((ub[j] - lb[j]) * np.random.rand() + lb[j])
                    Xnewf2 = Xf[i, j] + np.random.rand() * (Xbest_f[j] - Xf[i, j]) - np.random.rand() * (
                        Xworst_f[j] - Xf[i, j])
                    fit3 = fhd(Xnewf1)
                    fit4 = fhd(Xnewf2)
                    if fit3 < fit4:
                        Xnewffit[i,j] = Xnewf1
                    else:
                        Xnewffit[i,j] = Xnewf2
                # 选择产生的最佳个体
                Xnewf[i, :] = Xnewffit[i,:]
            fitness_new = np.array([fhd(x) for x in Xnew])
            save_csv(Xnew, fitness_new, 0, f"./data/BDS/after/after_BDS{t}.csv")
        else:
            # Exploitation Phase (Food Exists)
            if Temp > Threshold2:
                flag = ",Exploitation,hot";
                # Hot Mode
                for i in range(Nm):
                    Flag = np.random.choice([-1, 1])
                    Xnewm[i, :] = bestX + C3 * Flag * Temp * np.random.rand(dim) * (bestX - Xm[i, :])
                    # for j in range(dim):
                    #     Xnewm[i, j] = bestX[j] + C3 * Flag * Temp * np.random.rand() * (bestX[j] - Xm[i, j])

                for i in range(Nf):
                    Flag = np.random.choice([-1, 1])
                    Xnewf[i, :] = bestX + C3 * Flag * Temp * np.random.rand(dim) * (bestX - Xf[i, :])
                    # for j in range(dim):
                    #     Xnewf[i, j] = bestX[j] + Flag * C3 * Temp * np.random.rand() * (bestX[j] - Xf[i, j])
            else:
                flag = ",Exploitation,cold";
                # Cold Mode
                if np.random.rand() > 0.6:
                    flag = ",Exploitation,cold,fighting"
                    # Fight Mode
                    for i in range(Nm):
                        FM = np.exp(-fitness_f.min() / (fitness_m[i] + 1e-8))
                        Xnewm[i, :] = Xm[i, :] + C3 * FM * np.random.rand(dim) * (Q * Xbest_f - Xm[i, :])
                        # for j in range(dim):
                        #     FM = np.exp(-fitnessBest_f / (fitness_m[i] + np.finfo(float).eps))
                        #     Xnewm[i, j] = Xm[i, j] + C3 * FM * np.random.rand() * (Q * Xbest_f[j] - Xm[i, j])

                    for i in range(Nf):
                        FF = np.exp(-fitness_m.min() / (fitness_f[i] + 1e-8))
                        Xnewf[i, :] = Xf[i, :] + C3 * FF * np.random.rand(dim) * (Q * Xbest_m - Xf[i, :])
                        # for j in range(dim):
                        #     FF = np.exp(-fitnessBest_m / (fitness_f[i] + np.finfo(float).eps))
                        #     Xnewf[i, j] = Xf[i, j] + C3 * FF * np.random.rand() * (Q * Xbest_m[j] - Xf[i, j])
                else:
                    flag = ",Exploitation,cold,mating"
                    fitness_new = np.array([fhd(x) for x in Xnew])
                    save_csv(X, fitness, 0, f"./data/EOBL/before/before_MATING{t}.csv")
                    # Mating Mode
                    for i in range(Nm):
                        for j in range(dim):
                            Mm = np.exp(-fitness_f[i] / (fitness_m[i] + 1e-8))
                            Xnewm[i, j] = Xm[i, j] + C3 * np.random.rand() * Mm * (Q * Xf[i, j] - Xm[i, j])

                    for i in range(Nf):
                        for j in range(dim):
                            Mf = np.exp(-fitness_m[i] / (fitness_f[i] + 1e-8))
                            Xnewf[i, j] = Xf[i, j] + C3 * np.random.rand() * Mf * (Q * Xm[i, j] - Xf[i, j])
                    ## produce eggs
                    egg = np.random.choice([-1, 1])
                    if egg==1:
                        gworst=np.argmax(fitness_m)
                        Xnewm[gworst,:]=lb+np.random.randn(1,dim)*(ub-lb)
                        gworst = np.argmax(fitness_f)
                        Xnewf[gworst, :] = lb + np.random.randn(1, dim) * (ub - lb)
                    fitness_new = np.array([fhd(x) for x in Xnew])
                    save_csv(Xnew, fitness_new, 0, f"./data/EOBL/after/after_MATING{t}.csv")

        # ---
        ym = np.array([fhd(x) for x in Xnewm])
        yf = np.array([fhd(x) for x in Xnewf])
        # Elite reverse learning
        flag += ",EOBL"
        fitness_new = np.array([fhd(x) for x in Xnew])
        save_csv(Xnew, fitness_new, 0, f"./data/EOBL/before/before_EOBL{t}.csv")
        I1 = np.argsort(ym)

        a = round(0.1 * Nm)
        m_elite = Xnewm[I1[:a], :]

        for j in range(Nm):
            obl_Xnewm=np.zeros_like(Xnewm)
            for i in range(dim):
                # 反向解=k*(a+b)-x
                obl_Xnewm[j,i] = np.random.rand() * (np.max(m_elite[:, i]) + np.min(m_elite[:, i])) - Xnewm[j,i]
                ## 超出边界，就调整回来
                if obl_Xnewm[j,i] > ub[i] or obl_Xnewm[j,i] < lb[i]:
                    # rand*(a-b)+b
                    obl_Xnewm[j,i] = np.random.rand() * (np.max(m_elite[:, i]) - np.min(m_elite[:, i])) + np.min(m_elite[:, i])
            fit1 = fhd(obl_Xnewm[j,:])
            # 如果这个新的解比原先好
            if fit1 < ym[j]:
                Xnewm[j, :] = obl_Xnewm[j,:]

        I2 = np.argsort(yf)
        a = round(0.1 * Nf)
        f_elite = Xnewf[I2[:a], :]
        for j in range(Nf):
            obl_Xnewf=np.zeros_like(Xnewf)
            for i in range(dim):
                obl_Xnewf[j,i] = np.random.rand() * (np.max(f_elite[:, i]) + np.min(f_elite[:, i])) - Xnewf[j, i]
                if obl_Xnewf[j,i] > ub[i] or obl_Xnewf[j,i] < lb[i]:
                    obl_Xnewf[j,i] = np.random.rand() * (np.max(f_elite[:, i]) - np.min(f_elite[:, i])) + np.min(f_elite[:, i])
            fit1 = fhd(obl_Xnewf[j,:])
            if fit1 < yf[j]:
                Xnewf[j, :] = obl_Xnewf[j,:]

        fitness_new=np.array([fhd(x) for x in Xnew])
        save_csv(Xnew, fitness_new, 0, f"./data/EOBL/after/after_EOBL{t}.csv")

        # 更新雄性
        for j in range(Nm):
            # 将 Xnewm[j,:] 的元素限制在区间 [lb, ub] 内
            Xnewm[j, :] = np.clip(Xnewm[j, :], lb, ub)

            # 计算当前个体的适应度
            y = fhd(Xnewm[j, :])

            # 如果新的适应度更好，则更新适应度和位置
            if y < fitness_m[j]:
                fitness_m[j] = y
                Xm[j, :] = Xnewm[j, :]

        gbest1 = np.argmin(fitness_m)
        Ybest1 = fitness_m[gbest1]

        # 更新雌性
        for j in range(Nf):
            # 将 Xnewf[j,:] 的元素限制在区间 [lb, ub] 内
            Xnewf[j, :] = np.clip(Xnewf[j, :], lb, ub)

            # 计算当前个体的适应度
            y = fhd(Xnewf[j, :])

            # 如果新的适应度更好，则更新适应度和位置
            if y < fitness_f[j]:
                fitness_f[j] = y
                Xf[j, :] = Xnewf[j, :]

        gbest2 = np.argmin(fitness_f)
        Ybest2 = fitness_f[gbest2]

        # return the best solution
        if Ybest1 < fitnessBest_m:
            Xbest_m = Xm[gbest1, :]
            fitnessBest_m = Ybest1
        if Ybest2 < fitnessBest_f:
            Xbest_f = Xf[gbest2, :]
            fitnessBest_f = Ybest2
        if fitnessBest_m < fitnessBest_f:
            bestX = Xbest_m
            bestf = fitnessBest_m
        else:
            bestX = Xbest_f
            bestf = fitnessBest_f

        draw_plot(t, X, fitness, bestX, flag)
        curve.append(bestf)

    return curve, bestX, bestf

def scale_fitness_to_loss(fitness):
    """Normalize fitness to the range [0, 1]"""
    #return fitness
    loss=np.zeros_like(fitness)
    count = 0
    min=fitness.min()
    flag=0
    if(min<1e-10):flag=1
    for i in fitness:
        if(flag==1):
            i=i*1e+10
        loss[count]=-i
        count+=1
    min=loss.min()
    max=loss.max()
    count=0
    for i in loss:
        loss[count]=(i - min) / (max - min + 1e-8)
        print("i=",i,",loss=[",count,"]=",loss[count])
        count += 1
    return loss

def save_csv(X,fitness,t,url):
    # 创建 DataFrame
    gender_values=[]
    count=0
    for i in fitness:
        gender_values.append(count%2)
        count+=1
    loss=scale_fitness_to_loss(fitness)
    df = pd.DataFrame(X, columns=[f"x{i + 1}" for i in range(dim)])
    df['fitness'] = loss
    df['gender']=gender_values

    # 将每一代数据保存为 CSV
    if url == "":
        url=f"./data/each_iterations/population_gen_{t}.csv"
    df.to_csv(url, index=False)

def draw_plot(t,X,fitness, bestx, flag):
    loss = scale_fitness_to_loss(fitness)
    save_csv(X,fitness,t,"")
    # print(f"TBESO:Generation {t}" + flag)
    # 每10代更新一次图形
    if t==0:
        update_population_plot(X, loss, bestx, f"TBESO:Generation {t}" + flag)
    if t!=0 and t % 5 == 0:
        update_population_plot(X, loss, bestx, f"TBESO:Generation {t }" + flag)
    if t==20:
        update_population_plot(X, loss, bestx, f"TBESO:Generation {t }" + flag)

# 设置参数
popsize = 400 # 种群大小
maxgen = 25  # 最大迭代次数
dim = 2  # 变量维度
lb=-1
ub=1
#tent map 种群生成
r=1.99
#pos = np.random.uniform(lb, ub, (popsize, dim))  # 初始化位置

#pos = tent_map(popsize, dim, r=1.99, lb=lb, ub=ub)  # 使用Tent Map初始化种群

# 保存为 Excel 文件
#df = pd.DataFrame(pos, columns=[f"x_{i+1}" for i in range(dim)])  # 为每一列命名
#output_filename = "./data/tent_map_population.xlsx"
#df.to_excel(output_filename, index=False)

#print(f"数据已保存到: {output_filename}")

# # 读取保存的 Excel 文件
output_filename = "./data/tent_map_population.xlsx"
pos = pd.read_excel(output_filename)

lb = lb * np.ones(dim)  # 下界
ub = ub * np.ones(dim)  # 上界
# 运行TBESO算法并绘制动态结果
curve, bestx, bestf = TBESO(example_fhd, pos, popsize, maxgen, lb, ub, dim)

# 显示最终最优解
print(f"Best fitness: {bestf}")
