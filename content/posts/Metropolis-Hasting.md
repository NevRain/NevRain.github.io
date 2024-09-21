+++
title = 'Metropolis Hasting采样'
date = 2024-09-20T13:36:27+08:00
draft = false
+++

### 1 问题描述

生成任意概率密度函数$p(x)$的近似采样

### 2 原理

先假定目标概率函数是离散的，$x$的每一个取值称为一个状态，所有$p(x)$的取值可以构成概率分布向量。任意给定初始概率分布向量$\pi^0$，向量中的每一个下$\pi_i^0$表示当前处于状态$i$的概率。任意给定状态转移矩阵$P$，矩阵元素$P_{ij}$表示从状态$i$转移到状态$j$的概率。该状态转移矩阵定义了一条马尔可夫链$$\pi^0\rightarrow \pi^1 \rightarrow \pi^2 \rightarrow\cdots\rightarrow\pi^n$$

$\pi^i$到$\pi^j$的状态转移定义为如下过程$$\pi^j = \pi^i P$$

我们对每一个概率分布向量进行采样就可以得到样本集$x^0,x^1,\cdots,x^n$。

由**马尔可夫链收敛定理**得，从由$\pi^0$定义的初始状态经过$P$的无限多次状态转移后,概率分布向量将收敛到定值$\pi$，称为马尔可夫链的稳定分布。因此，从某一个足够大的样本开始，其后的所有样本所组成的样本集可以近似为对稳定分布$\pi$的采样。

问题被归结，找到合适的状态转移矩阵$P$，使得经过足够多次的状态转移，最终得到的概率分布向量恰好为目标概率密度函数即可。

由**马尔科夫链细致平稳条件**得，我们要找的$P$需要满足对任意状态$i,j$，有以下式子成立$$ \pi(i)P_{ij}=\pi(j)P_{ji} $$

Metropolis算法的思想是，如果任意确定状态转移矩阵$Q$，则通常不满足细致平稳条件，即$$\pi(i)Q_{ij} \neq \pi(j)Q_{ji}$$

但是可以通过添加合适的$\alpha$项，就能满足条件，即$$\pi(i)Q_{ij}\alpha(i,j) = \pi(j)Q_{ji}\alpha(j,i)$$

其中$$\alpha(i,j)=p(j)q(j,i)$$


由此可得Metropolis采样算法：
1. 预先确定任意的状态转移矩阵$Q$，其第$i$行$j$列的元素记作$q(i,j)$
2. 设当前状态为$i$，$Q$的第$i$行为下一状态的概率分布向量$\pi^j$（其所有元素和为1）
3. 对$\pi^j$进行采样，得到下一个状态$j$
4. 设$\alpha(i,j)=\pi(j)q(j,i)$为转移接受概率，即以$\alpha(i,j)$的概率接受这次转移，以$1-\alpha(i,j)$的概率拒绝这次转移
5. 拒绝这次转移后将回到第三步，直到接受转移。

Metropolis-Hasting采样算法是Metropolis采样算法的改进，它将接受概率进行放大，变成$$\pi(i)Q_{ij}\frac{\alpha(i,j)}{\alpha(j,i)} = \pi(j)Q_{ji}$$或者$$\pi(i)Q_{ij} = \pi(j)Q_{ji}\frac{\alpha(j,i)}{\alpha(i,j)}$$

只需要在第4步中，将接受概率修改为$$\alpha'(i,j)=min\{ \frac{p(j)q(j,i)}{p(i)q(i,j)},1 \}$$

进一步的，如果取用的$Q$为对称矩阵，则接受概率修改为$$\alpha'(i,j)=min\{ \frac{p(j)}{p(i)},1 \}$$

Metropolis-Hasting采样算法的精妙之处不仅在于其放大了接受概率从而减少被拒绝的次数，而且其最终的接受概率与状态转移矩阵无关，因此可以安全地将其推广到目标概率函数是**连续**的情况（连续的情况下，状态转移矩阵是无穷维）。


### 3 算法
```pseudocode
def MetropolisHastingSampler(TargetFunc f, int sampleNum):
    skip = 1000 # 任意较大的数，用于收敛到平稳态
    xi = 0
    samples = []

    while(len(samples) < sampleNum):
        while(True):
            xj = random() # 相当于矩阵Q的每一行都是均匀的概率分布向量
            pj = f(xj)
            pi = f(xi)
            a = min(pj/pi, 1)
            r = random()
            if(r <= a):
                xi = xj
                if(skip > 0):
                    skip = skip - 1
                else:
                    samples.append(xj)
                break
    
    return samples
```


### A 定理

#### 马尔可夫链收敛定理

如果一个非周期马尔科夫链具有概率转移矩阵$P$，且任意两状态联通，则$\lim_{n\to\infty}P_{ij}^{n}$存在且与初始概率分布无关，记为$\lim_{n\to\infty}P_{ij}^n=\pi_j$，有：

- $\lim_{n\to\infty}P^n=\begin{bmatrix}\pi_1&\pi_2&\dots&\pi_j\\\pi_1&\pi_2&\dots&\pi_j\\\pi_1&\pi_2&\dots&\pi_j\end{bmatrix}$
- $\pi_j=\sum\pi(i)P_{ij}$
- $\pi$是方程$\pi=\pi P$的唯一非负解

其中$\lim_{n\to\infty}P^{n}$的每一行被称为马尔可夫链的稳定分布。

#### 细致平稳条件

如果非周期马尔科夫链的转移矩阵$P$和分布$\pi(x)$对所有的$i,j$满足$\pi(i)P_{ij}=\pi(j)P_{ji}$，则由$P$确定的马尔科夫链的平衡分布就是$\pi(x)$。