

# MIMO推导
## A .Signal and Channel Model
$$\bf{y_j}=\sqrt{p_{u}}\sum_{i=1}^L\bf{H}_{ji}
x_{l}+n_{j} \tag{1}$$
>$L$ hexagonal cells  $K$ Single users  $M$antenna
>在这里xl服从复高斯分布是在一个cell内的传输信号
>$p_u$是这个cell内所有用户发射功率的平均值（开根号为幅值）
>$n_j$为接受处的噪音符合复高斯分布
>$H_{jl}$则表征了第l个cell到第j个接受机处的信道

接下来集中探索信道 $H_{jl}$的形式
$$\bold{H}_{jl}\triangleq[\bf{h}_{jl1},…,h_{jlk},..h_{jlK}] \tag{2}$$
>这个式子本质上在描述一个cell里面的单独一个用户与接受段的信道

并且我们对信道加入衰落以更好描述其性质
$$\bold{h}_{jlk}=\bold{\underline{h}}_{jkl}\beta_{jkl}^\frac{1}{2} \tag3$$
>在这个式子里$\beta$是大尺度衰落的系数
>$\underline{h}_{jkl}$是服从复高斯分布的
## B.Channel Estimation
在本节中我们会通过发送导频信号来估计信道信息，首先我们来定义uplink pilot sequence
$$\bold{s}_{lk}=\sqrt{\rho_{lk}}\underline{s}_{lk}\tag4$$
> $lk$指的是第l个cell里的第k个用户
> 导频矩阵是$\tau*1$的
> 发送信号时是乘上了幅值$\sqrt{\rho_{lk}}$的

但我们仍然难以把握导频信号的规律，那么不妨假设其具有正交性以方便计算
$$\bold{\underline{s}}_{lk}^H\bold{\underline{s}}_{lk}=1 \text{  and } \bold{\underline{s}}_{lk_1}^H\bold{\underline{s}}_{lk_2}=1 \quad \forall k1\neq k2 \tag{5}$$

>成立的条件是$\tau\geqq K$ 也就是导频信号的长度大于用户数

由于我们想去探寻最差情况下MIMO表现，我们不妨加入导频污染（用户之间的相互影响）
$$\bold{\underline{s}}_{ik}=\bold{\underline{s}}_{jk} \tag{6}$$
>也就是说不同cell内的第k个用户是发送了同样的（相互影响的导频）

并且我们可以列出此时第j个接受机里接受到的信号
$$ \bold{Y}_j=\sum_{l=1}^L \sum_{n=1}^K  \sqrt{\rho_{ln}}\bold{h}_{jln}\bold{\underline{s}}_{ln}^H+\bold{N}_j \tag{7}$$ 

>$N_j$是符合复高斯分布的接受端噪声
>其实就是每个cell中的每个用户通过和j相关的信道发送导频信号并且在接收端叠加高斯噪声
### B-1 LS信道估计法及其详细推导
LS信道估计法的表达式已经被定义
$$\bold{\widehat{h}}_{jjk}^{LS}\triangleq\frac{1}{\sqrt{\rho_{jk}}}\bold{Y}_j\bold{\underline{s}}_{jk} \tag8$$
现在将5-7带入8注意到导频的正交性，以及导频污染带来的复用

带入式子7

$$ \bold{\widehat{h}}_{jjk}^{LS}=\frac{1}{\sqrt{\rho_{jk}}}(\sum_{l=1}^L \sum_{i=1}^K  \sqrt{\rho_{ln}}\bold{h}_{jln}\bold{\underline{s}}_{ln}^H+\bold{N}_j)\bold{\underline{s}}_{jk} \tag{8.1} $$

$$ \bold{\widehat{h}}_{jjk}^{LS}=\sum_{l=1}^L \sum_{n=1}^K \frac{\sqrt{\rho_{ln}}}{\sqrt{\rho_{jk}}} \bold{h}_{jln}\bold{\underline{s}}_{ln}^H\bold{\underline{s}}_{jk}  +\frac{1}{\sqrt{\rho_{jk}}} \bold{N}_j\bold{\underline{s}}_{jk} \tag{8.2}$$
不难发现，此时可以利用导频信号的正交性与式子（6）进行化简
$$\left\{\begin{align}
\bold{\underline{s}}_{jn}^H\bold{\underline{s}}_{jk}=1&&n=k\\\bold{\underline{s}}_{jn}^H\bold{\underline{s}}_{jk}=0&&n \neq k\\
 \end{align}\right. \tag{8.3}$$
 $$\bold{\underline{s}}_{ik}=\bold{\underline{s}}_{jk} \tag{6}$$
 因此观察B.2发现第二个积分的范围内仅有$n=k$保留下来了($l～[1，L]$由于导频复用保留)，而且$\rho_{ln}$功率的下标改为$\rho_{lk}$
 
$$ \bold{\widehat{h}}_{jjk}^{LS}=\sum_{l=1}^L \frac{\sqrt{\rho_{lk}}}{\sqrt{\rho_{jk}}} \bold{h}_{jlk} +\frac{1}{\sqrt{\rho_{jk}}} \bold{N}_j\bold{\underline{s}}_{jk} \tag{8.4}$$
进一步不妨将$l=j$的部分单独提出，以得到$h_{jjk}$ 于是我们得到了式子9
$$ \bold{\widehat{h}}_{jjk}^{LS}=\bold{h}_{jjk}+\sum_{l\neq j}^L \frac{\sqrt{\rho_{lk}}}{\sqrt{\rho_{jk}}} \bold{h}_{jlk} +\frac{1}{\sqrt{\rho_{jk}}} \bold{N}_j\bold{\underline{s}}_{jk} \tag{9}$$
### B-2 MMSE信道估计法及其详细推导
$$\bold{\widehat{h}}_{MMSE}^{LS}\triangleq\mathbb{E}[\bold{\widehat{h}}_{jjk}(\bold{\widehat{h}}_{jjk}^{LS})^H](\mathbb{E}[\bold{\widehat{h}}_{jjk}^{LS}(\bold{\widehat{h}}_{jjk}^{LS})^H])^{-1}\bold{\widehat{h}}_{jjk}^{LS}  \tag{10}$$
接下来的核心内容便是利用无关性化简上述式子
#### B-2-1化简$\mathbb{E}[\bold{h}_{jjk}(\bold{\widehat{h}}_{jjk}^{LS})^H]$

$$ \bold{h}_{jjk}(\bold{\widehat{h}}_{jjk}^{LS})^H= \bold{h}_{jjk}(\bold{h}_{jjk}+\sum_{l\neq j}^L \frac{\sqrt{\rho_{lk}}}{\sqrt{\rho_{jk}}} \bold{h}_{jlk} +\frac{1}{\sqrt{\rho_{jk}}} \bold{N}_j\bold{\underline{s}}_{jk})^H\tag{10.1}$$
$$ \bold{h}_{jjk}(\bold{\widehat{h}}_{jjk}^{LS})^H= \bold{h}_{jjk}\bold{h}_{jjk}^H+\bold{h}_{jjk}\sum_{l\neq j}^L \frac{\sqrt{\rho_{lk}}}{\sqrt{\rho_{jk}}} \bold{h}_{jlk}^H +\frac{1}{\sqrt{\rho_{jk}}}\bold{h}_{jjk} \bold{\underline{s}}_{jk}^H\bold{N}_j^H\tag{10.2}$$
此时两边同时求期望
$$\mathbb{E}[\bold{h}_{jjk}(\bold{\widehat{h}}_{jjk}^{LS})^H]=\mathbb{E}[\bold{h}_{jjk}\bold{h}_{jjk}^H+\bold{h}_{jjk}\sum_{l\neq j}^L \frac{\sqrt{\rho_{lk}}}{\sqrt{\rho_{jk}}} \bold{h}_{jlk}^H +\frac{1}{\sqrt{\rho_{jk}}}\bold{h}_{jjk} \bold{\underline{s}}_{jk}^H\bold{N}_j^H\tag{10.3}]$$
$$ =\mathbb{E}[\bold{h}_{jjk}\bold{h}_{jjk}^H]+\mathbb{E}[\bold{h}_{jjk}\sum_{l\neq j}^L \frac{\sqrt{\rho_{lk}}}{\sqrt{\rho_{jk}}} \bold{h}_{jlk}^H]+\frac{1}{\sqrt{\rho_{jk}}}\mathbb{E}[\bold{h}_{jjk} \bold{\underline{s}}_{jk}^H\bold{N}_j^H] \tag{10.4}$$
由于信道衰落前的信号$\bold{h}_{jjk}$服从复高斯分布，因此其期望为0，方差为1矩阵，所以
$$\mathbb{D}[\bold{h}_{jjk}]=((\beta_{jjk})^{1/2})^2\mathbb{D}[\bold{\underline{h}}_{jjk}]$$
$$\mathbb{D}[\bold{h}_{jjk}]=\mathbb{cov}[\bold{h}_{jjk},\bold{h}_{jjk}]=\bold {\beta_{jjk}}I$$
$$\mathbb{cov}[\bold{h}_{jjk},\bold{h}_{jjk}]=\mathbb{E}[\bold{h}_{jjk}\bold{h}_{jjk}^H]+\mathbb{E}[\bold{h}_{jjk}]\mathbb {E}[\bold{h}_{jjk}^H]=\bold{\beta_{jjk}}I$$
$$\mathbb{E}[\bold{h}_{jjk}]=\mathbb {E}[\bold{h}_{jjk}^H]=\bold{0}$$
$$\mathbb{E}[\bold{h}_{jjk}\bold{h}_{jjk}^H]=\bold{\beta_{jjk}}I$$
而对于两个分别复从两个独立复高斯分布的随机变量
$$l\neq j$$
$$\mathbb{cov}[\bold{h}_{jjk},\bold{h}_{jlk}]=\mathbb{E}[\bold{h}_{jjk}\bold{h}_{jlk}^H]+\mathbb{E}[\bold{h}_{jjk}]\mathbb {E}[\bold{h}_{jlk}^H]=\bold{0}$$
$$\mathbb{E}[\bold{h}_{jjk}]=0$$
$$\mathbb{E}[\bold{h}_{jjk}\bold{h}_{jlk}^H]=0$$
由上述式子
>**两个独立高斯分布的内积的期望为0 【1】**
>**一个高斯分布与其自己的内积的期望等于其方差 [2]**

不再做证明，同理由于$\bold{\underline{s}}_{jk}^H\bold{N}_j^H$服从复高斯分布所以
$$\mathbb{E}[\bold{h}_{jjk}(\bold{\widehat{h}}_{jjk}^{LS})^H]=\bold {\beta_{jjk}}I$$
$$\mathbb{E}[\sum_{l\neq j}^L \frac{\sqrt{\rho_{lk}}}{\sqrt{\rho_{jk}}} \bold{h}_{jln}^H]=0$$
$$\mathbb{E}[\bold{h}_{jjk} \bold{\underline{s}}_{jk}^H\bold{N}_j^H]=0$$
$$\mathbb{E}[\bold{h}_{jjk}(\bold{\widehat{h}}_{jjk}^{LS})^H]=\bold {\beta_{jjk}}I \tag{10.5}$$
#### B-2-2 化简$\mathbb{E}[\bold{\widehat{h}}_{jjk}^{LS}(\bold{\widehat{h}}_{jjk}^{LS})^H]$

$$ \bold{\widehat{h}}_{jjk}^{LS}(\bold{\widehat{h}}_{jjk}^{LS})^H=( 
\bold{h}_{jjk}+\sum_{l\neq j}^L \frac{\sqrt{\rho_{lk}}}{\sqrt{\rho_{jk}}} \bold{h}_{jln} +\frac{1}{\sqrt{\rho_{jk}}} \bold{N}_j\bold{\underline{s}}_{jk})(\bold{h}_{jjk}^H+\sum_{l\neq j}^L \frac{\sqrt{\rho_{lk}}}{\sqrt{\rho_{jk}}} \bold{h}_{jln}^H +\frac{1}{\sqrt{\rho_{jk}}} \bold{\underline{s}}_{jk}^H\bold{N}_j^H)$$

$$ =\bold{h}_{jjk}\bold{h}_{jjk}^H+\bold{h}_{jjk}\sum_{l\neq j}^L \frac{\sqrt{\rho_{lk}}}{\sqrt{\rho_{jk}}} \bold{h}_{jln}^H +\frac{1}{\sqrt{\rho_{jk}}} \bold{h}_{jjk}\bold{\underline{s}}_{jk}^H\bold{N}_j^H+ 
\sum_{l\neq j}^L \frac{\sqrt{\rho_{lk}}}{\sqrt{\rho_{jk}}} \bold{h}_{jln} \bold{h}_{jjk}^H+\sum_{l_1\neq j}^L \frac{\sqrt{\rho_{l_1k}}}{\sqrt{\rho_{jk}}} \bold{h}_{jl_1n} \sum_{l_2\neq j}^L \frac{\sqrt{\rho_{l2k}}}{\sqrt{\rho_{jk}}} \bold{h}_{jl_2n}^H $$$$+\frac{1}{\sqrt{\rho_{jk}}} \sum_{l\neq j}^L \frac{\sqrt{\rho_{lk}}}{\sqrt{\rho_{jk}}} \bold{h}_{jln} \bold{\underline{s}}_{jk}^H\bold{N}_j^H+\frac{1}{\sqrt{\rho_{jk}}} \bold{N}_j\bold{\underline{s}}_{jk}\bold{h}_{jjk}^H+\frac{1}{\sqrt{\rho_{jk}}} \bold{N}_j\bold{\underline{s}}_{jk}\sum_{l\neq j}^L \frac{\sqrt{\rho_{lk}}}{\sqrt{\rho_{jk}}} \bold{h}_{jln}^H +\frac{1}{\sqrt{\rho_{jk}}} \bold{N}_j\bold{\underline{s}}_{jk}\frac{1}{\sqrt{\rho_{jk}}} \bold{\underline{s}}_{jk}^H\bold{N}_j^H$$
由于展开式子以及结论 **【1】【2】**其中第一项来源于上式第一项，第二项来自于上式第五项**【3】**，第三项来自于上式最后一项
$$\mathbb{E}[\bold{\widehat{h}}_{jjk}^{LS}(\bold{\widehat{h}}_{jjk}^{LS})^H]=（\beta_{jjk} +\sum_{l\neq j}^L\frac{\rho_{lk}}{\rho_{jk}}\beta_{jlk}+\frac{1}{\rho_{jk}})I \tag{10.6}$$
>**[3]** 当且仅当$l_1=l_2$时为一个高斯分布自己的内积的期望，此时式子保留结果为$\beta I$，其余的交叉项均为0，类似于正交性
#### B-2-3综合上面的式子
将10.6 10.5带入10
$$\mathbb{E}[\bold{\widehat{h}}_{jjk}(\bold{\widehat{h}}_{jjk}^{LS})^H](\mathbb{E}[\bold{\widehat{h}}_{jjk}^{LS}(\bold{\widehat{h}}_{jjk}^{LS})^H])^{-1}\bold{\widehat{h}}_{jjk}^{LS}=\bold {\beta_{jjk}}I[(\beta +\sum_{l\neq j}^L\frac{\rho_{lk}}{\rho_{jk}}+\frac{1}{\rho_{jk}})I]^{-1} \bold{\widehat{h}}_{jjk}^{LS}$$
$$=\bold {\beta_{jjk}}II^{-1}\frac{1}{\beta_{jjk} +\sum_{l\neq j}^L\frac{\rho_{lk}}{\rho_{jk}}\beta_{jlk}+\frac{1}{\rho_{jk}} }\bold{\widehat{h}}_{jjk}^{LS}$$
$$=\frac{\bold {\beta_{jjk}}\rho_{jk}}{\sum_{l=1}^L{\rho_{lk}}\beta_{jlk}+1}\bold{\widehat{h}}_{jjk}^{LS}\tag{11}$$
综上所述，目标达成






