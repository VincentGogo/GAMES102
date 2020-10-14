# 作业 1

> 2020/10/14

## 问题

**Input**：已知平面内 n 个点 $P_j(x_j,y_j), j=1,2,\dots,n$。

**Output**: 拟合这些点的函数。

要求：实现不同的拟合方法，并进行比较。输入点集可以进行交互式鼠标指定，或者其他方法生成。

![https://i.bmp.ovh/imgs/2020/10/c10f7ef259809aa0.png](https://i.bmp.ovh/imgs/2020/10/c10f7ef259809aa0.png)

## 一、插值型拟合方法

#### 1. 使用多项式函数（幂基函数的线性组合）

$f(x)=\sum_{i=0}^{n-1}\alpha_i B_i(x)$ 插值 $\{P_j\}$，其中 $B_i(x)=x^i$ 

用矩阵向量表示即：
$$
\begin{equation}
\left[
\begin{matrix}
x_0^n & x_0^{n-1} & x_0^{n-2} & \cdots & x_0^1 &1 \\
x_1^n & x_1^{n-1} & x_1^{n-2} & \cdots & x_1^1 &1 \\
\vdots &\vdots    &\vdots &\ddots &\vdots &\vdots \\
x_n^n & x_n^{n-1} & x_n^{n-2} & \cdots & x_n^1 &1 \\
\end{matrix}
\right]
\left[
\begin{matrix}
\alpha_0\\
\alpha_1\\
\vdots\\
\alpha_n\\
\end{matrix}
\right]=
\left[
\begin{matrix}
y_0\\
y_1\\
\vdots\\
y_n
\end{matrix}
\right]
\end{equation}
$$
![image-20201014171938604](C:\Users\lch\AppData\Roaming\Typora\typora-user-images\image-20201014171938604.png)

#### 2. 使用 Gauss 基函数的线性组合 

$f(x)=b_0 + \sum_{i=1}^{n}b_i g_i(x)$  插值 $\{P_j\}$，其中
$$
g_i(x)=\exp\left(-\frac{(x-x_i)^2}{2\sigma^2}\right)
$$
即对称轴在插值点上，$i=1,\dots,n$，缺省设 $\sigma =1$ 

用矩阵向量表示：
$$
\begin{equation}
\left[
\begin{matrix}
1 & g_1(x_1) & g_2(x_1)  & \cdots & g_n(x_1)  \\
1 & g_1(x_2) & g_2(x_2)  & \cdots & g_n(x_2)  \\
\vdots &\vdots    &\vdots &\ddots &\vdots  \\
1 & g_1(x_n) & g_2(x_n)  & \cdots & g_n(x_n)  \\
\end{matrix}
\right]
\left[
\begin{matrix}
b_0\\
b_1\\
\vdots\\
b_n\\
\end{matrix}
\right]=
\left[
\begin{matrix}
y_1\\
y_2\\
\vdots\\
y_n
\end{matrix}
\right]
\end{equation}
$$
$\sigma = 10.995  图像如下

![](C:\Users\lch\AppData\Roaming\Typora\typora-user-images\image-20201014172243368.png)

$\sigma = 5.062$  图像如下

![](C:\Users\lch\AppData\Roaming\Typora\typora-user-images\image-20201014172154035.png)

$\sigma = 20.027$  图像如下

![](C:\Users\lch\AppData\Roaming\Typora\typora-user-images\image-20201014172433540.png)

$\sigma = 50$  图像如下

![](C:\Users\lch\AppData\Roaming\Typora\typora-user-images\image-20201014172511619.png)



## 二、逼近型拟合方法

#### 最小二乘法

 固定幂基函数的最高次数m (m<n)，使用最小二乘法：$\min E$，其中 $E(x)=\sum_{i=0}^{n}(y_i-f(x_i))^2$ 拟合 $\{P_j\}$。

固定最高次数后多项式函数为m (m<n) ：$f(x) = a_0 + a_1x^1 + a_2x^2\cdots + a_mx^m$ ，则有：
$$
E(x)=\sum_{i=0}^{n}(y_i-f(x_i))^2$=\sum_{i=0}^{n}(y_i - \sum_{j=0}^{m}{a_j}{x_i^j})^2
$$
对 $a_k(k\in\{0,1,\cdots,m\})$ 求导：
$$
\begin{aligned} 
\frac {\partial E}{\partial a_k} &=  \frac{\partial}{\partial a_k}\sum_{i=0}^{n}(y_i-\sum_{j=0}^{m}{a_j}{x_i^j})^2 \\
 &= \sum_{i=0}^{n}\frac{\partial}{\partial a_k}(y_i - \sum_{j=0}^{m}{a_j}{x_i^j})^2 \\ 
 &= \sum_{i=0}^{n}2(y_i-\sum_{j=0}^{m}{a_j}{x^j}) \frac{\partial}{\partial a_k}({y_i - \sum_{j=0}^{m}{a_j}{x_i^j}})\\
 &= -2\sum_{i=0}^{n}{x_i^k}(y_i - \sum_{j=0}^{m}{a_j}{x_i^j})
\end{aligned}
$$
要求得最小 $\min E$ , 对每个 $a_k$ 的偏导数需要取值为0，即
$$
\frac{\partial E}{\partial a_k} = x_0^k(y_0-\sum_{j=0}^{m}{a_j}{x_0^j}) + x_1^k(y_1-\sum_{j=0}^{m}{a_j}{x_1^j}) + \cdots + x_n^k(y^n-\sum_{j=0}^{m}{a_j}{x_n^j}) = 0
$$
写成矩阵形式：
$$
\begin{pmatrix}
\frac{\partial E}{\partial a_0} \\
\frac{\partial E}{\partial a_1} \\
\vdots\\
\frac{\partial E}{\partial a_m} \\
\end{pmatrix}
= 
\begin{pmatrix}
1 & 1 &1 &\cdots & 1\\
x_0 & x_1 & x_2 &\cdots & x_n\\
x_0^2 &x_1^2 &x_2^2 &\cdots &x_n^2\\
\vdots &\vdots &\vdots &\ddots &\vdots\\
x_0^m &x_1^m &x_2^m &\cdots &x_n^m\\
\end{pmatrix}
\begin{pmatrix}
y_0-\sum_{j=0}^{m}{a_j}{x_0^j}\\
y_1-\sum_{j=0}^{m}{a_j}{x_0^j}\\
\vdots\\
y_n-\sum_{j=0}^{m}{a_j}{x_0^j}\\
\end{pmatrix}
= 
\begin{pmatrix}
0\\
0\\
\vdots\\
0
\end{pmatrix}
$$
令
$$
X = 
\begin{pmatrix}
1 &x_0 &x_0^2 &\cdots &x_0^m \\
1 &x_1 &x_1^2 &\cdots &x_1^m \\
1 &x_2 &x_2^2 &\cdots &x_1^m \\
\vdots &\vdots &\vdots &\ddots &\vdots\\
1 &x_n &x_n^2 &\cdots &x_n^m \\
\end{pmatrix}
Y = \begin{pmatrix}
y_0\\
y_1\\
y_2\\
\vdots\\
y_n\\
\end{pmatrix}
A = \begin{pmatrix}
a_0\\
a_1\\
a_2\\
\vdots\\
a_m\\
\end{pmatrix}
$$
其中 $X$ 的维数为：$n*(m+1)$， $Y和A$ 的维数为：$n*1$, 则上面的矩阵可以写为：
$$
\begin{aligned}
&X^T(Y-XA) = 0 \\
&\Rightarrow Y = XA \\
&\Rightarrow A = (X^T X)^{-1}(X^T)Y
\end{aligned}
$$
即最高次数固定为$m$是的系数列矩阵为 $A$, 然后可以用 $f(x) = a_0 + a_1x^1 + a_2x^2\cdots + a_mx^m$ 拟合 $\{P_j\}$

下面是 不同m值拟合时的图形：

![](C:\Users\lch\AppData\Roaming\Typora\typora-user-images\image-20201014172612639.png)

![](C:\Users\lch\AppData\Roaming\Typora\typora-user-images\image-20201014172643074.png)



![image-20201014172748359](C:\Users\lch\AppData\Roaming\Typora\typora-user-images\image-20201014172748359.png)





#### 岭回归（Ridge Regression）

对上述最小二乘法误差函数增加 $E_1$ 正则项，参数 $\lambda$，$\min (E+\lambda E_1)$，其中 $E_1=\sum_{i=1}^n\alpha_i^2$ 
$$
\begin{aligned}
E + E_1 &=\sum_{i=0}^{n}(y_i-f(x_i))^2 =\sum_{i=0}^{n}(y_i - \sum_{j=0}^{m}{a_j}{x_i^j})^2 + \lambda \sum_{i=0}^{n}a_i^2 \\
 \frac {\partial (E+E_1)}{\partial a_k} &=  -2\sum_{i=0}^{n}{x_i^k}(y_i - \sum_{j=0}^{m}{a_j}{x_i^j}) + 2\lambda {a_k} \\
 

\end {aligned}
$$
 写成矩阵形式：
$$
\begin{pmatrix}
\frac{\partial (E+E_1)}{\partial a_0} \\
\frac{\partial (E+E_1)}{\partial a_1} \\
\vdots\\
\frac{\partial (E+E_1)}{\partial a_m} \\
\end{pmatrix}
= -2
\begin{pmatrix}
1 & 1 &1 &\cdots & 1\\
x_0 & x_1 & x_2 &\cdots & x_n\\
x_0^2 &x_1^2 &x_2^2 &\cdots &x_n^2\\
\vdots &\vdots &\vdots &\ddots &\vdots\\
x_0^m &x_1^m &x_2^m &\cdots &x_n^m\\
\end{pmatrix}
\begin{pmatrix}
y_0-\sum_{j=0}^{m}{a_j}{x_0^j}\\
y_1-\sum_{j=0}^{m}{a_j}{x_0^j}\\
\vdots\\
y_n-\sum_{j=0}^{m}{a_j}{x_0^j}\\
\end{pmatrix}
+ 2 \lambda
 \begin{pmatrix}
a_0\\
a_1\\
a_2\\
\vdots\\
a_m\\
\end{pmatrix}
= 
\begin{pmatrix}
0\\
0\\
\vdots\\
0
\end{pmatrix}
$$
同最小二乘法，令：
$$
X = 
\begin{pmatrix}
1 &x_0 &x_0^2 &\cdots &x_0^m \\
1 &x_1 &x_1^2 &\cdots &x_1^m \\
1 &x_2 &x_2^2 &\cdots &x_1^m \\
\vdots &\vdots &\vdots &\ddots &\vdots\\
1 &x_n &x_n^2 &\cdots &x_n^m \\
\end{pmatrix}
Y = \begin{pmatrix}
y_0\\
y_1\\
y_2\\
\vdots\\
y_n\\
\end{pmatrix}
A = \begin{pmatrix}
a_0\\
a_1\\
a_2\\
\vdots\\
a_m\\
\end{pmatrix}
$$
则：
$$
\begin{aligned}
&-2X^T(Y-XA) + 2\lambda A = 0 \\
&\Rightarrow X^T(Y - XA) = \lambda A \\
&\Rightarrow X^TY = (X^TX + \lambda I) A \\
&\Rightarrow A = (X^TX + \lambda I)^{-1}X^TY
\end{aligned}
$$
下面时 m = 7， 取不同 $\lambda$ 时的图形：

![](C:\Users\lch\AppData\Roaming\Typora\typora-user-images\image-20201014173329255.png)

![](C:\Users\lch\AppData\Roaming\Typora\typora-user-images\image-20201014173418169.png)

![image-20201014173453112](C:\Users\lch\AppData\Roaming\Typora\typora-user-images\image-20201014173453112.png)

![image-20201014173527480](C:\Users\lch\AppData\Roaming\Typora\typora-user-images\image-20201014173527480.png)



可以看出 $\lambda$ 取值需要合适，不能太大也不能太小



#### 四种方法图形比较

![image-20201014174602605](C:\Users\lch\AppData\Roaming\Typora\typora-user-images\image-20201014174602605.png)



## 总结

1. 多项式函数插值 在少量数据时可以得到较好的结果，但数据增多时，图形容易出现震荡现象。
2. 高斯基函数插值 图形和 参数 $\sigma$ 取值有关，当 $\sigma$ 较小时，图形会出现突刺现象，当 $\sigma$ 合适时，图形会变得更平滑
3. 固定最高次数多项式拟合 在阶数偏低时拟合效果不够，在阶数增加后，拟合效果越来越好，当阶数达到n-1时，拟合图形和多项式插值一样。
4. 岭回归中 $\lambda$ 取值对图形影响较大，当 $\lambda = 0$ 时，函数即是固定最高次数的最小二乘法得出的图形 

