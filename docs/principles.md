> 中文编写 英语不好 有英语好的可以翻译下提个PR 感谢

其实原理很简单

首先我们有一个球面 $ S_0 $ 半径为 $ R $

然后球上一点 $ P_0(\theta _0, \varphi _0) $ 当作初始位置

然后 我们从这个点作另一个球面 $ S $ 两球面交圆为 $ C $

该圆上的点到 $ S_0 $ 距离都为 $ S $的半径 $ r $

容易得到我们的实际绘图的距离中点的距离为 $ \rho = arcsin(r) $

现在我们将过点$ P_0 $ 与 $ S_0 $ 相切 且z轴正半轴相交的线作为我们绘图的y轴

设 $$ \theta_p, \varphi_p $$作为参考参数

容易得到

$$
r = \sqrt{sin^2(\theta_p) + sin^2(\varphi_p)} \\
\rho = arcsin(\sqrt{sin^2(\theta_p) + sin^2(\varphi_p)}) \\
\theta = arctan(\frac{sin(\varphi_p)}{sin(\theta_p)})
$$

代入求解

$$
\begin{aligned}
x &= \rho \cdot cos(\theta) \\
&= arcsin(\sqrt{sin^2(\theta_p) + sin^2(\varphi_p)}) \cdot cos(arctan(\frac{sin(\varphi_p)}{sin(\theta_p)})) \\
&= \frac{arcsin(\sqrt{sin^2(\theta_p) + sin^2(\varphi_p)})}{\sqrt{1+(\frac{sin(\varphi_p)}{sin(\theta_p)})^2}} \\
&= \frac{arcsin(\sqrt{sin^2(\theta_p) + sin^2(\varphi_p)})}{\sqrt{sin(\varphi_p)^2 + sin(\theta_p)^2}} \cdot sin(\theta_p) \\
&= \rho \cdot sin(\theta) \\
&= \frac{arcsin(\sqrt{sin^2(\theta_p) + sin^2(\varphi_p)})}{\sqrt{sin(\varphi_p)^2 + sin(\theta_p)^2}} \cdot sin(\varphi_p) \\
\end{aligned}
$$

得到最后的方程

$$
\begin{cases}
x = \frac{arcsin(\sqrt{sin^2(\theta_p) + sin^2(\varphi_p)})}{\sqrt{sin(\varphi_p)^2 + sin(\theta_p)^2}} \cdot sin(\theta_p) \\
y = \frac{arcsin(\sqrt{sin^2(\theta_p) + sin^2(\varphi_p)})}{\sqrt{sin(\varphi_p)^2 + sin(\theta_p)^2}} \cdot sin(\varphi_p)
\end{cases}
$$