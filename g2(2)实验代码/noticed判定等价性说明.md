# noticed 判定与三元组中心化子判定的等价性说明

本文说明下面两个条件为什么等价，并解释为什么代码里改用三元组中心化子的判定是合理的。下面所有数学内容都直接按公式显示。

## 1. 记号与问题

设有一个 theta-稳定约化子代数，以及其中的一个位于非紧部分的幂零元。记

$$
e\in \mathfrak l\cap \mathfrak p.
$$

再取一个与该幂零元相关联的 sl2 三元组

$$
(x,e,f)\subset \mathfrak l,\qquad [x,e]=2e,\ [x,f]=-2f,\ [e,f]=x.
$$

我们关心两个条件：

### 条件 A：定义型判定

对任意位于紧部分中的元素 $z$，若

1. $z$ 是半单元；
2. $z$ 与 $e$ 对易；

则必有

$$
z\in \mathfrak z(\mathfrak l).
$$

这正是图片里 noticed 的定义。

### 条件 B：三元组中心化子判定

$$
\mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f)=\mathfrak z(\mathfrak l)\cap \mathfrak k.
$$

代码里实际检查的是两边维数相等：

$$
\dim \mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f)=\dim (\mathfrak z(\mathfrak l)\cap \mathfrak k).
$$

因为总有包含关系

$$
\mathfrak z(\mathfrak l)\cap \mathfrak k\subseteq \mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f).
$$

所以在这里“维数相等”等价于“子代数相等”。

---

## 2. 为什么中心一定包含在三元组中心化子里

若某个元素位于中心与紧部分的交中，那么它与整个李代数中的所有元素都对易，特别与三元组中的三个元素都对易，因此

$$
z\in \mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f).
$$

所以自然有

$$
\mathfrak z(\mathfrak l)\cap \mathfrak k\subseteq \mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f).
$$

这就是代码里只比较维数就足够的根本原因。

---

## 3. 关键结构事实：三元组中心化子是单个幂零元中心化子的 Levi 因子

这是 Jacobson--Morozov 理论中的标准事实。

更准确地说，对与该幂零元相关联的 sl2 三元组，有：

$$
\mathfrak c_{\mathfrak l}(e)=\mathfrak c_{\mathfrak l}(x,e,f)\ltimes \mathfrak u.
$$

其中上式中的 $\mathfrak u$ 表示幂零根基，而三元组中心化子记作

$$
\mathfrak c_{\mathfrak l}(x,e,f).
$$

而上面出现的三元组中心化子是一个约化子代数，也就是单个幂零元中心化子的 Levi 因子。

把这个三元组中心化子再与紧部分相交后，得到

$$
\mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f).
$$

控制了

$$
\mathfrak c_{\mathfrak l\cap\mathfrak k}(e).
$$

中全部半单方向。换句话说：

- 中央化这个幂零元的半单元，都落在某个 Levi 因子里；
- 而这里最自然的 Levi 因子正是三元组中心化子。

因此，定义中“所有中央化该幂零元的半单元都在中心里”这个条件，可以转写为“这个 Levi 因子除了中心之外没有别的部分”。

---

## 4. 条件 B 推出条件 A

现在假设

$$
\mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f)=\mathfrak z(\mathfrak l)\cap \mathfrak k.
$$

取任意位于紧部分中的半单元 $z$，并假设它与幂零元对易。

由于这样的 $z$ 落在单个幂零元的中心化子里，而三元组中心化子正是其中控制半单方向的 Levi 因子，所以这里所有半单方向都已经包含在

$$
\mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f)
$$

里。

但根据假设，这个三元组中心化子恰好就是

$$
\mathfrak z(\mathfrak l)\cap \mathfrak k.
$$

所以这个半单元只能属于中心，即

$$
z\in \mathfrak z(\mathfrak l).
$$

这就是条件 A。

---

## 5. 条件 A 推出条件 B

反过来假设条件 A 成立，即任意位于紧部分中的半单元只要与该幂零元对易，就一定属于中心。

考虑

$$
\mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f).
$$

这是一个约化李代数。记它的导出代数为

$$
\mathfrak s=\big[\mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f),\ \mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f)\big].
$$

那么 $\mathfrak s$ 是半单李代数，而

$$
\mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f)
=
\mathfrak z\!\big(\mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f)\big)\oplus \mathfrak s.
$$

又因为

$$
\mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f)\subseteq \mathfrak c_{\mathfrak l\cap\mathfrak k}(e).
$$

所以其中任意半单元也都中央化该幂零元。由条件 A，这些半单元都必须落在中心与紧部分的交中。

现在来处理你指出的关键点：不仅要说明半单元落在

$$
\mathfrak z(\mathfrak l)\cap\mathfrak k,
$$

还要说明整个约化李代数本身就只能等于这个中心。

如果 $\mathfrak s\neq 0$，那么作为非零半单李代数，它必含有非零半单元；例如取 $\mathfrak s$ 的一个 Cartan 子代数，其中的元素都是半单元。由于

$$
\mathfrak s\subseteq \mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f)\subseteq \mathfrak c_{\mathfrak l\cap\mathfrak k}(e),
$$

这些半单元都中央化该幂零元，于是由条件 A 必须属于

$$
\mathfrak z(\mathfrak l)\cap\mathfrak k.
$$

但 $\mathfrak s$ 是半单李代数，它的中心为零，所以 $\mathfrak s$ 不可能含有这样的非零中心元。这与 $\mathfrak s\neq 0$ 矛盾。因此只能有

$$
\mathfrak s=0.
$$

于是

$$
\mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f)
=
\mathfrak z\!\big(\mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f)\big)
$$

是一个交换约化李代数；对代数群意义下的约化李代数，这个中心是 toral 的，因此其中不会再出现额外的幂零方向。再结合前面已经得到的包含关系

$$
\mathfrak z(\mathfrak l)\cap\mathfrak k
\subseteq
\mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f),
$$

以及“其中每个半单元都属于 $\mathfrak z(\mathfrak l)\cap\mathfrak k$”，便可推出整个中心化子只能有

$$
\mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f)=\mathfrak z(\mathfrak l)\cap\mathfrak k.
$$

这就是条件 B。

---

## 6. 为什么旧代码里的 Levi(C/Z)=0 不够精确

旧代码使用的是：

$$
C=\mathfrak c_{\mathfrak l\cap\mathfrak k}(e),\qquad Z=\mathfrak z(\mathfrak l)\cap\mathfrak k.
$$

然后检查

$$
\operatorname{Levi}(C/Z)=0.
$$

这个条件看起来很像上面的条件，但它仍然不如条件 B 直接，原因是：

- 它只是在商代数层面检查“是否还有非零 Levi 因子”；
- 但定义本身谈的是“所有中央化该幂零元的半单元是否都在中心里”；
- 最贴近这个定义的对象其实是三元组中心化子，也就是中央化子的 Levi 部分本身，而不是商以后再看 Levi 分解。

因此从实现角度，直接检查

$$
\mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f)=\mathfrak z(\mathfrak l)\cap\mathfrak k.
$$

更稳妥，也更不容易误判。

---

## 7. 对应到当前代码

当前代码采用的是下面这个等价判据：

$$
\dim \mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f)=\dim(\mathfrak z(\mathfrak l)\cap\mathfrak k).
$$

对应实现位置：

- `Subalgebra(L, [t[1], t[2], t[3]])`：构造由三元组生成的子代数；
- `LieCentralizer(L, triple_subalg)`：求其中心化子；
- 再与紧部分相交；
- 最后比较它与中心和紧部分交的维数。

也就是说，代码已经从“通过单个幂零元中心化子的商代数间接判断”改成了“通过三元组中心化子直接判断”。

---

## 8. 最终结论

在当前使用的 associated sl2 三元组语境下，下面两个条件是等价的。第一个条件中的 $z$ 预先假定为半单元，因此公式里不再额外写文字说明：

$$
\forall z\in \mathfrak l\cap\mathfrak k,\ [z,e]=0\Rightarrow z\in \mathfrak z(\mathfrak l).
$$

与

$$
\mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f)=\mathfrak z(\mathfrak l)\cap\mathfrak k.
$$

因此，代码中把 noticed 判定改写成“三元组中心化子恰好等于中心”的形式，是符合定义的。

---

## 9. 一个实现层面的备注

严格说，数学上最本质的是“子代数相等”；代码里使用“维数相等”只是因为已经知道包含关系

$$
\mathfrak z(\mathfrak l)\cap \mathfrak k\subseteq \mathfrak c_{\mathfrak l\cap\mathfrak k}(x,e,f).
$$

所以两边维数一旦相等，就自动得到子代数相等。
