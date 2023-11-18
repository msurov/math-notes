# Segment point distance
The distance between a point $q$ and line segment $p(t)$:

$$p=\left(p_{2}-p_{1}\right)t+p_{1}$$

is the minimum of the function

$$ d^{2}=\left(\left(p_{2}-p_{1}\right)t+p_{1}-q\right)^{T}\left(\left(p_{2}-p_{1}\right)t+p_{1}-q\right) $$

The neseccary condition is

$$ \frac{dd^{2}}{dt}=2\left(\left(p_{2}-p_{1}\right)t_{*}+p_{1}-q\right)^{T}\left(p_{2}-p_{1}\right)=0 $$

express $t_{*}$ from the above

$$t_{*}=\frac{\left(q-p_{1}\right)^{T}\left(p_{2}-p_{1}\right)}{\left(p_{2}-p_{1}\right)^{T}\left(p_{2}-p_{1}\right)} $$

Finally, clamp the value if $t_{*}$ by $[-1,1]$ and substitute into the distance formula

$$d_{\*}=\left\Vert \left(p_{2}-p_{1}\right)t_{*}+p_{1}-q\right\Vert $$

# Plane Closest Point

We work in Euclidean space with the scalar product defined as $\left\langle v,w\right\rangle \equiv v^{T}Pw$.

1. There is a plane defined in implicit form $ax+by+cz+d=0$ or $\left\langle n,p-p_{0}\right\rangle =0$. $n$ is a unit vector which is normal to the plane: $\left\langle n,n\right\rangle =1$.
2. The relation betwen $n$ and coefficients $\left(a,b,c\right)$ is

$$n=\lambda P^{-1}\left(\begin{array}{c}
a\\
b\\
c
\end{array}\right) \quad \text{where} \quad \lambda=\pm\frac{1}{\sqrt{\left(\begin{array}{c}
a\\
b\\
c
\end{array}\right)^{T}P^{-1}\left(\begin{array}{c}
a\\
b\\
c
\end{array}\right)}} $$

   is responsible for normalization. This relation follows from the plane equation 

$$\left\langle \lambda P^{-T}\left(\begin{array}{c}
a\\
b\\
c
\end{array}\right),p-p_{0}\right\rangle =\lambda\left(\begin{array}{c}
a\\
b\\
c
\end{array}\right)^{T}P^{-1}P\left(p-p_{0}\right)=\lambda\left(\begin{array}{c}
a\\
b\\
c
\end{array}\right)^{T}\left(p-p_{0}\right)=0$$

3. The closest point of the plane $q_{*}$ shich is closest to the given point $q$ is

$$ q_{*} = \left[I-vv^{T}P\right]\left(q-p_{0}\right)+p_{0} $$
