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
