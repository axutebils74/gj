# 浏览器需要支持 Fetch Webassembly  
# 引用 
# giac.wasm v1.9.0
[https://www-fourier.ujf-grenoble.fr/~parisse/giac/](https://www-fourier.ujf-grenoble.fr/~parisse/giac/)
[http://www-fourier.ujf-grenoble.fr/~parisse/giac.html](http://www-fourier.ujf-grenoble.fr/~parisse/giac.html)
# LZMA
[https://github.com/jcmellado/js-lzma](https://github.com/jcmellado/js-lzma)
# katex
[https://github.com/KaTeX/KaTeX](https://github.com/KaTeX/KaTeX)
# 使用方法
[http://www-fourier.ujf-grenoble.fr/~parisse/giac/doc/en/casinter/index.html](http://www-fourier.ujf-grenoble.fr/~parisse/giac/doc/en/casinter/index.html) 
# 不会的可以使用help命令.
|     |     |     |     |
| --- | --- | --- | --- |
| 函数  | 作用  | 例子  | 结果  |
| **改写** |     |     |     |
| simplify | 简化表达式 | simplify(4*atan(1/5)-atan(1/239)) | 1/4 * π |
| normal | 展开表达式 | normal((2*x+1)^2) | 4\*x^2+4\*x+1 |
| factor | 因式分解 | factor(x^4-1)<br>factor(x^4-4,sqrt(2)) | $(x-1)\cdot (x+1) (x^{2}+1)$<br>$(x-\sqrt{2})\cdot (x+\sqrt{2}) (x^{2}+2)$ |
| partfrac | 部分分式 | partfrac(x/(x^2-4)) | $\frac{1}{2\cdot (x-2)}+\frac{1}{2\cdot (x+2)}$ |
| subst | 值代换 | subst(x^2+y^2,x=4,y=0)<br>或者<br>subst(x^2+y^2,\[x=4,y=0\])<br>或者<br>subst(x^2+y^2,\[x,y\],\[4,0\]) | 16  |
| reorder | 变量重排 | reorder(x^2+y^3,\[y,x\]) | y^3+x^2 |
| texpand | 对数，指数，三角函数展开 | texpand(sin(2\*x)+exp(x+y)+ln(x\*y)) | 2\*cos(x)\*sin(x)+exp(x)*exp(y)+ln(y)+ln(x) |
| lin | 指数函数线性化 | lin（$\left(e^{x}+e^{x^{2}}\right)^{2}$） | $e^{2\cdot x}+2 e^{x^{2}+x}+e^{2\cdot x^{2}}$ |
| tlin | 三角函数线性化用于积分 | tlin(sin(x)^3)<br>tlin(cos(x)*cos(y)) | 3/4\*sin(x)-1/4\*sin(3*x)<br>1/2\*cos(x-y)+1/2\*cos(x+y) |
| tcollect | 三角函数线性化用于值域 | tcollect(sin(x)+cos(x)) | $\sqrt{2} \cos\left(x-\frac{1}{4}\cdot \pi \right)$ |
| trigsin | 用sin x改写函数 ,改写优先级 sin x > cos x > tan x | trigsin(cos(x)^4+sin(x)^2) | $\sin\left(x\right)^{4}-\sin\left(x\right)^{2}+1$ |
| trigcos | 用cos x改写函数 ,改写优先级 cos x >sin x > tan x | cos(1) | $\cos\left(x\right)^{4}-3 \cos\left(x\right)^{2}+2$ |
| trigtan | 用tan x改写函数 ,改写优先级 tan x >sin x > tan x | sin(x)^4 + sin(x)^2 | $\frac{(2 \tan\left(x\right)^{4}+\tan\left(x\right)^{2})}{(\tan\left(x\right)^{4}+2 \tan\left(x\right)^{2}+1)}$ |
| halftan | 半角公式改写 | halftan(sin(x)) | $\frac{2 \tan\left(\frac{x}{2}\right)}{(\tan\left(\frac{x}{2}\right)^{2}+1)}$ |
| trig2exp | 欧拉公式将三角函数改写为指数 | trig2exp(sin(x)) | $\frac{(e^{i\cdot x}-e^{-i\cdot x})}{2*i}$ |
| exp2trig | 指数改写三角函数 | exp2trig(e^(i*x)) | $\cos\left(x\right)+i \sin\left(x\right)$ |
| lncollect | 对ln函数进行合并因式 | lncollect(ln(x)+ln(y)) | ln(x*y) |
| exp2pow | 改写形如 e^(ln(x) * 3) 这样的式子 | exp2pow(e^(ln(x) * 3)) | x^3 |
| pow2exp | exp2pow逆过程 | pow2exp(x^3) | e^(3*ln(x)) |
| **计算** |     |     |     |
| limit | 求极限 limit(fn,a,b,d)<br>fn: 表达式<br>a: 变量 默认x<br>b: 值 默认 0<br>d: 方向  (0 两侧) (1 +)  (-1 -) 默认两侧 | limit(sin(x)/x)<br>limit(sin(x*y)/y,y)<br>limit(arctan(x),x,+inf)<br>limit(arctan(1/x),x,0,-1) | 0<br>x<br>π / 2<br>-π / 2 |
| diff | 求导数 diff(fn,a) 或者用∂<br>fn: 表达式<br>a: 变量 默认x | diff(ln(x))<br>diff(ln(x-y),y)<br>(x^2+3x)'<br>∂(x^2+x)<br>∂(x^2+y,y) | 1/x<br>-1/(x-y)<br>2*x+3<br>2*x<br>1 |
| integrate | 求积分 integrate(fn,a,l,r) 或者用 ∫<br>fn: 表达式<br>a: 变量 默认x<br>l:下限 <br>r:上限<br>或者使用int 为了区分取整函数，必须填写变量旧版本因为没有micropython时不需要填写 | integrate(cos(x))<br>integrate(cos(x-y),y)<br>integrate(1/(1+x^2),x,-inf,+inf)<br>∫(cos(x))<br>∫(cos(x-y),y)<br>int(16.5,x)<br>int(16.5) | sin(x)<br>-sin(x-y)<br>π<br>sin(x)<br>-sin(x-y)<br>16.5x<br>16 |
| series | 求给定变量值附近的级数展开式<br>series(fn,a=b,c,d) 或者 series(fn,a,b,c,d) 或者<br>series(fn,a,b,c,d)<br>fn: 表达式<br>a: 变量 默认x<br>b: 值<br>c:阶数 默认5<br>d:方向  (0 两侧) (1 +)  (-1 -) 默认两侧<br>如果只有一个变量，类似于taylor,<br>series可以用来求多元函数的泰勒展开<br>series(fn,\[a...\],\[b...\],c,d) 或者<br>series(fn,\[a=b,....\],c,d)<br>不会出现余项order size<br>series最后一个参数为polynom时不会出现余项order size | (1) series(1/(1-x)) = series(1/(1-x),x)<br>(2) series(1/(1-y),y=2)<br>(3) series(1/(1-x),x=0,2)<br>(4) series(arctan(1/x),x=0,2,1)<br>(5) series(1/(1+x+y),\[x,y\],\[0,0\],2)<br>(6) series(1/(1+x+y),\[x=0,y=0\],2)<br>(7) series(1/(1-x),polynom)<br>(8) series(1/(1-x),polynom)(6)<br>(9) series(1/(1-x-y),\[x=0,y=0\],2)(x=1) | $(1)1-x+x^{2}-x^{3}+x^{4}-x^{5}+x^{6} \mathrm{order\_size}\left(x\right)$<br>$(2)-1+y-2-\left(y-2\right)^{2}+\left(y-2\right)^{3}-\left(y-2\right)^{4}+\left(y-2\right)^{5}+\left(y-2\right)^{6} \mathrm{order\_size}\left(y-2\right)$<br>$(3)1+x+x^{2}+x^{3} \mathrm{order\_size}\left(x\right)$<br>$(4)\frac{1}{2} \cdot \pi -x+x^{3} \mathrm{order\_size}\left(x\right)$<br>$(5)1-x-y+x^{2}+2\cdot x\cdot y+y^{2}$<br>$(6)1-x-y+x^{2}+2\cdot x\cdot y+y^{2}$<br>$(7)1+x+x^{2}+x^{3}+x^{4}+x^{5}$<br>$(8)9331$<br>$(9)2+y+1+2\cdot y+y^{2}$ |
| sum | 求和 sum(fn,a,c,d,e) <br>或者∑(fn,a,c,d)<br>fn: 表达式<br>a: 变量 默认x<br>c:下标<br>d:上标<br>e:步长<br>相当于 $\sum_{a=c}^{d}fn$ | sum(x,x,1,6)<br>sum(x,x,5,1) = -sum(x,x,2,4) <br>sum(x,x) = sum(x,x,1,x-1)+c<br>sum(\[1,2,3,4,5\])<br>a:=makelist(x->x,1,5);sum(a) | 21<br>9<br>x*(x-1)/2<br>15<br>\[1,2,3,4,5\],15 |
| desolve | 求解微分方程<br>desolve(a,b,c)<br>a: 表达式<br>b: 自变量默认为x 可省略<br>c: 因变量默认为y | desolve(y'+y=x)<br>desolve(x''+x=0,x)<br>desolve(x''+x=y,y,x)<br>desolve(\[y''+y=0,y(0)=1,y'(0)=1\]) | c_0*exp(-x)+x-1<br>c\_0\*cos(x)+c\_1\*sin(x)<br>c\_0\*cos(y)+c\_1\*sin(y)+y<br>cos(y)+sin(y) |
| ilaplace | 逆拉普拉斯 | ilaplace(1/(x^2))<br>ilaplace(1/(s^2),s)<br>ilaplace(1/(s^2),s,x) | x<br>s<br>x |
| laplace | 拉普拉斯用于微分方程 | laplace(sin(x))<br>laplace(sin(s),s)<br>laplace(sin(x)) | 1/(1+x^2)<br>1/(1+s^2)<br>1/(1+x^2) |
| odesolve | 常微分求解 | odesolve(x+y,\[x,y\],\[1,2\],2)<br>相当于 desolve(\[y'=x+y,y(1)=2\])(2) | 7.87312731383 |
| ifft | 逆傅里叶变换 | ifft(10,-2+2\*i,-2,-2-2\*i) | \[1,2,3,4\] |
| fft | 傅里叶变换 | fft(\[1,2,3,4\])<br>ifft(fft(\[1,2,3\])+fft(\[2,3,4\])) | \[10,-2+2\*i,-2,-2-2\*i\]<br>\[3,5,7\] |
| invztrans | 逆z变换 | invztrans(x/x-2) | 2^x |
| ztrans | z变换 | ztrans(e^x) | x/(-e+x) |
| grad | 求梯度<br>可以看做把x当变量，y当变量，求导数 | grad(x^2+x*y,\[x,y\]) | \[2*x+y,x\] |
| hessian | 求海森矩阵<br>可以看做先对第i个变量求偏导，再对第j个变量求偏导 | hessian(x*sin(y),\[x,y\]) | \[\[0,cos(y)\],\[cos(y),-x*sin(y)\]\] |
| curl | 求旋度 | curl(\[x,y,z*x\],\[x,y,z\]) | \[0,-z,0\] |
| divergence | 求散度<br>x+y 对x导数 + x*y 对y导数 | divergence(\[x+y,x*y\],\[x,y\]) | x+1 |
| vpotential | 找一个函数使得该函数的旋度为所给的值 | vpotential(\[2\*x\*y+3,x^2-4\*z,-2\*y*z\],\[x,y,z\]) | \[0,-2\*x\*y\*z,-1/3\*x^3+4\*x\*z+3*y\] |
| potential | 找一个函数使得该函数的梯度为所给的值 | potential(\[2*x+y,x\],\[x,y\]) | x^2+x*y |
| **求解** |     |     |     |
| solve | 解方程 | solve(x^2=3)<br>solve(\[y-z=0,z-x=0,x-y=0,x-1+y+z=0\],\[x,y,z\]) | \[-sqrt(3),sqrt(3)\]<br>list\[\[1/3,1/3,1/3\]\] |
| **画图** |     |     |     |
| plot | 函数画图 | plot(x^2) | 建议 help(plot) |
| implicitplot | 隐函数画图 | implicitplot(x^2+y^2=1) | 建议 help(implicitplot) |
| plotparam | 参数方程画图 | plotparam(\[sin(t)+1,cos(t)\]) 或者plotparam(sin(t)+1+i*cos(t)) | 建议 help(plotparam) |
| plotpolar | 极坐标画图 | plotpolar(1-sin(x)) | 建议 help(plotpolar) |
| plotfield | 向量场 | plotfield(x+y) | 建议 help(plotfield) |
| plotode |     | plotode(sin(t*y),\[t,y\],\[0,1\]) | 建议 help(plotode) |
| plotseq | 蛛网图 | plotseq(2*x,1,3); | 建议 help(plotseq) |
| plotlist | 线段  | plotlist(\[\[0,3\],\[2,1\],\[3,4\],\[0,3\]\]) | 建议 help(plotlist) |
| boxwhisker | 画盒图 | boxwhisker(\[\[6,0,1,3,4,2,5\],\[0,1,3,4,2,5,6\],\[1,3,4,2,5,6,0\],\[3,4,2,5,6,0,1\],\[4,2,5,6,0,1,3\],\[2,5,6,0,1,3,4\]\]) | 建议 help(boxwhisker) |
| histogram | 直方图 | histogram(\[1,2,1,1,2,1,2,4,3,3\]); | 建议 help(histogram) |
| scatterplot | 散点图 | scatterplot(\[\[1,2,3\],\[2,0,1\],\[-1,2,3\]\]) | 建议 help(scatterplot) |
| polygonplot |     | polygonplot(\[\[1,2,3\],\[2,0,1\],\[-1,2,3\]\]) | 建议 help(polygonplot) |
| linear\_regression\_plot | 线性回归 | linear\_regression\_plot(\[\[0.0,0.0\],\[1.0,1.0\],\[2.0,4.0\],\[3.0,9.0\],\[4.0,16.0\]\]) | 建议 help(linear\_regression\_plot) |
| bar_plot | 柱状图 | bar_plot(\[3/2,2/3,5/4,4/5,7/6,6/7,9/8,8/9,11/10\]) | 建议 help(bar_plot) |
| plotcdf |     | plotcdf(binomial,10,0.5)<br>plotcdf(normald,0.0,1.0)<br>plotcdf(\[1,3,4,3,5,6\]) | 建议 help(plotcdf) |
| plotcontour | 等高线 | plotcontour(x^2+y^2) | 建议 help(plotcontour) |
| plotarea |     | plotarea(sin(x),x=0..pi); | 建议 help(plotarea) |
| bezier | 贝塞尔曲线 | bezier(1,1+i,2+i,3-i,plot) | 建议 help(bezier) |
| **几何** |     |     |     |
| point | 画点  | point(1,2) 或 point(1+2*i) |     |
| segment | 画线段 | segment(point(i),point(1+3*i))<br>或 segment(i,1+3*i) |     |
| line | 画直线 | line(point(i),point(1+3*i)) 或 line(i,1+3*i)<br>line(2*x+y=0) |     |
| vector | 画向量<br>如果没有第一个参数默认为原点 | vector(point(-i),point(1+i)) 或 vector(-i,1+i) |     |
| equation | 无法使用 |     |     |
| parameq | 返回该函数参数方程 | parameq(circle(0,1)) | e^it (注：即x = cos(t),y=sin(t)) |
| midpoint | 画中点 | midpoint(-2,2i) |     |
| circle | 画圆<br>表达式画圆<br>过圆心画圆<br>过直径两点画圆<br>可以画弧 | circle(x^2+y^2=1)<br>circle(0,point(2))<br>circle(0,2)<br>circle(0,i,pi/4,pi/2); |     |
| isobarycenter | 画质心 | isobarycenter(-1,1-i,i) |     |
| center | 画圆心 | center(x^2+y^2=2)<br>center(circumcircle(0,1,1+i)) |     |
| circumcircle | 三点画圆 | circumcircle(0,1,1+i) |     |
| inter | 画交点 | inter(line(i,1-i),circle(0,1)) |     |
| triangle | 画三角形 | triangle(1,2+i,3) |     |
| rectangle | 画矩形，包括旋转的矩形 | rectangle(-i,1,2) |     |
| square | 画正方形，不包括旋转 | square(i,i-1) |     |
| isopolygon | 正多边形作图<br>isopolygon(a,b,c)<br>c为正数,过ab两点的多边形作图<br>c不为正数，以a为中心，b为正多边形其中一点 | isopolygon(0,1,5)<br>isopolygon(0,1,-5) |     |
| polygon | 定点集画多边形 | polygon(0,1+i,2) |     |
| parallel | 过一点画与已知直线平行的直线 | parallel(1,line(2,3+i)) |     |
| perpendicular | 过一点画与已知直线垂直的直线 | perpendicular(1,line(2,3+i)) |     |
| perpen_bisector | 画线段的平分线，或者两点的平分线 | perpen_bisector(1-i,i) |     |
| bisector | 给定三点ABC(AB,AC)画角BAC平分线 | bisector(0,1,i) |     |
| altitude | 给定三点ABC，画过点A与BC垂直的直线 | altitude(-1,1-i,i) |     |
| tangent | 给定曲线和一点，画出该曲线的切线且经过该点 | tangent(circle(i,1+i),point(-2)) |     |
| distance | 求两点的距离 | distance(point(0),point(1+i)) | sqrt(2) |
| **算术** |     |     |     |
| gcd | 最大公因数 | gcd(16,28,24)<br>gcd(x*(x-2),(x-2)*(x-6)) | 4<br>x-2 |
| iquo | 欧式除法  <br>相当于 a // b | iquo(7,2)<br>7//2<br>iquo(7,2+i) | 3<br>3<br>3-i |
| quo | 多项式的欧式除法 | quo(x^3+2x^2+3x+4,-x+2) <br>quo(\[1,2,3,4\],\[-1,2\])<br>quo(t^3+2t^2+3t+4,-t+2,t) | -poly1\[-1,-4,-11\]<br>-x^2-4*x-11<br>-t^2-4*t-11 |
| irem | 欧式余数<br>相当于a%b | irem(125,15)<br>irem(-7,3)<br>irem(25+12\*i,5+7\*i) | 5<br>2<br>-4+i |
| rem | 多项式的欧式余数 | rem(\[1,2,3,4\],\[-1,2\])<br>rem(x^3+2x^2+3x+4,-x+2)<br>rem(t^3+2t^2+3t+4,-t+2,t) | poly1\[26\]<br>26<br>26 |
| iabcuv | 给定三个整数a,b,c，找到2个整数u，v使得<br>a\*u+b\*v=c<br>u，v组合有无数种,<br>c必须为a，b的最大公因数的整数倍 | iabcuv(21,28,7)<br>iabcuv(21,28,1) | \[-1,-1\]<br>环中无解 |
| abcuv | 给定三个多项式a,b,c，找到2个多项式u，v使得<br>a\*u+b\*v=c | abcuv(x^2+2*x+1,x^2-1,x+1)<br>abcuv(X^2+2*X+1,X^2-1,X^3+1,X) | \[1/2,-1/2\]<br>\[1/2*(-X+2),1/2\*3\*X\] |
| ifactor | 整数分解 | ifactor(4095) | 3^2\*5\*7*13 |
| lcm | 最小公倍数 | lcm(6,4,2)<br>lcm(x*(x+1),(x+1)*(x+2)) | 12<br>x*(x+1)*(x+2) |
| ichinrem | 中国余数定理<br>x%7 = 2<br>x%5 = 3<br>返回a,v<br>使得a+b*k对任意整数k都成立<br>其中a是最小整数解<br>b是最小公倍数 | ichinrem(\[2,7\],\[3,5\]) | \[23,35\] |
| isprime | 判断是否为质数 | isprime(2**31-1) | True |
| nextprime | 得到大于某个数的下一个质数 | nextprime(2^32-1) | 4294967311 |
| powmod | 快速幂求模，给定三个数a,b,c<br>求a^b % c<br>快于a^b % c | powmod(5,2,13) | 12  |
| mod | 求模<br>相当于 %<br>1 mod 3 | 7 mod 3 | 1   |
| GF  | 求伽罗瓦域 |     |     |
| lcoeff | 求多项式的最高阶 | lcoeff(11\*x^3+2\*x+8) | 11  |
| coeff | 返回多项式各个阶数的系数，如果最后有参数k则返回k阶的值 | coeff(x*3+2)<br>coeff(5*y^2-3,y)<br>coeff(5*y^2-3,y,2) | poly1\[3,2\]<br>poly1\[5,0,-3\]<br>5 |
| degree | 返回多项式阶数 | degree(x^3+1) | 3   |
| horner | 用霍纳方法计算多项式值<br>可以理解为把3带入x^2+1中 | horner(x^2+1,3) | 10  |
| ptayl | 给定多项式P和a<br>返回泰勒多项式Q<br>使得P(x)=Q(x-a)<br>可以理解为把x+a带入多项式中 | ptayl(x,1) | x+1 |
| proot | 多项式求根<br>建议使用solve(x^2+x-2=0)<br>或者<br>fsolve(x^2+x-2=0) | proot(x^2+x-2)<br>proot(\[1,1,-2\]) | \[-2,1\]<br>\[-2,1\] |
| **线性代数** |     |     |     |
| dot | 两个向量的点乘<br>相当于* | dot(\[1,2\],\[3,4\])<br>\[1,2\]*\[3,4\] | 11<br>11 |
| cross | 两个向量的叉乘 | cross(\[1,2\],\[3,4\])<br>cross(\[1,2,3\],\[3,4,5\]) | \[0,0,-2\]<br>\[-2,4,-2\] |
| tran | 矩阵的转置或者trn | tran(\[\[1,2\],\[3,4\]\]) | matrix\[\[1,3\],\[2,4\]\] |
| identity | 返回给定维数的单位矩阵<br>可简写为idn | identity(2)<br>idn(2) | matrix\[\[1,0\],\[0,1\]\]<br>matrix\[\[1,0\],\[0,1\]\] |
| makemat | 制作矩阵makemat(fn,x,y),<br>相当于matrix(x,y,fn)<br>x : 行<br>y：列<br>fn：lambda 表达式默认为 i和j | makemat((x,y)->(x-y),3,3)<br>matrix(3,3,(x,y)->(x-y)) | matrix\[\[0,-1,-2\],\[1,0,-1\],\[2,1,0\]\]<br>matrix\[\[0,-1,-2\],\[1,0,-1\],\[2,1,0\]\] |
| hilbert | 返回n维希尔伯特矩阵<br>相当于makemat((x,y)->(1/(x+y+1)),3,3) | hilbert(3)<br>makemat((x,y)->(1/(x+y+1)),3,3) | matrix\[\[1,1/2,1/3\],\[1/2,1/3,1/4\],\[1/3,1/4,1/5\]\] |
| vandermonde | 给定向量\[a,b,c,..\] 返回范德蒙德行列式<br>相当于makemat((x,y)->(\[a,b,c\]\[x\]^y),3,3) | vandermonde(\[1,2,3\])<br>makemat((x,y)->(\[1,2,3\]\[x\]^y),3,3) | matrix\[\[1,1,1\],\[1,2,4\],\[1,3,9\]\]<br>matrix\[\[1,1,1\],\[1,2,4\],\[1,3,9\]\] |
| rref | 化成行最简形矩阵<br>经过初等变换使得每列的第一个非零数字都是1，而且每列的第一个非零数字的右方都是零 | rref(\[\[3,1,-2\],\[3,2,2\]\]) | \[\[1,0,-2\],\[0,1,4\]\] |
| ker | 求矩阵的零空间<br>满足AX = 0 | ker(\[\[1,2\],\[3,6\]\])<br>\[\[1,2\],\[3,6\]\] * trn(\[\[2,-1\]\]) | \[\[2,-1\]\]<br>\[0\]<br>matrix\[\[0\],\[0\]\] |
| det | 行列式的值 | det(\[\[1,2\],\[3,4\]\]) | -2  |
| linsolve | 线性方程求解<br>建议用solve(\[x+y=1,x-y=2,\[x,y\]\]) | linsolve(\[x+y=1,x-y=2\],\[x,y\])<br>solve(\[x+y=1,x-y=2,\[x,y\]\]) | \[3/2,-1/2\]<br>list\[\[3/2,-1/2\]\] |
| lu  | 给定矩阵A进行LU分解<br>将矩阵分解为下三角L，上三角U两种形式 | lu(\[\[6,12,18\],\[5,14,31\],\[3,8,18\]\])<br>\[\[1,0,0\],\[2,1,0\],\[5/3,-1/6,1\]\]*\[\[3,8,18\],\[0,-4,-18\],\[0,0,-2\]\]<br>\[\[0,0,1\],\[1,0,0\],\[0,1,0\]\]*\[\[6,12,18\],\[5,14,31\],\[3,8,18\]\] | \[2,0,1\],\[\[1,0,0\],\[2,1,0\],\[5/3,-1/6,1\]\],\[\[3,8,18\],\[0,-4,-18\],\[0,0,-2\]\]<br>\[\[3,8,18\],\[6,12,18\],\[5,14,31\]\]<br>\[\[3,8,18\],\[6,12,18\],\[5,14,31\]\] |
| qr  | 给定矩阵A进行QR分解<br>q是酉矩阵没有虚数时可认为是正定矩阵，R是上三角矩阵 | simplify(qr(\[\[1,2\],\[3,4\]\]))<br>a:=\[\[√10/10,3*√10/10\],\[3*√10/10,-√10/10\]\];simplify(a*tran(a))<br>a*\[\[√10,7*√10/5\],\[0,√10/5\]\] | matrix\[\[√10/10,3*√10/10\],\[3*√10/10,-√10/10\]\],\[\[√10,7*√10/5\],\[0,√10/5\]\]<br>matrix\[\[1,0\],\[0,1\]\]<br>\[\[1,2\],\[3,4\]\] |
| cholesky | cholesky分解<br>找到一矩阵L使得A=L*tran(L) | cholesky(\[\[3,1\],\[1,4\]\])<br>simplify( \[\[3*√3/3,0\],\[√3/3,11/3*√33/11\]\]*trn(\[\[3*√3/3,0\],\[√3/3,11/3*√33/11\]\])) | \[\[3*√3/3,0\],\[√3/3,11/3*√33/11\]\]<br>\[\[3,1\],\[1,4\]\] |
| cond | 请用help(cond) |     |     |
| egv | 求可对角化矩阵（对角矩阵的相似矩阵，可由初等变换变得）中的可逆矩阵P<br>也就是使得P^-1AP为对角矩阵的P | Q:=egv(\[\[-2,-2,1\],\[-2,1,-2\],\[1,-2,-2\]\])<br>Q^-1*\[\[-2,-2,1\],\[-2,1,-2\],\[1,-2,-2\]\]*Q | \[\[1,-3,-3\],\[-2,0,-3\],\[1,3,-3\]\]<br>\[\[3,0,0\],\[0,-3,0\],\[0,0,-3\]\] |
| egvl | 请用help(egvl) |     |     |
| jordan | 请用help(jordan) |     |     |
| pcar | 请用help(pcar) |     |     |
| pmin | 返回矩阵a的最小多项式 | pmin((sqrt(5)-1)/2)<br>solve(x**2+x-1) | x**2+x-1<br>list\[1/2*(-√5-1),1/2*(√5-1)\] |
| rank | 求矩阵的秩 | rank(\[\[1,1\],\[1,1\]\]) | 1   |
