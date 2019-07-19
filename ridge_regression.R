library("lars")
library(ElemStatLearn)
library(plyr)

###mchol函数将对称方阵分解为一个下三角矩阵乘以该矩阵转置的形式,函数返回值为下三角矩阵
#输入：欲分解的矩阵x
#输出：cholesky分解所得矩阵L
mchol <- function(x)
{
  #求矩阵x的行列数,m为行数,n为列数
  mn <- dim(x)
  m <- mn[1]
  n <- mn[2]
  
  #检验x是否为方阵
  if(m != n) 
  {
    return ("Wrong dimensions of matrix!")
  }
  
  #检验x是否为对称矩阵
  if(sum(t(x) != x) > 0) 
  {
    return ("Input matrix is not symmetrical!")
  }
  
  #L为与x行列数相等的零矩阵，用于存放分解所得下三角矩阵
  L <- matrix(0, m, m)
  
  #循环每进行一次,求解一列矩阵L的元素
  #矩阵x第i列和第i行之前的元素不再使用，相当于矩阵x减少一个维数，故下述将循环所至第i列记为当前矩阵x和矩阵L的第一列
  for(i in 1:m)
  {
    #L的主对角线上第一个元素为x的主对角线上第一个元素开方
    L[i,i] <- sqrt(x[i,i])
    if(i < m)
    {
      #求当前矩阵L的第一列除第一个元素外的其他元素
      L[(i+1):m,i] <- x[(i+1):m,i]/L[i,i]
      
      #矩阵L第一列（除第一个元素）乘以它的转置得到TLM用于更新矩阵x，效果同TLM%*%TLM
      TLV <- L[(i+1):m,i]                               #记录已求出第一列除第一个元素外剩下元素
      TLM <- matrix(TLV, m-i, m-i)                      #TLV按列复制成矩阵
      TLM <- sweep(TLM, 2, TLV, "*")                    #sweep(x， MARGIN， STATS， FUN=”-“， …) 对矩阵进行运算
      #MARGIN为1，表示行的方向上进行运算，为2表示列的方向上运算(是指将参数从列的方向移下去算)
      #STATS是运算的参数，FUN为运算函数，默认是减法
      
      #减少一个维数的矩阵x更新为原来对应位置上的元素减去TLM，为下一次循环做准备
      x[(i+1):m,(i+1):m] <- x[(i+1):m,(i+1):m] - TLM
    }
  }
  #矩阵的返回值为我们要求的下三角矩阵L
  L  
}

###mforwardsolve函数求解线性方程租Lx=b，其中L为下三角矩阵
#输入：下三角矩阵L，向量b
#输出：线性方程组的解x
mforwardsolve <- function(L, b)
{
  #求L的行列数,m为L的行数,n为L的列数
  mn <- dim(L)
  m <- mn[1]
  n <- mn[2]
  
  #判断L是否为方阵
  if(m != n) 
  {
    return ("Wrong dimensions of matrix L!")
  }
  
  #判断L是否为下三角矩阵
  for (i in 1:(m-1))
  {
    if(sum(L[i,(i+1):m] != 0) > 0)#逐行判断上三角是否全为0元素
    {
      return ("Matrix L must be a lower triangular matrix!")
    }
  }
  
  #判断L的行数与b的长度是否相等
  if(m != length(b))
  {
    return ("Wrong dimensions of matrix L or vector b!")
  }
  
  #0向量记录求解结果
  x=rep(0, m)
  
  #循环每进行一次,求解一个x中的元素，看作矩阵L向量x向量b的维数减一
  #故下述将矩阵L的第i列记为当前矩阵L第一列，将向量x向量b的第i个元素记为当前向量第一个元素
  for(i in 1:m)
  {
    #求当前循环中x的第一个元素
    x[i] <- b[i] / L[i,i]
    #降维后的b向量为原来位置上的元素减去当前矩阵L的第一列的乘积
    if(i < m) 
    {
      b[(i+1):m] <- b[(i+1):m] - x[i]*L[(i+1):m,i]
    }      
  }
  #函数返回的x向量即为线性方程组的解
  x  
}

###mbacksolve函数求解线性方程租Lx=b，其中L为上三角矩阵
#输入：上三角矩阵L，向量b
#输出：线性方程组的解x
mbacksolve <- function(L, b)
{
  #求L的行列数,m为L的行数,n为L的列数
  mn <-dim(L)
  m <- mn[1]
  n <- mn[2]
  
  #判断L是否为方阵
  if(m != n)
  {  
    return ("Wrong dimensions of matrix L!")
  }
  
  #判断L是否为上三角矩阵
  for (i in 2:m)
  {
    if(sum(L[i,1:(i-1)] != 0) > 0)
    {
      return ("Matrix L must be a upper triangular matrix!")
    }
  }
  
  #判断L的行数与b的列数是否相等
  if(m != length(b))
  {
    return ("Wrong dimensions of matrix L or vector b!")
  }  
  
  x <- rep(0, m)
  
  #循环每进行一次,求解一个x中的元素，看作矩阵L向量x向量b的维数减一
  #故下述将矩阵L的第i列记为当前矩阵L最后一列，将向量x向量b的第i个元素记为当前向量最后一个元素
  for(i in m:1)
  {
    #求当前循环中x的最后一个元素
    x[i] <- b[i] / L[i,i]
    #降维后的向量b为原来位置上的元素减去刚才求出的x元素与当前上三角矩阵L最后一列（除最后一个元素）的乘积
    if(i > 1) 
    {
      b[(i-1):1] <- b[(i-1):1] - x[i]*L[(i-1):1,i]
    }      
  }
  #函数返回值x向量即为线性方程组的解
  x  
}

###ridgereg函数用于实现岭回归参数beta的估计，参数x和y分别为回归方程的自变量和因变量,lambda为L2正则项的调节参数
#此函数求解线性方程租(t(x)%*%x+lambada)%*%beta=t(x)%*%y,将t(x)%*%x+lambada进行cholesky分解为R%*%t(R),forwardsolve求解L%*%d=t(x)%*%y,其中d=t(R)%*%beta,backsolve求解t(R)%*%beta=d,即得参数beta的估计值 
#输入：自变量x，因变量y，调节参数lambda
#输出：回归系数beta的估计值
ridgereg <- function(lambda, x, y)
{
  #y=data[,m]; x=data[,-m]
  #n为自变量矩阵行数,即n个样本,p为自变量矩阵列数,即p个参数
  np <- dim(x)
  n <- np[1]
  p <- np[2]
  scale_data = scale_fun (x,y)
  #对x,y进行标准化
  y = scale(y,center=T,scale=T) 
  x = scale(x,center=T,scale=T) 
  #储存y，x的均值方差
 
  #标准化的数据没有截距项
  x = as.matrix(x)
  
  
  #利用cholesky分解求取回归方程的参数beta的估计值  
  V <- t(x)%*%x + diag( rep(lambda, p))              #t(x)%*%x+lambda作为线性方程组的系数矩阵V
  U <- as.vector(t(x)%*%y)                           
  R <- mchol(V)                                           #调用mchol函数将系数矩阵V进行cholesky分解,V=R%*%t(R)
  M <- mforwardsolve(R, U)                                #使用前代法求解R%*%M=t(x)%*%y,其中M=t(R)%*%beta
  betai=mbacksolve(t(R), M)*()                                     #使用回代法求解t(R)%*%beta=M,即可得beta的估计值
  beta0=
}

scale_fun <- function(x, y)
{
  #储存y，x的均值方差
  y_mean = mean(y)
  y_var = var(y)
  x_mean = colMeans(x)
  x_var =apply(x,2,var)
  scale_data = list( y_mean, y_var,x_mean, x_var)
  scale_data
}

pred <- function(b, nx)
{
  #nx=prostate[1:2,1:8]
  b <- as.vector(b)
  p <- length(b) - 1
  
  #将数据矩阵nx重新排列，每一行为一个样品，重排矩阵的原因是下面例子中调用的数据原结构为dataframe
  nx <- as.matrix(nx, ncol <- p)
  n <- dim(nx)[1]
  
  #计算预测值
  apply(t(nx)*b[2:(p+1)], 2, sum) + b[1]  
}

###mridge函数用于实现删去一个样品的岭回归
mridge=function(i,lambda,x,y) 
{
  ridgereg(lambda,x[-i,],y[-i])
}

###cvridgeregerr函数用交叉验证实现岭回归，参数依次为调节参数lambda,自变量x(数据矩阵),因变量y,返回值为测试均方误差   
#输入：超参数lambda,自变量x(数据矩阵),因变量y
#输出：删一交叉验证岭回归测试均方误差
cvridgeregerr<-function(lambda,x,y){  
  #lambda=1
  np<-dim(x)
  n<-np[1]
  p<-np[2]
  #矩阵中的元素作为第一个参数输入mridge，表示去掉的数据编号，结果第i行为删去第i个样本的岭回归系数估计值
  coe<-t(apply(as.matrix(1:n,ncol=1),1,mridge,lambda,x,y))
  #储存每次删除后数据的均值方差
  scale_data<-apply(as.matrix(1:n,ncol=1),1,scale_fun,x,y)
  
  #coe第i行和数据矩阵第i个样本做点对点相乘，对行求和，计算测试均方误差
  mean((apply(coe*x,1,sum)-y)^2)
}

###ridgeregerr函数用于计算训练均方误差  
#输入：岭回归超参数lambda，数据矩阵x，因变量y
#输出：训练均方误差
ridgeregerr=function(lambda,x,y)
{
  mean((pred(ridgereg(lambda,x,y),x)-y)^2)
}   

###############################
###在不同的lambda下,比较训练均方误差和测试均方误差，以选取合适的调节参数lambda
library(ElemStatLearn)
x <- as.matrix(prostate[ ,1:8])
y <- as.vector(prostate[ ,9])
LAM <- seq(0.001, 10, len=50)
#计算岭回归50个模型的训练均方误差，将结果从list展开成向量
err <- unlist(lapply(LAM, ridgeregerr, x, y))
#计算岭回归50个模型的测试均方误差，将结果从list展开成向量
pe <- unlist(lapply(LAM, cvridgeregerr, x, y))
x <- rep(1:50, 2)
plot(pe)
#取交叉验证中使测试均方误差最小的lambda
lam=LAM[which.min(pe)]
xsd = apply(x,2,var)
