library("lars")
library(ElemStatLearn)
library(plyr)

###mchol�������ԳƷ���ֽ�Ϊһ�������Ǿ�����Ըþ���ת�õ���ʽ,��������ֵΪ�����Ǿ���
#���룺���ֽ�ľ���x
#�����cholesky�ֽ����þ���L
mchol <- function(x)
{
  #�����x��������,mΪ����,nΪ����
  mn <- dim(x)
  m <- mn[1]
  n <- mn[2]
  
  #����x�Ƿ�Ϊ����
  if(m != n) 
  {
    return ("Wrong dimensions of matrix!")
  }
  
  #����x�Ƿ�Ϊ�Գƾ���
  if(sum(t(x) != x) > 0) 
  {
    return ("Input matrix is not symmetrical!")
  }
  
  #LΪ��x��������ȵ���������ڴ�ŷֽ����������Ǿ���
  L <- matrix(0, m, m)
  
  #ѭ��ÿ����һ��,���һ�о���L��Ԫ��
  #����x��i�к͵�i��֮ǰ��Ԫ�ز���ʹ�ã��൱�ھ���x����һ��ά������������ѭ��������i�м�Ϊ��ǰ����x�;���L�ĵ�һ��
  for(i in 1:m)
  {
    #L�����Խ����ϵ�һ��Ԫ��Ϊx�����Խ����ϵ�һ��Ԫ�ؿ���
    L[i,i] <- sqrt(x[i,i])
    if(i < m)
    {
      #��ǰ����L�ĵ�һ�г���һ��Ԫ���������Ԫ��
      L[(i+1):m,i] <- x[(i+1):m,i]/L[i,i]
      
      #����L��һ�У�����һ��Ԫ�أ���������ת�õõ�TLM���ڸ��¾���x��Ч��ͬTLM%*%TLM
      TLV <- L[(i+1):m,i]                               #��¼�������һ�г���һ��Ԫ����ʣ��Ԫ��
      TLM <- matrix(TLV, m-i, m-i)                      #TLV���и��Ƴɾ���
      TLM <- sweep(TLM, 2, TLV, "*")                    #sweep(x�� MARGIN�� STATS�� FUN=��-���� ��) �Ծ����������
      #MARGINΪ1����ʾ�еķ����Ͻ������㣬Ϊ2��ʾ�еķ���������(��ָ���������еķ�������ȥ��)
      #STATS������Ĳ�����FUNΪ���㺯����Ĭ���Ǽ���
      
      #����һ��ά���ľ���x����Ϊԭ����Ӧλ���ϵ�Ԫ�ؼ�ȥTLM��Ϊ��һ��ѭ����׼��
      x[(i+1):m,(i+1):m] <- x[(i+1):m,(i+1):m] - TLM
    }
  }
  #����ķ���ֵΪ����Ҫ��������Ǿ���L
  L  
}

###mforwardsolve����������Է�����Lx=b������LΪ�����Ǿ���
#���룺�����Ǿ���L������b
#��������Է�����Ľ�x
mforwardsolve <- function(L, b)
{
  #��L��������,mΪL������,nΪL������
  mn <- dim(L)
  m <- mn[1]
  n <- mn[2]
  
  #�ж�L�Ƿ�Ϊ����
  if(m != n) 
  {
    return ("Wrong dimensions of matrix L!")
  }
  
  #�ж�L�Ƿ�Ϊ�����Ǿ���
  for (i in 1:(m-1))
  {
    if(sum(L[i,(i+1):m] != 0) > 0)#�����ж��������Ƿ�ȫΪ0Ԫ��
    {
      return ("Matrix L must be a lower triangular matrix!")
    }
  }
  
  #�ж�L��������b�ĳ����Ƿ����
  if(m != length(b))
  {
    return ("Wrong dimensions of matrix L or vector b!")
  }
  
  #0������¼�����
  x=rep(0, m)
  
  #ѭ��ÿ����һ��,���һ��x�е�Ԫ�أ���������L����x����b��ά����һ
  #������������L�ĵ�i�м�Ϊ��ǰ����L��һ�У�������x����b�ĵ�i��Ԫ�ؼ�Ϊ��ǰ������һ��Ԫ��
  for(i in 1:m)
  {
    #��ǰѭ����x�ĵ�һ��Ԫ��
    x[i] <- b[i] / L[i,i]
    #��ά���b����Ϊԭ��λ���ϵ�Ԫ�ؼ�ȥ��ǰ����L�ĵ�һ�еĳ˻�
    if(i < m) 
    {
      b[(i+1):m] <- b[(i+1):m] - x[i]*L[(i+1):m,i]
    }      
  }
  #�������ص�x������Ϊ���Է�����Ľ�
  x  
}

###mbacksolve����������Է�����Lx=b������LΪ�����Ǿ���
#���룺�����Ǿ���L������b
#��������Է�����Ľ�x
mbacksolve <- function(L, b)
{
  #��L��������,mΪL������,nΪL������
  mn <-dim(L)
  m <- mn[1]
  n <- mn[2]
  
  #�ж�L�Ƿ�Ϊ����
  if(m != n)
  {  
    return ("Wrong dimensions of matrix L!")
  }
  
  #�ж�L�Ƿ�Ϊ�����Ǿ���
  for (i in 2:m)
  {
    if(sum(L[i,1:(i-1)] != 0) > 0)
    {
      return ("Matrix L must be a upper triangular matrix!")
    }
  }
  
  #�ж�L��������b�������Ƿ����
  if(m != length(b))
  {
    return ("Wrong dimensions of matrix L or vector b!")
  }  
  
  x <- rep(0, m)
  
  #ѭ��ÿ����һ��,���һ��x�е�Ԫ�أ���������L����x����b��ά����һ
  #������������L�ĵ�i�м�Ϊ��ǰ����L���һ�У�������x����b�ĵ�i��Ԫ�ؼ�Ϊ��ǰ�������һ��Ԫ��
  for(i in m:1)
  {
    #��ǰѭ����x�����һ��Ԫ��
    x[i] <- b[i] / L[i,i]
    #��ά�������bΪԭ��λ���ϵ�Ԫ�ؼ�ȥ�ղ������xԪ���뵱ǰ�����Ǿ���L���һ�У������һ��Ԫ�أ��ĳ˻�
    if(i > 1) 
    {
      b[(i-1):1] <- b[(i-1):1] - x[i]*L[(i-1):1,i]
    }      
  }
  #��������ֵx������Ϊ���Է�����Ľ�
  x  
}

###ridgereg��������ʵ����ع����beta�Ĺ��ƣ�����x��y�ֱ�Ϊ�ع鷽�̵��Ա����������,lambdaΪL2������ĵ��ڲ���
#�˺���������Է�����(t(x)%*%x+lambada)%*%beta=t(x)%*%y,��t(x)%*%x+lambada����cholesky�ֽ�ΪR%*%t(R),forwardsolve���L%*%d=t(x)%*%y,����d=t(R)%*%beta,backsolve���t(R)%*%beta=d,���ò���beta�Ĺ���ֵ 
#���룺�Ա���x�������y�����ڲ���lambda
#������ع�ϵ��beta�Ĺ���ֵ
ridgereg <- function(lambda, x, y)
{
  #y=data[,m]; x=data[,-m]
  #nΪ�Ա�����������,��n������,pΪ�Ա�����������,��p������
  np <- dim(x)
  n <- np[1]
  p <- np[2]
  scale_data = scale_fun (x,y)
  #��x,y���б�׼��
  y = scale(y,center=T,scale=T) 
  x = scale(x,center=T,scale=T) 
  #����y��x�ľ�ֵ����
 
  #��׼��������û�нؾ���
  x = as.matrix(x)
  
  
  #����cholesky�ֽ���ȡ�ع鷽�̵Ĳ���beta�Ĺ���ֵ  
  V <- t(x)%*%x + diag( rep(lambda, p))              #t(x)%*%x+lambda��Ϊ���Է������ϵ������V
  U <- as.vector(t(x)%*%y)                           
  R <- mchol(V)                                           #����mchol������ϵ������V����cholesky�ֽ�,V=R%*%t(R)
  M <- mforwardsolve(R, U)                                #ʹ��ǰ�������R%*%M=t(x)%*%y,����M=t(R)%*%beta
  betai=mbacksolve(t(R), M)*()                                     #ʹ�ûش������t(R)%*%beta=M,���ɵ�beta�Ĺ���ֵ
  beta0=
}

scale_fun <- function(x, y)
{
  #����y��x�ľ�ֵ����
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
  
  #�����ݾ���nx�������У�ÿһ��Ϊһ����Ʒ�����ž����ԭ�������������е��õ�����ԭ�ṹΪdataframe
  nx <- as.matrix(nx, ncol <- p)
  n <- dim(nx)[1]
  
  #����Ԥ��ֵ
  apply(t(nx)*b[2:(p+1)], 2, sum) + b[1]  
}

###mridge��������ʵ��ɾȥһ����Ʒ����ع�
mridge=function(i,lambda,x,y) 
{
  ridgereg(lambda,x[-i,],y[-i])
}

###cvridgeregerr�����ý�����֤ʵ����ع飬��������Ϊ���ڲ���lambda,�Ա���x(���ݾ���),�����y,����ֵΪ���Ծ������   
#���룺������lambda,�Ա���x(���ݾ���),�����y
#�����ɾһ������֤��ع���Ծ������
cvridgeregerr<-function(lambda,x,y){  
  #lambda=1
  np<-dim(x)
  n<-np[1]
  p<-np[2]
  #�����е�Ԫ����Ϊ��һ����������mridge����ʾȥ�������ݱ�ţ������i��Ϊɾȥ��i����������ع�ϵ������ֵ
  coe<-t(apply(as.matrix(1:n,ncol=1),1,mridge,lambda,x,y))
  #����ÿ��ɾ�������ݵľ�ֵ����
  scale_data<-apply(as.matrix(1:n,ncol=1),1,scale_fun,x,y)
  
  #coe��i�к����ݾ����i����������Ե���ˣ�������ͣ�������Ծ������
  mean((apply(coe*x,1,sum)-y)^2)
}

###ridgeregerr�������ڼ���ѵ���������  
#���룺��ع鳬����lambda�����ݾ���x�������y
#�����ѵ���������
ridgeregerr=function(lambda,x,y)
{
  mean((pred(ridgereg(lambda,x,y),x)-y)^2)
}   

###############################
###�ڲ�ͬ��lambda��,�Ƚ�ѵ���������Ͳ��Ծ�������ѡȡ���ʵĵ��ڲ���lambda
library(ElemStatLearn)
x <- as.matrix(prostate[ ,1:8])
y <- as.vector(prostate[ ,9])
LAM <- seq(0.001, 10, len=50)
#������ع�50��ģ�͵�ѵ���������������listչ��������
err <- unlist(lapply(LAM, ridgeregerr, x, y))
#������ع�50��ģ�͵Ĳ��Ծ������������listչ��������
pe <- unlist(lapply(LAM, cvridgeregerr, x, y))
x <- rep(1:50, 2)
plot(pe)
#ȡ������֤��ʹ���Ծ��������С��lambda
lam=LAM[which.min(pe)]
xsd = apply(x,2,var)