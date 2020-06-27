import math

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt


def GLf(f,a,b,n): #GLf: Gauss-Legendre quadrature for f
    """
    Compute numerical approximation using quadrature Gauss-Legendre.
    Weights and nodes are obtained with table for n=0,1,2,3,4
    Args:
        f (function): function expression of integrand
        a (float): left point of interval
        b (float): right point of interval
        n (float): number of subintervals
    Returns:
        sum_res (float): numerical approximation to integral of f in the interval a,b
    """   
    W_n0=[2]
    N_n0=[0]
    W_n1=[1,1]
    N_n1=[(1/3)**.5,-(1/3)**.5]
    W_n2=[5/9,8/9,5/9]
    N_n2=[-(3/5)**.5,0,(3/5)**.5]
    W_n3=[0.347855,0.652145,0.652145,0.347855]
    N_n3=[-0.861136,-0.339981,0.339981,0.861136]
    W_n4=[0.236927,0.478629,0.568889,0.478629,0.236927]
    N_n4=[-0.90618,-0.538469,0,0.538469,0.90618]
    W_vector=[]
    N_Vector=[]

    if n==0:
        W_Vector=W_n0
        N_Vector=N_n0
    elif n==1:
        W_Vector=W_n1
        N_Vector=N_n1
    elif n==2:
        W_Vector=W_n2
        N_Vector=N_n2
    elif n==3:
        W_Vector=W_n3
        N_Vector=N_n3
    elif n==4:
        W_Vector=W_n4
        N_Vector=N_n4
    else:
        W_Vector=W_n4
        N_Vector=N_n4

    part_sum=0
    const=(b-a)/2
    for i in np.arange(n+1):
        arg=.5*((b-a)*N_Vector[i]+a+b)
        part_sum+=W_Vector[i]*f(arg)
    sum_res=const*part_sum
    
    return sum_res

def relative_absolute_error(aprox, obj):
    """
    Add docstring
    Args:

    Return:

    """
    if(np.abs(obj) > 0 ):
        return np.abs(aprox-obj)/np.abs(obj)
    else:
        return np.abs(aprox-obj)


if __name__ == "__main__":
    #Gauss-Legendre
    f = lambda t: np.exp(-t**2/2)
    a = 0
    b = 1
    n = 4
    GL_approx = GLf(f,a,b,n)
    GL_obj, err = quad(f,a,b)
    #print('GL_approx',GL_approx)
    print('GL_obj:',GL_obj,'\n')
    n_Vector=[]
    rel_err_Vector=[]
    GL_aprox_Vector=[]
    rel_err_Vector_log=[]
    for i in np.arange(1,n+1):
      n_Vector.append(i)
      GL_aprox_Vector.append(GLf(f,a,b,i))
    
    GL_aprox_Vector = np.array(GL_aprox_Vector)
    rel_err_Vector = relative_absolute_error(GL_aprox_Vector, GL_obj)
    rel_err_Vector_log = np.log(rel_err_Vector)
    
    print('# of subintervals n:',n_Vector,'\n\n','GL approximation for each n',GL_aprox_Vector,'\n\n',\
    'Relative error for each n:',rel_err_Vector,'\n\n','Log10 Relative error for each n:',rel_err_Vector_log,'\n\n')

    #PLOT
    plt.subplot(2,1,1)
    plt.plot(n_Vector,rel_err_Vector,'o-')
    plt.xlabel('n subintervals');plt.ylabel('Relative_Error')
    plt.legend('Relative_error',loc=0)
    plt.grid(True)

    plt.subplot(2,1,2)
    plt.plot(n_Vector,rel_err_Vector_log,'s-')
    plt.xlabel('n subintervals');plt.ylabel('Log (Relative_Error)')
    plt.legend('Relative_error',loc=0)
    plt.grid(True)

    #plt.show()
    plt.savefig("25_june_class_feedback/Gauss_Legendre_Rel_Error_Natural_Log.png")

    