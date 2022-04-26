import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

class diff:
    def __init__(self, a,J, N, deltax, deltat):
        self.a,self.N, self.J, self.deltat, self.deltax=a,N, J, deltat, deltax
        self.u=np.zeros([N,J])
        self.u[:,0]=100
        
        #for i in range(1,N-1):
        #    self.u[i,0]=np.sin(np.pi*i/(N-1))
        #    print(i, self.u[i,0])

        self.u[0,:]=0
        self.u[self.N-1,:]=0
        print(self.u)
        '''
 m is number of iterations
'''
    def solver1(self, m):
        u=self.u
        q=self.a**2 *self.deltat/(self.deltax)**2
        
        for i in range(m):
            
            for (n,j),x in np.ndenumerate(u):
                if j> self.J-2 or n == 0 or n > self.N-2:
                    continue

                else:
                    u[n,j+1]=u[n,j]+q*(u[n+1,j] + u[n-1,j]-2*u[n,j])


    def solver2(self,m):
        
        z=np.copy(self.u)
        q=self.a**2 *self.deltat/(self.deltax)**2
        
        
        a=np.repeat((-1),self.N-1)
        c=np.copy(a)
        b=np.repeat((2+2/q),self.N-1)
        d=np.zeros(self.N-1)
        y=np.copy(d)
        alpha=np.copy(d)
        beta=np.copy(d)

        
        for k in range(m):
            for n in range(1, self.J):
                u=self.u
                z=self.u
                d=np.zeros(self.N-1)
                y=np.copy(d)
                alpha=np.copy(d)
                beta=np.copy(d)
                for i in range(1, self.N-3, 1):
                    d[i]=u[i,n-1]+(-2+2/q)*u[i+1,n-1]+u[i+2,n-1]
                d[0]=u[0,n]+(-2+2/q)*u[1,n-1]+u[2,n-1]+u[0,n-1]
                d[self.N-3]=u[self.N-3,n-1]+(-2+2/q)*u[self.N-2,n-1]+u[self.N-1,n-1]+u[self.N-1,n]
                y[0]=y[0]+b[0]
                alpha[0]=-c[0]/y[0]
                beta[0]=d[0]/y[0]
                
                for i in range(1,self.N-3):
                    y[i]=b[i]+a[i]*alpha[i-1]
                    alpha[i]=-c[i]/y[i]
                    beta[i]=(d[i]-a[i]*beta[i-1])/y[i]
                y[self.N-3]=b[self.N-3]+a[self.N-3]*alpha[self.N-4]
                beta[self.N-3]=( d[self.N-3]-a[self.N-3]*beta[self.N-1-3])/y[self.N-3]

                u[self.N-2,n]=beta[self.N-3]
                for i in range(self.N-3, 0, -1):
                    u[i,n]=alpha[i-1]*z[i+1,n]+beta[i-1]
                
  

    def plotter(self):
        X = np.arange(0,self.J*self.deltax,self.deltax)
        T = np.arange(0,self.N*self.deltat,self.deltat)
        X, T = np.meshgrid(X, T)
        

        fig=plt.figure()
        ax=fig.add_subplot(projection='3d')
        surf = ax.plot_surface(X, T, self.u, cmap=cm.coolwarm, alpha = 0.5,
                       linewidth=0, antialiased=False)
        plt.xlabel('t')
        plt.ylabel('x')
        ax.set_title('T(x,t)')
        fig.colorbar(surf, ax=ax, fraction=0.02, pad=0.1, label='Temperature')
        plt.show()

        fig, ax = plt.subplots()
        CS = ax.contour(X, T, self.u)
        
        ax.clabel(CS, inline=True, fontsize=10)
        ax.set_title('equipotential lines  of temperature')
        fig.show()


x=diff(7, 100, 10, 100, 10)

x.solver2(1)
x.plotter()
