import numpy as np

class BSM:
    
    def __init__(self,K=True,r=True,q=True,option_type=True,sigma=True,S_max=True,T=True,delta_s=1,delta_t=1):
        self.K = K
        self.r = r
        self.q = q
        self.option_type = option_type
        self.sigma = sigma
        self.S_max = S_max
        self.T = T
        self.delta_s = delta_s
        self.delta_t = delta_t
        self.grid = np.zeros((int(S_max/delta_s)+1,int(T/delta_t)+1))
        
    def solve(self):
        #sweep the grid from right to left
        #initialize containers
        A = np.zeros((self.grid.shape[0],self.grid.shape[0]))
        #do the sweep
        for j in range(self.grid.shape[1]-1,-1,-1):
            if j == self.grid.shape[1]-1:
                for i in range(self.grid.shape[0]):
                    if i == 0:
                        self.grid[i,j] = self.get_BC('bottom',i)
                    elif i == self.grid.shape[0]-1:
                        self.grid[i,j] = self.get_BC('top',i)
                    else:
                        self.grid[i,j] = self.get_BC('right',i)
            else:
                b = self.grid[:,j+1]
                for i in range(np.size(b)):
                    if i == 0 or i == np.size(b)-1:
                        A[i,:] = 0
                        A[i,i] = 1
                    else:
                        A[i,i-1] = (0.5*(self.r - self.q)*self.delta_t*i) - (0.5*self.sigma**2*i**2*self.delta_t)
                        A[i,i] = 1 + (self.delta_t*self.r) + (self.sigma**2*i**2*self.delta_t)
                        A[i,i+1] = (-0.5*self.delta_t*i*(self.r - self.q)) - (0.5*self.sigma**2*i**2*self.delta_t)
                self.grid[:,j] = np.linalg.solve(A,b)
            for i in range(np.size(self.grid[:,j])):
                if self.grid[i,j] < 0:
                    self.grid[i,j] = 0
    
    def get_BC(self,loc,i):
        if loc == 'bottom':
            if self.option_type == 'call':
                return 0
            elif self.option_type == 'put':
                return self.K
        elif loc == 'top':
            if self.option_type == 'call':
                return self.S_max
            elif self.option_type == 'put':
                return 0
        elif loc == 'right':
            if self.option_type == 'call':
                return np.max((self.delta_s*i) - self.K,0)
            elif self.option_type == 'put':
                return np.max(self.K - (self.delta_s*i),0)
            
    def update_sigma(self,sigma):
        self.sigma = sigma