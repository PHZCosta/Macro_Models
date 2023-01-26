import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


# Now, we want to implement these methods into a class, in order to make it easier to use:
class Ramsey:
    def __init__(self, A=1, alpha=0.3, delta=0.05, rho=0.02, theta=0.5, fine=60, far=1.5):
        self.A, self.alpha, self.delta, self.rho, self.theta = A, alpha, delta, rho, theta
        self.fine, self.far = fine, far
        
        # Finding the steady state:
        self.k_ss = (A*alpha/(delta+rho))**(1/(1-alpha))
        self.c_ss = A*self.k_ss**alpha - delta*self.k_ss
        
        # Defining the grid for consumption and capital:
        self.k_grid = np.linspace(0, self.k_ss*far, fine)
        self.c_grid = np.linspace(0, self.c_ss*far, fine)
        
        # Setting the nullcine equations:
        self.c_null = (A*alpha/(delta + rho))**(1/(1-alpha))
        self.k_null = A*self.k_grid**alpha - delta*self.k_grid
        
        # Setting the grids to be used in the phase diagram:
        self.M1, self.M2 = np.meshgrid(self.c_grid, self.k_grid)
        
        # Evaluating the differential equations in the grid:
        # Ps: This uses a function that I defined below.
        self.dc_grid, self.dk_grid = self.dc_dk(self.M1, self.M2)
        
       # Some arbitrary initial conditions to plot the explosive paths:
        self.arbitrary_init = np.array([(self.k_ss*0.1,self.c_ss*0.1),
                                (self.k_ss*0.2,self.c_ss*0.4),
                                (self.k_ss*1.3,self.c_ss*1.3),
                                (self.k_ss*1.4,self.c_ss*1.2)])
        
        
    # Defining a function that plots the nullcine equations:
    def simple_graph(self):
        plt.plot([self.c_null]*self.fine, self.c_grid, color = 'green', label='$\\dot c = 0$', linewidth=2)
        plt.plot(self.k_grid, self.k_null, color = 'red', label='$\\dot k = 0$', linewidth=2)
        plt.xlabel('Capital ($k$)')
        plt.ylabel('Consumption ($c$)')
        plt.ylim(0, self.c_grid[-1])
        plt.xlim(0, self.k_grid[-1])
        plt.legend()
        plt.grid()
        
    # Defining the differential equations:
    def dc_dk(self, c,k):
        # Ps: we have some zeros in the grid and this will cause Python to report a warning. 
        # To avoid that we will use the np.errstate function:
        with np.errstate(divide='ignore', invalid='ignore'):
            dc = c*(self.A*self.alpha*k**(self.alpha-1) - self.delta - self.rho)/self.theta
            dk = self.A*k**self.alpha - self.delta*k - c
        return dc, dk
    
    # Function to plot the phase diagram using the quiver function:
    def graph_quiver(self, n=3):
        # I load the simple graph function:
        self.simple_graph()
        # Because the grid is usually too dense, we only plot the n-th arrow:
        plt.quiver(self.M2[::n], self.M1[::n], self.dk_grid[::n], self.dc_grid[::n], 
                   angles = 'xy', scale_units = 'xy', scale = 2.5)
        
    # Function to plot the phase diagram using the streamplot function:
    def graph_streamplot(self):
        # I load the simple graph function:
        self.simple_graph()
        # Plotting the streamplot:
        plt.streamplot(self.k_grid, self.c_grid, np.transpose(self.dk_grid), np.transpose(self.dc_grid), density=1)
    
    # Function to plot the phase diagram with a saddle path:
    def graph_saddle_path(self):        
        # Plotting the streamplot:
        plt.streamplot(self.k_grid, self.c_grid, np.transpose(self.dk_grid), np.transpose(self.dc_grid), start_points=[ 
                                                                                                (self.k_ss,self.c_ss), 
                                                                                                (self.k_ss+1e-2,self.c_ss+1e-2)], integration_direction='backward',
                                                                                        broken_streamlines=False, color="orange", arrowsize=2, linewidth=2)
        # Plloting the old streamplot:
        plt.streamplot(self.k_grid, self.c_grid, np.transpose(self.dk_grid), np.transpose(self.dc_grid), density=0.5, arrowsize=1, linewidth=0.5)
        plt.plot([], [], '-', color="orange", label='Saddle Path')
         # I load the simple graph function:
        self.simple_graph()
    
    # Function to plot the phase diagram with some explosive paths:
    def graph_explosive_paths(self, init_cond):
        # Adding explosive paths:
        plt.streamplot(self.k_grid, self.c_grid, np.transpose(self.dk_grid), np.transpose(self.dc_grid), start_points=init_cond,
                       color='purple', integration_direction='forward', broken_streamlines=False, linewidth=2)
        # Also plotting the points of the initial conditions:
        plt.scatter(init_cond[:,0], init_cond[:,1], color='purple')
        # Plotting the saddle path:
        plt.streamplot(self.k_grid, self.c_grid, np.transpose(self.dk_grid), np.transpose(self.dc_grid), 
                       start_points=[ (self.k_ss,self.c_ss),(self.k_ss+1e-3,self.c_ss+1e-3)], 
                       integration_direction='backward',broken_streamlines=False, color="orange", arrowsize=2, linewidth=2)
        plt.plot([], [], '-', color="orange", label='Saddle Path')
        plt.plot([], [], '-', color="purple", label='Explosive Paths')
        self.simple_graph()
    
    
        

