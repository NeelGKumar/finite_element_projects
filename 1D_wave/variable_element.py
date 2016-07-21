import pylab as py # library for matrix operation

class Wave:
    
    def __init__(self, Nel, l, A, density, E, delT, nT, Freq, debug_status): 
        '''
        member variables of the class Wave
        '''
        
        self._debug = debug_status;
        #Problem conditions
        self.n_element = Nel; #Number of elements
        self.n_node = Nel+1; #Number of nodes
        self.L_element = py.zeros(self.n_element); #Length
        self.X = py.zeros(self.n_node) #node coordinates  
        self.LengthOfElements() #position of the nodes and element length
        
        self.area = A; #Crossection area
        self.density = density; #Density
        self.E = E; #Material constant
        self.DelT = delT; #Time step
        self.N_step = nT; #Number of time steps
        self.Ifreq= Freq; #Output Frequency         
              

        #find Lumped mass
        self.RM = py.zeros(self.n_node) # Lumped mass
        self.LumpedMass()
              
        #set boundary condtions- everything zero
        self.D = py.zeros(self.n_node) #Node displacement
        self.V = py.zeros(self.n_node) #Node Velocity
        self.F_ext = py.zeros(self.n_node) #External Force
        self.F_int = py.zeros(self.n_node) #Internal Force
        self.BoundaryConditions() #Is any vairable is NonZero
        
        #variables for storing the results        
        self.sigma = py.zeros(self.N_step)
        self.displacement = py.zeros(self.N_step)
        self.velocity = py.zeros(self.N_step)
        
        self.timesteps = py.zeros(self.N_step)#time at different steps for plotting

        if(self._debug):
            print 'Lumped mass matrix \n', self.RM, '\n'
            print 'Initial Nodal Displacement \n', self.D
            print 'Initial Nodal Velocity \n', self.V
         
        
        self.Solve()
        
    def LumpedMass(self):
        '''
        Form Lumped Mass Matrix- Only Diagonal
        '''
        
        
        
        for i in range(self.n_node-1):
            k = 0.5*self.density*self.area*self.L_element[i]
            self.RM[i] = self.RM[i] + k;
            self.RM[i+1] = self.RM[i+1] + k;
          
        

    def LengthOfElements(self):
        '''
        Variable Size
        '''
        
        for i in range(30):
            self.X[i+1] = self.X[i] + .5
            self.L_element[i] = .5
            
        for i in range(5):
            self.X[i+31] = self.X[i+30] + 1
            self.L_element[i+30] = 1
        print 'X of nodes \n', self.X, '\n'
        print 'Length of Elements \n', self.L_element, '\n'
        
           
           
    def BoundaryConditions(self):
        '''
        Set displacement, velocity and External Force Boundary Conditions
        '''
        self.F_ext[0] = 100.0;
            
       
    def InternalForce(self):
        
        #zero internal force vectors
        self.F_int = 0.0*self.F_int
        
        #Loop over Elements    
        for i in range(self.n_element):  
            RL = self.X[i+1] - self.X[i]
            
            #Compute Strain and Stress
            Strain = (self.D[i+1]-self.D[i])/RL
            Stress = self.E*Strain
            
            #Assemble Contribution into Iternal Force vector
            F1 = -Stress*self.area
            F2 = Stress*self.area
            self.F_int[i] = self.F_int[i] + F1
            self.F_int[i+1] = self.F_int[i+1] + F2
        
        #if(self._debug):
            #print 'F_int \n', self.RM, '\n' 
         
    
    
    def Solve(self):
        
        #Integration Loop
        for i in range(self.N_step):
            
            #Get internal Forces from last time step
            self.InternalForce()
            
            #Update displacements using central difference method
            for j in range(self.n_node):
                D_old = self.D[j]
                self.D[j] = self.DelT**2*(self.F_ext[j] - self.F_int[j])/self.RM[j] + self.D[j] + self.DelT*self.V[j]
                
                #Enforce zero displacement boundary condition at buil-in end                
                if ((j+1) == self.n_node):
                    self.D[j] = 0.0;   
                self.V[j] = (self.D[j] - D_old)/self.DelT
            
            #print 'Displacement \n', self.D
            self.sigma[i] = self.E*(self.D[21] - self.D[20])/self.L_element[20]
            self.displacement[i] = self.D[20] + 0.5*(self.D[21] - self.D[20])
            self.velocity[i] = self.V[20] + 0.5*(self.V[21] - self.V[20])
            
            self.timesteps[i] = i*self.DelT*1.0e3
        


# in-class testing        
if __name__ == '__main__':
    Nel = 35; l = .5; A = 1.0; density = 7.4e-4; 
    E = 30e6; delT = 2.4e-6; nT = 400; Freq = 1.0;
    debug_status = True
    
    W = Wave( Nel, l, A, density, E, delT, nT, Freq, debug_status)
    #print 'stress \n', W.sigma
    py.figure()    
    py.plot(W.timesteps, W.sigma)
    s = 'Stress in psi at x=9.75 '
    py.title(s)
    py.xlabel('Time(ms)')
    
    py.figure()    
    py.plot(W.timesteps, W.displacement)
    s = 'Dsiplacement at x=9.75'
    py.title(s)
    py.xlabel('Time(ms)')
    
    py.figure()    
    py.plot(W.timesteps, W.velocity)
    s = 'Velocity at x=9.75'
    py.title(s)
    py.xlabel('Time(ms)')
    