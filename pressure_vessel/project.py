import pylab as py # library for matrix operation
import matplotlib.pyplot as plt

class project:

    def __init__(self,E,v,N,a,b,pressure, debug_status,error_analysis): 
        '''
        member variables of the class Project
        '''
        self.debug = debug_status #ture if in debugging mode
        self.E = E #Youngs modulus
        self.v = v #Poisson Ratio
        self.N=N #Mesh resolution
        self.a=a #Inner radius
        self.b=b #Outer radius
        self.P = pressure #Pressure       
        self.n_nodes = (N+1)**2  #No. of nodes      
        self.xy = py.zeros(((N+1)**2,2)) #Coordinates of nodes
        self.n_el = self.N**2*2 #No. of elements    
        self.element = py.zeros((self.n_el,3), int) #Connectivity information, nodes indices for each elements
        self.load = py.zeros(2*self.n_nodes)    #load vector
        self.k_global = py.zeros((self.n_nodes*2,self.n_nodes*2)) #global stiffness matrix for plane stress condition
        self.E_matrix = self.E/(1-v*v)*py.array([[1,v,0],[v,1,0],[0,0,(1-v)/2]]) #Stores constitutive matrix
        self.displacement = py.zeros(2*self.n_nodes**2) #displacement vector
        self.displ_scaler = py.zeros(self.n_nodes**2)   #radical displacment, scaler
        self.xy_final = py.zeros(2*self.n_nodes**2) #nodes coordinates after displacement
        self.stress_vector = py.zeros((len(self.element),3))    #stress vector from simulation
        self.strain_vector = py.zeros((len(self.element),3))    #strain vector from simulation
        self.von = py.zeros(len(self.element)) #Von Mises stress
        self.cond=0
        self.residual=0

        '''
        member functions intialization
        '''
        self.Mesh() #Mesh generator, outout: node coordinates and connectivity information
        self.LoadVector()   #Calculate nodal loads from distributive loads
        self.Global_Stiffness() #Calculate element stiffness matrix and assemble them together
        self.Dirichlet_BC() #Add boundary condition, to both global stiffness matrix and loads
        self.Global_displacement()  #Solve displacement by linear solver
        self.Strain_stress()    #Calculate the strain, stress and Von Mises stress from displacement
        
        if error_analysis: #Error analysis by comparing to analytical solution
            [self.error_displ, self.error_se,self.rerror_displ,self.rerror_se, self.error_normal_avg,self.error_normal,self.rerror_normal_avg,self.rerror_normal]=self.Error() 
            #Displacement error, Relative displacemnt error, Strain energy error, relative strain energy error, 
            #normal stress error over whole mesh, normal stress error,  relative normal stress error over whole mesh, relative normal stress error 

        if error_analysis==False:   
            self.Plot_lines()   #plot all the results including:
                                            #(1) von Mises stress, 
                                            #(2) mesh and displacement, 
                                            #(3) normal stress error, 
                                            #(4) nodal displacement error
      
    def Mesh(self):
        '''
        Generates mesh 
        output: 
            node coodinates: xy
            connectivity information: element
        '''
        # generate node coodinates: x and y
        xy = py.zeros(((self.N+1)**2,2))
        for k in range (0,self.N+1):
            for j in range (0,self.N+1):
                #equations for x and y coordinate from handout
                x=(self.a+(self.b-self.a)/self.N*k)*py.cos(py.pi*j/2/self.N)
                y=(self.a+(self.b-self.a)/self.N*k)*py.sin(py.pi*j/2/self.N)
                xy[j+k*(self.N+1)]=[float(x),float(y)]        
        self.xy = xy
        
        # connectivity information: element       
        element = py.zeros((self.N**2*2,3), int)   
        m=0
        for i in range (0,(self.N+1)*self.N):
            if (i+1)%(self.N+1) !=0:
               #define elements, two at a time 
               element[m]=[i,i+self.N+1,i+1] 
               element[m+1]=[i+self.N+1,i+self.N+2,i+1]
               m+=2
        self.element = element
        
        if self.debug:
            print '\n', 'in mesh: nodes','\n', xy
            print 'elements', '\n',element
    
    def LoadVector(self):
        '''
        Generates work equivalent load vector for thick walled cylindrical pressure vessel.
        The pressure is uniform inside the cylinder in the radial direction. Also we found theoratically that it would be distributed equally among two inside nodes for every element.                                                                                      
        '''
        #Node Coordinates and Elements copied
        xy = self.xy
        el_vector = self.element
        
        #initialize a load vector
        load = py.zeros(2*self.n_nodes)
        
        #Variables to check total X and Y load         
        PX = 0; PY = 0;
        
        #Load calculation by iterating over the inner surface elements
        for i in range(self.N):
            el = el_vector[2*i,:] #element with pressure
            side = xy[el[0]]-xy[el[2]] #inside edge in local coordinates 
            l = py.norm(side) #length of the edge
            
            #Load in X and Y direction. Halved since equally shared by two vertices
            Px = 0.5*self.P*abs(side[1]) 
            Py = 0.5*self.P*abs(side[0])
            
            #Load Componenets added in the global load vector
            load[2*i] = load[2*i] + Px;
            load[2*i+2] = load[2*i+2] + Px;
            load[2*i+1] = load[2*i+1] + Py;
            load[2*i+3] = load[2*i+3] + Py;
            
            if self.debug:
                print 'element side facing inside', el[0],el[2]
                print 'length', l
                print 'Px', Px, 'Py', Py   
            #total x and y load components
            PX += 2*Px; PY += 2*Py;
            
        if self.debug:
            print 'PX', PX, 'PY', PY, 'Total', py.math.sqrt(PX*PX+PY*PY)
            print 'load' , load
            
        self.load = load # copied to class variable
        
    def B(self,el):
        '''
        Generates Strain Matrix for a Constant Strain Triangle element, Assuming first node lies at origin.
        '''
        
        #Relative Element x,y- Coordinates with first node at origin
        x2 = self.xy[el[1]][0]-self.xy[el[0]][0]
        y2 = self.xy[el[1]][1]-self.xy[el[0]][1]
        x3 = self.xy[el[2]][0]-self.xy[el[0]][0]
        y3 = self.xy[el[2]][1]-self.xy[el[0]][1]
        
        #2*Area        
        J =x2*y3-x3*y2
        
        #Strain Matrix directly from formula
        B = py.array([[(y2-y3)/J, 0, y3/J, 0, -y2/J, 0], [0, (x3-x2)/J, 0, -x3/J, 0, x2/J], [(x3-x2)/J, (y2-y3)/J, -x3/J, y3/J, x2/J, -y2/J]])  
        
        return J,B
  
    
    def Global_Stiffness(self):
        '''
        Generates Global Stiffness Matrix for the plane structure
        '''
        elem = self.element;
        B = py.zeros((6,6))
        for i in range (0,py.size(elem,0)): 
            #for each element find the stifness matrix
            K = py.zeros((self.n_nodes*2,self.n_nodes*2))            
            el = elem[i]
            
            #nodes formatted for input            
            [node1, node2, node3] = el;
            node1x = 2*(node1);node2x = 2*(node2);node3x = 2*(node3);
            node1y = 2*(node1)+1;node2y = 2*(node2)+1;node3y = 2*(node3)+1;
            #Area, Strain Matrix and E Matrix multiplied to get element stiffness            
            [J,B] = self.B(el)
            local_k =0.5*abs(J)*py.dot(py.transpose(B),py.dot(self.E_matrix,B))
            
            if self.debug:            
                print 'K for elem', el, '\n', local_k
            #Element K-Matrix converted into Global K-Matrix format 
            K[py.ix_([node1x,node1y,node2x,node2y,node3x,node3y],[node1x,node1y,node2x,node2y,node3x,node3y])] += local_k

            #Adding contibution into Global Stiffness           
            self.k_global += K
            
        if self.debug: 
                print 'Global Stiffness','\n', self.k_global, '\n', 'size', py.shape(self.k_global), '\n Symmetry test' , py.dot(py.inv(self.k_global),self.k_global)    
            
    def Dirichlet_BC(self):
        '''
        Apply grounding conditions for i-th DOF by making K[i,:]and K[:,i] zero and K[i,i] = 1.
        Also i-th load is made zero.
        '''
        
        # Copy Global Stiffness, Global Load and number of nodes along each polar coordinate
        K = self.k_global
        l = self.load
        n = self.N+1
        
        if self.debug:
            print '\n','\n','Dirichlet_BC', K, '\n','\n','\n load', l
        #iteratively apply boundary conditions on all the grounded degrees of freedom 
        for i in range(0,n):
            if self.debug:
                print '\n every iteration' , 2*n*i+1, 2*n*(i+1)-2
            #Y_displacement = 0 for nodes on X axis
            K[2*n*i+1,:] = py.zeros((1,self.n_nodes*2))
            K[:,2*n*i+1] = py.zeros((1,self.n_nodes*2))
            K[2*n*i+1,2*n*i+1] = 1.0
            l[2*n*i+1] = 0.0
            
            #X_displacement = 0 for nodes on Y axis
            K[2*n*(i+1)-2,:] = py.zeros((1,self.n_nodes*2))
            K[:,2*n*(i+1)-2] = py.zeros((1,self.n_nodes*2))
            K[2*n*(i+1)-2,2*n*(i+1)-2] = 1.0
            l[2*n*(i+1)-2] = 0.0
            
        self.k_global = K # copy new global
        self.load = l # copy new load
        if self.debug:
            print '\n','Dirichlet_BC','\n', K, '\n','\n','\n load', l
        
       
    def Strain_stress(self):
        '''
        Strain is obtained by Strain_Matrix (B)X Displacement (u). Stress is Stiffess_Matrix(E) x Strain.
        Von mises stress is also computed to study convergence.
        '''
        
        # copy displacement vector after reshaping with two columns for x and y         
        d = self.displacement.reshape(self.n_nodes,2) 
        
        #stress-strain calculation for each element
        for i in range(self.n_el): 
            
            el = self.element[i] # present element
            
            #Displacement formatted for an element
            disp=py.array([d[el[0]][0],d[el[0]][1],d[el[1]][0],d[el[1]][1],d[el[2]][0],d[el[2]][1]])
        
            #Element Strain vector = Product of Strain Matrix and Displacement
            [J,B] = self.B(el)
            strain = py.dot(B,disp.T)
            self.strain_vector[i] = strain
        
            #Element Stress vector = Product of Element K-Matrix and Strain Vector
            stress = py.dot(self.E_matrix,strain)
            self.stress_vector[i] = stress

            #von-mises stress for plotting
            self.von[i] = py.math.sqrt(0.5*((stress[0]-stress[1])**2 + stress[0]**2 + stress[1]**2 + 6*(stress[2])**2))
              
               
    def Global_displacement(self):
        '''
        Final displacement as sum of intitial nodal values and caluclated displacement
        '''
        
        #Solve linear system for displacement and copy to class variable
        e=py.linalg.eigvals(self.k_global)
        self.cond=e[0]/e[-1]  #conditioning number
        disp = py.linalg.solve(self.k_global,self.load) #solving linear equation
        self.residual=max(self.load-py.dot(self.k_global,disp)) #residual calculation
        self.displacement = py.array([disp]).T
        displ_scaler = py.zeros((len(disp)/2,1))
        
        #calculating distance moved by each node
        for i in range (len(disp)/2):
            displ_scaler[i] = (disp[2*i]**2+disp[2*i+1]**2)**0.5
        
        self.displ_scaler=displ_scaler   #copy displacement values  
        
        #final nodal position
        self.xy_final = self.xy + disp.reshape(self.n_nodes,2)
        
        return self.xy_final
    
    def Analytical_displacement(self):
        #Analytical_displacement
        r=py.array([(self.xy[:,0]**2+self.xy[:,1]**2)**0.5]).T
        return self.P*(1-self.v**2)*self.a**2/self.E/(self.b**2-self.a**2)*(r/(1+self.v)+self.b**2/(1-self.v)/r) 
    
    def Analytical_Strain_Energy(self):
        #Analytical_Strain_Energy
        return py.pi*(self.P**2)*(self.a**2)/4/(self.b**2-self.a**2)/self.E*((1-self.v)*self.a**2+(1+self.v)*self.b**2)
        
    def Analytical_Normal_Stress(self):
        #Analytical_Normal_Stress
        n_stress = 2*(self.a**2)/(self.b**2-self.a**2)*self.P
        normal_vec = n_stress*py.ones(self.n_el)
        return normal_vec
     
    def Error(self):
        '''
        Displacement error, strain energy error and normal stress error compare to Analytical Solution
        '''    
        #dispalcement error and relative error
        error_displ = self.Analytical_displacement()-self.displ_scaler
        rerror_displ = error_displ/self.Analytical_displacement()
        
        #strain energy error and relative error
        error_se = self.Analytical_Strain_Energy()-py.dot(py.dot(self.displacement.T,self.k_global),self.displacement)/2.0        
        rerror_se = error_se/self.Analytical_Strain_Energy()
        
        #normal stress error - average, relative         
        error_normal_stress_avg = sum(self.Analytical_Normal_Stress() - self.stress_vector[:,0]-self.stress_vector[:,1])/self.n_el
        error_normal = self.stress_vector[:,0]+self.stress_vector[:,1]-self.Analytical_Normal_Stress()
        rerror_normal_stress_avg = error_normal_stress_avg /self.Analytical_Normal_Stress()[0]
        rerror_normal = error_normal/self.Analytical_Normal_Stress()
        
        return error_displ, error_se,rerror_displ, rerror_se, error_normal_stress_avg,error_normal,rerror_normal_stress_avg,rerror_normal 
        
        
    def Plot_lines(self):
        xy = self.xy
        
        #(1) plot von Mises stress, 
        plt.figure(1)
        plt.grid(True)
        plt.gca().set_aspect('equal')
        d = self.xy
        m = min(self.von); M = max(self.von); l = M - m;

        for i in range(self.n_el):
            k =  0.5*(self.von[i]-m)/l
            el = self.element[i]
            x=py.array([d[el[0]][0],d[el[1]][0],d[el[2]][0]])
            y=py.array([d[el[0]][1],d[el[1]][1],d[el[2]][1]])
            plt.fill(x,y,color=(0,.25+k,.75-k))       
        s = 'plot of von mises stress: Max ' + str(M) + ', Min ' + str(m)
        plt.title( s)
        plt.xlabel('X - axis')
        plt.ylabel('Y - axis')
        
        #(2)plot mesh and displacement,
        plt.figure(2)        
        for el in self.element:
                vx=[xy[el[0],0],xy[el[1],0]]
                vy=[xy[el[0],1],xy[el[1],1]]
                plt.plot(vx,vy,'k') 
                vx=[xy[el[2],0],xy[el[1],0]]
                vy=[xy[el[2],1],xy[el[1],1]]
                plt.plot(vx,vy,'k') 
                vx=[xy[el[0],0],xy[el[2],0]]
                vy=[xy[el[0],1],xy[el[2],1]]
                plt.plot(vx,vy,'k')
        
        xy = self.xy_final
        for el in self.element:
                vx=[xy[el[0],0],xy[el[1],0]]
                vy=[xy[el[0],1],xy[el[1],1]]
                plt.plot(vx,vy,'c') 
                vx=[xy[el[2],0],xy[el[1],0]]
                vy=[xy[el[2],1],xy[el[1],1]]
                plt.plot(vx,vy,'c') 
                vx=[xy[el[0],0],xy[el[2],0]]
                vy=[xy[el[0],1],xy[el[2],1]]
                plt.plot(vx,vy,'c')
        plt.grid(True)
        plt.gca().set_aspect('equal')
        stress_err = self.Analytical_Normal_Stress() - self.stress_vector[:,0]-self.stress_vector[:,1]
        
        if error_analysis:
            #(3) normal stress error,
            plt.figure(3)
            plt.grid(True)
            plt.gca().set_aspect('equal')
            M = max(stress_err); m = min(stress_err); 
            for i in range(self.n_el):                
                k =  stress_err[i]
                el = self.element[i]
                x=py.array([d[el[0]][0],d[el[1]][0],d[el[2]][0]])
                y=py.array([d[el[0]][1],d[el[1]][1],d[el[2]][1]])
                
                if(k>0):
                    plt.fill(x,y,color=(0,k/M,0))
                elif(k<0):
                    plt.fill(x,y,color=(0,0,k/m))
                
            s = 'plot of error in normal stress: Max ' + str(M) + ', Min ' + str(m)
            plt.title( s)
            plt.xlabel('X - axis')
            plt.ylabel('Y - axis')

        
if __name__ == '__main__':  
    
    #Control Debugging mode
    debug_status = False
    error_analysis = True #True to turn on error analysis
    if debug_status:
        print '\n', '----In DEBUG MODE----', '\n'
        
    #program input
    E=10.0    #Young`s Modules   
    v=0.25 #possison ratio
    a=0.5   #inner radius
    b=1.0   #outter radius
    P=1.0 #pressure     
    N=5     #resolution
     
    if error_analysis==False:
        testcase=project(E,v,N,a,b,P,debug_status,error_analysis)
                
   
   #for error analysis purpose only
    if error_analysis:
        startN=2 #resolution, multiple N is used for convergence analysis
        endN=12 
        # Arrays for convergence analysis
        numoftest=endN-startN    
        numofels=py.zeros((numoftest,1)) 
        
        maxvonarray=py.zeros((numoftest,1))   
        minvonarray=py.zeros((numoftest,1))  
        
        Maxdisplacementerrorarray=py.zeros((numoftest,1)) 
        Maxdisplacementrerrorarray=py.zeros((numoftest,1)) 
        
        Strainenergyerrorarray=py.zeros((numoftest,1)) 
        Strainenergyrerrorarray=py.zeros((numoftest,1)) 
        
        Normalstresserrorarray=py.zeros((numoftest,1)) 
        Normalstressrerrorarray=py.zeros((numoftest,1)) 
        
        avgNormalstresserrorarray=py.zeros((numoftest,1)) 
        
        Condarray=py.zeros((numoftest,1)) 
        residualarray=py.zeros((numoftest,1)) 
        
        # Run simulation for multiple resolution
        for i in range(startN,endN):  
            testcase=project(E,v,i,a,b,P,debug_status,error_analysis)
            numofels[i-startN]=2*(i**2)   
            
            maxvonarray[i-startN]=max(testcase.von)   
            minvonarray[i-startN]=min(testcase.von) 
            
            Maxdisplacementerrorarray[i-startN]=max(abs(testcase.error_displ))
            Maxdisplacementrerrorarray[i-startN]=max(abs(testcase.rerror_displ))
            
            Strainenergyerrorarray[i-startN]=testcase.error_se
            Strainenergyrerrorarray[i-startN]=testcase.rerror_se
            
            Normalstresserrorarray[i-startN]=max(testcase.error_normal)
            Normalstressrerrorarray[i-startN]=max(testcase.rerror_normal)
            
            avgNormalstresserrorarray[i-startN]=testcase.error_normal_avg
            
            Condarray[i-startN]=testcase.cond
            residualarray[i-startN]=testcase.residual
            # Print and plot result for highest resolution
            if i==1+numoftest:
                testcase.Plot_lines()
                print '\nvon Mises stress \n max:', max(testcase.von),'min: ',min(testcase.von)
                print '\nMax displacement error and relative displacement error: \n', max(abs(testcase.error_displ)),max(abs(testcase.rerror_displ))
                print '\nStrain energy error and relative strain energy error:\n', testcase.error_se, testcase.rerror_se
                print '\nMax Normal stress Error and max relative error in normal stress\n', max(testcase.error_normal), max(testcase.rerror_normal)
                print '\nAverage Normal stress Error\n', testcase.error_normal_avg 
                
            
        #plot convergence result
        py.figure()    
        py.plot(numofels,maxvonarray)
        s = 'Max von Mises stress'
        py.title(s)
        py.xlabel('No. of element')
        py.figure()    
        py.plot(numofels,minvonarray)
        s = 'Min von Mises stress'
        py.title(s)
        py.xlabel('No. of element')
        
        
        py.figure()    
        py.plot(numofels,Maxdisplacementerrorarray)
        s = 'Max displacement error'
        py.title(s)
        py.xlabel('No. of element')
        py.figure()    
        py.plot(numofels,Maxdisplacementrerrorarray)
        s = 'Max relative displacement error'
        py.title(s)
        py.xlabel('No. of element')
        
        py.figure()  
        py.plot(numofels,Strainenergyerrorarray)
        s = 'Strain energy error'
        py.title(s)
        py.xlabel('No. of element')
        py.figure()  
        py.plot(numofels,Strainenergyrerrorarray)
        s = 'Relative strain energy error'
        py.title(s)
        py.xlabel('No. of element')
        
        py.figure()  
        py.plot(numofels,Normalstresserrorarray)
        s = 'Max Normal Stress Error'
        py.title(s)
        py.xlabel('No. of element')
        py.figure()  
        py.plot(numofels,Normalstressrerrorarray)
        s = 'Max relative Normal Stress Error'
        py.title(s)
        py.xlabel('No. of element')
        
        py.figure()  
        py.plot(numofels,avgNormalstresserrorarray)
        s = 'Average Normal Stress Error'
        py.title(s)
        py.xlabel('No. of element')
        
        py.figure()  
        py.plot(numofels,Condarray)
        s = 'Condition number for stiffness matrix'
        py.title(s)
        py.xlabel('No. of element')
        
        py.figure()  
        py.plot(numofels,residualarray)
        s = 'Residual'
        py.title(s)
        py.xlabel('No. of element')


    

    
    
    
               
                
                
                
            
