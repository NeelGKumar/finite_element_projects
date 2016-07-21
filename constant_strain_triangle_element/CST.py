import pylab as py # library for matrix operation

class CST:
    
    def __init__(self,E,v,iflag,xy,element,displacement, debug_status): 
        '''
        member variables of the class BarElementFEM
        '''
        self.debug = debug_status #ture if in debugging mode
        self.E = E #Youngs modulus
        self.v = v #Poisson Ratio
        self.iflag = iflag #Plane stress-strain identifier
        self.xy = py.array(xy) # copying xy values of the sorted nodes
        self.element = py.array(element) # copying connectivity
        self.displacement = py.array(displacement) # copying displacement
        self.n_element = py.size(self.element,0) #number of elements
        self.n_nodes = py.shape(xy)[0] #number of nodes      
        self.k_global = py.zeros((self.n_nodes*2,self.n_nodes*2)) #global stiffness matrix
        self.E_matrix = self.Constitutive_matrix() #Stores constitutive matrix
        
    def Constitutive_matrix(self):
        '''
        Generates E matrix for plane stress or plane strain
        '''
        v = self.v
        if self.iflag == 1: #plane strain
            return self.E/(1-2*v)/(1+v)*py.array([[1-v,v,0],[v,1-v,0],[0,0,0.5-v]])
        elif self.iflag == 0: #plane stress
            return self.E/(1-v*v)*py.array([[1,v,0],[v,1,0],[0,0,(1-v)/2]])
 
    def B(self,el):
        '''
        Generates Strain Matrix for a Constant Strain Triangle element, Assuming first node lies at origin.
        '''
        
        #Relative Element x,y- Coordinates with first node at origin
        x2 = self.xy[el[1]-1][0]-self.xy[el[0]-1][0]
        y2 = self.xy[el[1]-1][1]-self.xy[el[0]-1][1]
        x3 = self.xy[el[2]-1][0]-self.xy[el[0]-1][0]
        y3 = self.xy[el[2]-1][1]-self.xy[el[0]-1][1]
        
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
            node1x = 2*(node1-1);node2x = 2*(node2-1);node3x = 2*(node3-1);
            node1y = 2*(node1-1)+1;node2y = 2*(node2-1)+1;node3y = 2*(node3-1)+1;
            
            #Area, Strain Matrix and E Matrix multiplied to get element stiffness            
            [J,B] = self.B(el)
            local_k =0.5*abs(J)*py.dot(py.transpose(B),py.dot(self.E_matrix,B))
            if self.debug:            
                print 'K for elem', el, '\n', local_k
                
            #Element K-Matrix converted into Global K-Matrix format 
            K[py.ix_([node1x,node1y,node2x,node2y,node3x,node3y],[node1x,node1y,node2x,node2y,node3x,node3y])] = K[py.ix_([node1x,node1y,node2x,node2y,node3x,node3y],[node1x,node1y,node2x,node2y,node3x,node3y])]+local_k

            #Adding contibution into Global Stiffness           
            self.k_global = self.k_global + K
            
        if self.debug:            
                print 'Global Stiffness','\n', self.k_global        
        
       
    def Strain_stress(self,el):
        
        #Displacement formatted for an element
        d=py.array([self.displacement[el[0]-1][0],self.displacement[el[0]-1][1],self.displacement[el[1]-1][0],self.displacement[el[1]-1][1],self.displacement[el[2]-1][0],self.displacement[el[2]-1][1]])
        
        #Element Strain vector = Product of Strain Matrix and Displacement
        [J,B] = self.B(el)        
        strain = py.dot(B,d.T)
        
        #Element Stress vector = Product of Element K-Matrix and Strain Vector
        stress = py.dot(self.E_matrix,strain)
        return strain, stress
        
    def Global_load(self):
        #Format displacement into a vector
        disp = self.displacement.reshape(2*self.n_nodes)
        
        #Global Load = Product of global stifness matrix and displacement vector
        return py.dot(self.k_global,disp)
             
      
if __name__ == '__main__':  
    
    #Control Debugging mode
    debug_status = False 
    if debug_status:
        print '\n', '----In DEBUG MODE----', '\n'
    
    f = open('case1.txt', 'r') #handle to file 'case1.txt' in read only mode
    i = 0
    for line in f:
        if line.strip():            
            t = line.split()
            
            #read number of nodes, elements and flag for plain stress/strain
            if i==0:             
                n_node = int(t[0]); n_elem = int(t[1]); 
                xy = py.zeros((n_node,2))
                disp = py.zeros((n_node,2))
                element = py.zeros((n_elem,3), int)
                    
            #read Youngs Modulus, Poisson's ratio and plane stress/strain
            elif i == 1:
                YoungM = float(t[0])
                Poisson = float(t[1])
                Plane = float(t[2])
                        
            #read node coordinate values
            elif i <= n_node+1:               
                xy[int(t[0])-1] = [float(t[1]),float(t[2])]
                           
            #read connectivity values            
            elif i <= n_node+n_elem+1:
                element[i - n_node-2] = [int(t[0]),int(t[1]),int(t[2])]

            #read displacement values
            elif i <= n_node+n_elem+n_node+1:
                disp[int(t[0])-1] = [float(t[1]),float(t[2])]
            i= i+1 
            
    if debug_status:
        print 'xy', xy, '\n', 'element', element, '\n', 'disp', disp 
    
    #Create a new object "testcase" for CST class   
    testcase=CST(YoungM,Poisson,Plane,xy,element,disp, debug_status)    

    #Print non-zero entities in Global Stiffness Matrix
    testcase.Global_Stiffness()    
    print '\n','Nonzero entries in the global stiffness matrix'
    for i in range(testcase.n_nodes*2):
        for j in range(testcase.n_nodes*2):
            if abs(testcase.k_global[i][j])> 1e-10 :            
                print '(',i+1,',',j+1,',',testcase.k_global[i][j],')'
    
    #Print non-zero entities in Global Stiffness Matrix    
    print '\n','Stress and Strain in each element'    
    for i in range(testcase.n_element):
        [strain,stress] = testcase.Strain_stress(testcase.element[i])
        print 'element ',i+1, ': strain', strain,'; stress ', stress 
    
    #Print Global Load Vector
    global_load = testcase.Global_load()
    print '\n', 'Global Load Vector', '\n', global_load