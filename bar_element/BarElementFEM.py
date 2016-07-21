import pylab as py # library for matrix operation

class BarElementFEM:
    
    def __init__(self,xy,connectivity,displacement): 
        '''
        member variables of the class BarElementFEM
        '''
        self._debug=False
        self._xy = py.array(xy) # copying xy values of the nodes
        self._connectivity = py.array(connectivity) # copying connectivity
        self._displacement = py.array(displacement)
        self.n_element = py.size(self._connectivity,0)
        self.n_nodes = py.shape(xy)[0] #number of node
        #global stiffness matrix intialized with zeros
        self.k_global = py.zeros((self.n_nodes*2,self.n_nodes*2))
        self.L = py.zeros(self.n_element)
        self.strain_global=[]
        
    def Transform2D(self,startXY,endXY):
        '''
        Function to obtain transformation matrix for an oriented bar element
        '''
        #find angle of orientation, atan2 takes care of the quadrant
        theta = py.math.atan2(endXY[1]-startXY[1],endXY[0]-startXY[0])
        if self._debug:        
            print 'theta', theta
        ct = py.math.cos(theta)#cosine 
        st = py.math.sin(theta)#sine
        T = py.array(([ct,0.],[st,0.],[0,ct],[0,st])) #transformation matrix
        return T
        
    def StiffnessAndStrain(self): 
        '''
        Function to obtain stiffness matrix and strains in each element
        '''
        for el in range(self.n_element): #enumerate through all the elements
            K = py.zeros((self.n_nodes*2,self.n_nodes*2))#Full stiffness matrix for one element
            if self._debug:            
                print 'bar element:', el+1
            startNode = self._connectivity[el,0]-1 #start and end node of the element
            endNode = self._connectivity[el,1]-1
            startXY = self._xy[startNode]        #start and end coordinates of the element
            endXY = self._xy[endNode]
            L = py.norm(endXY - startXY)    # Length of the element
            self.L[el]= L
            #new start and end coordinates of the element
            newstartXY=self._xy[startNode]+self._displacement[startNode] 
            newendXY=self._xy[endNode]+self._displacement[endNode] 
            newL=py.norm(newendXY - newstartXY) # Deformed Length of the element
            self.strain_global.append((newL-L)/L);     # Axial Strain in element
            T = py.zeros((4,2))
            T = self.Transform2D(startXY,endXY) #get the Transformation matrix for the element
            if self._debug: 
                print startNode, endNode, startXY, endXY
                print 'T'
                print T
            L_inv = 1.0/L
            k0=py.array([[L_inv, -L_inv],[-L_inv, L_inv]])#Stiffness matrix in 1D
            k_element = py.dot(T,py.dot(k0,py.transpose(T)) )#Stiffness matrix in 2D
            #Convert from 2x2 local K matrix to a full stiffness matrix but for single element 
            K[py.ix_([2*startNode,2*startNode+1,2*endNode,2*endNode+1],[2*startNode,2*startNode+1,2*endNode,2*endNode+1])] = K[py.ix_([2*startNode,2*startNode+1,2*endNode,2*endNode+1],[2*startNode,2*startNode+1,2*endNode,2*endNode+1])]+k_element 
            if self._debug:             
                print 'element stiffness', 'L_inv', L_inv
                print k_element
                print 'K'
                print K
            self.k_global = self.k_global + K  # Assemble all the stiffness into one Global Stiffness        
                    
        return self.k_global,self.strain_global

# in-class testing        
if __name__ == '__main__':
    
    f = open('case1.txt', 'r') #handle to file 'case1.txt' in read only mode
    i = 0
    for line in f:
        if line.strip():            
            t = line.split()
            #read number of nodes and elements
            if i==0:             
                n_node = int(t[0])
                n_elem = int(t[1])   
            #read node coordinate values
            elif i == 1:
                xy0 = [float(t[1]),float(t[2])]
                xy = [xy0]
            elif i <= n_node:
                xy.append([float(t[1]),float(t[2])])
            #read connectivity values            
            elif i == n_node+1:
                connection0 = [int(t[1]),int(t[2])]
                connectivity = [connection0]
            elif i <= n_node+n_elem:
                connectivity.append([int(t[1]),int(t[2])])    
            #read displacement values
            elif i == n_node+n_elem + 1:
                disp0 = [float(t[1]),float(t[2])]
                displacement = [disp0]
            elif i <= n_node+n_elem+n_node:
                displacement.append([float(t[1]),float(t[2])])
            i= i+1    
    '''
    xy =[[0.0,0.0],[1.0,0.0],[0.5,0.7]]
    connectivity = [[1,3],[3,2]]
    displacement = [[2.0,0.0],[2.0,-0.2],[2.0,1.4]]
    '''
    Truss = BarElementFEM(xy,connectivity,displacement)
    Truss.StiffnessAndStrain()

    print 'Nonzero entries in the global stiffness matrix'
    for i in range(Truss.n_nodes*2):
        for j in range(Truss.n_nodes*2):
            if Truss.k_global[i][j]!=0:            
                print i+1,j+1,Truss.k_global[i][j]
    print 'Axial strain in each element'
    for i in range(Truss.n_element):
        print i+1,Truss.strain_global[i]