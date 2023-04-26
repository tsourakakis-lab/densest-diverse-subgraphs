class FibonacciHeap:
    def __init__( self ):
        self.min_node = None
        self.n = 0
        self.root_list = None
    
    # internal class -- Node of a tree
    class Node:
        def __init__(self,key, value = None):
            self.key = key
            self.value = value
            self.parent = self.child = self.left = self.right = None
            self.degree = 0
            self.mark = False
            
    def insert( self, key, value ):
        x = self.Node( key, value ) 
        x.left = x.right = x
        self.insert_root_list( x )
        if not self.min_node: # if list does not have a min node
            self.min_node = x
        else:
            # insert x into  H's root list
            if x.key < self.min_node.key:
                self.min_node = x
        self.n += 1
        return x
            
    # insert node to root's doubly linked list        
    def insert_root_list( self, x ):
        if self.root_list == None:
            self.root_list = x
        else:
            x.left = self.root_list
            x.right = self.root_list.right
            self.root_list.right.left = x
            self.root_list.right = x
        
    def find_min( self ):
        return self.min_node
    
    def unite( self, h2 ):
        # to implement
        return None
    
    # traverse doubly linked list 
    def traverse( self, head ):
        node = stop = head
        flag = True
        while True:
            if node == stop and flag == False:
                break
            elif node == stop:
                flag = False
            yield node
            node = node.right
    
    def remove_from_root_list( self, x ):            
        if x == self.root_list:
            self.root_list = x.right
        x.left.right = x.right
        x.right.left = x.left
        return
    
    def extract_min( self ):
        z = self.find_min()
        if z:
            if z.child:
                children = [x for x in self.traverse( z.child )]
                for x in children:
                    self.insert_root_list( x )
                    x.parent = None
            # remove z from the root list of H
            self.remove_from_root_list( z )
            if z == z.right:
                self.min_node = self.root_list = None
            else:
                self.min_node = z.right
                self.consolidate()
            self.n -= 1
        return z
                    
    def consolidate( self ):
        A = [None for i in range(self.n)]
        rootList_nodes = [w for w in self.traverse ( self.root_list )]
        for w in rootList_nodes:
            x = w
            d = x.degree
            while A[d]:
                y = A[d]
                if x.key > y.key:
                    temp = x
                    x, y = y, temp
                self.heap_link(y, x)
                A[d] = None
                d += 1
            A[d] = x
        self.min_node = None
        for i in range(len(A)):
            if A[i] is not None:
                if self.min_node == None:
                    #self.insert_root_list( A[i ]) # here
                    self.min_node = A[i]
                else:
                    #self.insert_root_list( A[i ])
                    if A[i].key < self.min_node.key:
                        self.min_node = A[i]
        return
                
    def heap_link( self, y, x):
        self.remove_from_root_list( y ) # remove y from root list of H
        y.left = y.right = y
        # make y a child of x
        self.insert_parent_list(x,y)
        y.parent = x
        x.degree += 1
        y.mark = False
        
    def insert_parent_list( self, parent, y ):
        if parent.child == None:
            parent.child = y
        else:
            y.right = parent.child.right
            y.left = parent.child
            parent.child.right.left = y
            parent.child.right = y
            
    def decrease_key( self, x, k):
        if k > x.key:
            print(" Error: New key greater than current key")
            return
        x.key = k
        y = x.parent
        if y and x.key < y.key:
            self.cut(x,y)
            self.cascading_cut(y)
        if x.key < self.min_node.key:
            self.min_node = x
            
    def cut(self, x, y):
        self.remove_from_child_list(y,x)
        y.degree -= 1
        self.insert_root_list( x )
        x.parent = None
        x.mark = False
        
    def cascading_cut(self, y ):
        z = y.parent
        if z:
            if y.mark == False:
                y.mark = True
            else:
                self.cut(y,z)
                self.cascading_cut( z )
                
    def remove_from_child_list(self,parent,x):
        if parent.child == parent.child.right: # x is parent's only child
            parent.child = None
        elif parent.child == x:
            parent.child = x.right
            x.right.parent = parent
        x.left.right = x.right
        x.right.left = x.left
        
    def delete_node( self, x):
        self.decrease_key(x, -99999999)
        self.extract_min()
