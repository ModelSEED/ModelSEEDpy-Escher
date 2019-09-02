class EscherNode():
    
    def __init__(self, x, y, node_type, primary = True):
        self.x = x
        self.y = y
        self.primary = primary
        self.node_type = node_type
    
    @property
    def x(self):
        return self.x
    
    @property
    def y(self):
        return self.y