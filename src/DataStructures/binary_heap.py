#!/usr/bin/env

# Binary heap class
# Public functions:
#           constructor
#           insert
#           pop
#           peek
#           __len__
# Private functions:
#           parent, leftChild, rightChild           
#           heapifyUp
#           heapifyDown
#           swap
#           __setitem__
#           __getitem__

class BinaryHeap:
    # Constructor
    def __init__(self,item_list=[]):
        self.array = []
        for it in item_list:
            self.insert(it)
    # Insert item            
    def insert(self,it):
        self.array.append(it)
        self.heapifyUp()
    # Return parent index        
    def parent(self,index):
        return (index-1)>>1
    # Return left child index
    def leftChild(self,index):
        return (index<<1)+1
    # Return right child index
    def rightChild(self,index):
        return (index<<1)+2
    # Peek at highest priority element
    def peek(self):
        return self[0]
    # Set item value, given index
    def __setitem__(self,i,val):
        self.array[i] = val
    # Get item value, given index
    def __getitem__(self,i):
        return self.array[i]
    def __len__(self):
        return len(self.array)
    # Heapify up after adding element to heap
    def heapifyUp(self):
        i = len(self)-1
        while i != 0:
            p_index = self.parent(i)
            if self[i] < self[p_index]:
                self.swap(i,p_index)
                i = p_index
            else:
                break
    # Heapify down after removing most preferred element            
    def heapifyDown(self):
        highest_priority = 0 
        prev = 0
        l = self.leftChild(prev)
        r = self.rightChild(prev)
        while prev < len(self)-1 and \
                (l < len(self) and self[prev] >= self[l]) or \
                (r < len(self) and self[prev] >= self[r]):
            if l < len(self) and self[highest_priority] > self[l]:
                highest_priority = l
            if r < len(self) and self[highest_priority] > self[r]:
                highest_priority = r
            if highest_priority != prev:                
                self.swap(highest_priority,prev)
                prev = highest_priority
                l = self.leftChild(prev)
                r = self.rightChild(prev)
            else:
                break
    # Given two indices, swap their values
    def swap(self,index0,index1):
        self[index0],self[index1] = self[index1],self[index0]
    # Extract first element        
    def pop(self):
        if len(self) == 0:
            return None
        self.swap(0,len(self)-1)
        priority_item = self.array.pop()
        self.heapifyDown()
        return priority_item
