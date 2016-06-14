#!/usr/bin/env

# Simple module to double marker ids and check if they are the same.

def readUndoubledMarker(u_marker,location=-1):
    head = 2*int(u_marker)
    tail = head-1
    if location < 0:
        return map(str,[tail,head])
    else:
        return [(x,location) for x in map(str,[tail,head])]

def readDoubledMarker(d_marker):
    m,loc = d_marker.split('.')
    return (m,int(loc))

def isMate(ext1,ext2):
    m1 = int(ext1)
    m2 = int(ext2)
    if m1%2 == 1:
        m1 += 1
    if m2%2 == 1:
        m2 += 1
    if m1 == m2:
        return True
    return False


def makeMate(marker):
    if marker.__class__ == 'string'.__class__:
        if '.' in marker:
            m = int(marker[:marker.find('.')])
            id = marker[marker.find('.'):]                        
        else:        
            m = int(marker)
            id = ''
        if m%2 == 1:
            m += 1
        else:
            m -= 1
        return str(m)+id
    elif marker.__class__ == (1,2).__class__:        
        m = int(marker[0])
        id = marker[1]
        if m%2 == 1:
            m += 1
        else:
            m -= 1
        return (str(m),id)            
    return TypeError

