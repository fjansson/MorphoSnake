#!/usr/bin/env python3 
#
# MorphoSnake - Python Morphometry program for 2D shapes
#
# Fredrik Jansson, Pirom Konglerd, Catherine Reeb 2014-2017
# Released under the GNU General Public License v3.0.
# https://github.com/fjansson/MorphoSnake
#
# written for Python3

import mahotas as mh
import matplotlib as mpl
import matplotlib.colors
import matplotlib.pyplot as plt
import numpy as np
from skimage import morphology
import networkx as nx
import json
import sys
import os.path
import pickle


# list of all keys measuring lengths - these are scaled with px_mm
length_keys = ['dia', 'branchlength', 'branchlength_e', 'W_mean', 'W_max', 'W_min', 'W_std', 'apicaldist']

# list of all keys measuring angles - these are converted to degrees
angle_keys = ['alpha', 'alpha_e']


# colors and other visual settings
background_color = (0.2,0.2,0.2)
leaf_color = (0.1,0.5,0)
skel_color = (0.2,0.9,0)
node_color = (1.0,0.2,1)
node_color_fertile = (1.0,1.0,0.0)
node_color_carac = (1.0,0.0,0.0)
node_alpha = 0.8
node_edge_color = 'magenta'
edge_color = 'magenta'
edge_width = 4
edge_alpha = 0.6
terminal_disk_color = (1,0,0.4)
terminal_disk_alpha = 0.2
root_disk_color = 'y'
root_disk_alpha = .5

# color maps for imshow
leaf_colors = matplotlib.colors.ListedColormap([leaf_color, background_color]) # can also give triples of rgb floats: (1,.3,.1)
skel_colors = matplotlib.colors.ListedColormap([skel_color])

# from http://rosettacode.org/wiki/Zhang-Suen_thinning_algorithm
def neighbours(i, x, y):
    '''Return 8-neighbours of point p1 of picture, in order'''
    x1, y1, x_1, y_1 = x+1, y-1, x-1, y+1
    #print((x,y))
    return [i[y+1][x], i[y+1][x+1], i[y][x+1], i[y-1][x+1],  # P2,P3,P4,P5
            i[y-1][x], i[y-1][x-1], i[y][x-1], i[y+1][x-1]]  # P6,P7,P8,P9

# from http://rosettacode.org/wiki/Zhang-Suen_thinning_algorithm
# count the number of 0->1 transitions in a list
def transitions(neighbours):
    n = neighbours + neighbours[0:1]    # P2, ... P9, P2
    return sum((n1, n2) == (0, 1) for n1, n2 in zip(n, n[1:]))

# squared distance between two points p and q, given as tuples
def distsq(p, q):
    return (p[0]-q[0])**2 + (p[1]-q[1])**2

def dist(p, q):
    return np.sqrt(distsq(p, q))

# return the shortest distance from p to any point in s.
# p is an (x, y) tuple, s is a list of (x, y) tuples .
def distToSet(s, p):
    d = min([distsq(p, q) for q in s])
    return np.sqrt(d)

# p is an (x, y) tuple, s is a list of (x, y) tuples .
# return the point in s closest to p.
def findClosest(s, p):
    ds = [distsq(p, q) for q in s]
    return s[ds.index(min(ds))]
   


# find terminals and junctions, by counting 01 crossings in the neighborhood.
# 1 crossing       - terminal
# 2 crossings      - part of line, not reported
# 3 or 4 crossings - a junction
#
# Also, a pixel with 5 neighbors is probably a junction,
# but these will be missed by the crossing rules.
# However, 5-neighbor pixels can also appear between two closely spaced
# well-behaved junctions.
#
# All pixels with 5 neighbors are printed (as a warning)
# and added to the junction list, if further than 3 px from an existing junction

def terminals(skel):
    terminals = []
    junctions = []
    Y,X = skel.shape
    for y in np.arange(1,Y-1):
        for x in np.arange(1,X-1):
            if skel[y][x]:
                neigh = neighbours(skel,x,y)
                nn =  sum(neigh)
                n = transitions(neigh)
                if n == 1:
                    terminals += [(x,y)]
                if n >= 3:
                    junctions += [(x,y)]
                if nn >= 5:
                    print("Found a pixel with %d neighbors at (%d, %d), assuming a junction"%(nn,x,y))
                    if distToSet(junctions, (x,y)) < 3:
                        print("...but it is too close to an existing junction")
                    else:
                        junctions += [(x,y)]
                if n > 4:
                    print("Strange junction at (%d,%d) with %d transitions."%(x,y),n)
    return terminals, junctions
    
# offsets to the 8 neighbor pixels around one pixel
NeighborOffsets = [(1,0), (1,1), (0,1), (-1,1), (-1,0), (-1,-1), (0, -1), (1, -1)]


# Walk the skeleton starting at p.  Stop when meeting a junction or a
# terminal, or when it's impossible to continue because there are no
# unvisited neighbors.  Returns the final point, a string indicating
# the reason for stopping, and a length estimate.

# The walk assumes that the skeleton is as thin as possible, i.e. no
# staircases.  Every pixel in the non-junction parts, is assumed to
# have exactly two neighbors.  This seems to work with the skeleton
# returned by skicit image's thinning-based skeletonizer.

def skeletonWalk(skel, visited, p, startJunction, junctions, terminals, dmap):
    x,y = p
    l = 0 # distance estimate
    thickness = []
    skelPixels = []
    # print ("skeletonWalk at (%d,%d)"%p)
    while True:
        visited[y][x] = True
        thickness.append(dmap[y][x])
        skelPixels.append((x,y))
        c = 0
        nxt = (0,0)
        dl = 0
        for dx,dy in NeighborOffsets:
            if skel[y+dy][x+dx]  and distsq( (x+dx,y+dy), startJunction) > 2:
                if (x+dx,y+dy) in junctions:
                    # we found a new junction
                    #print ("junction", (x+dx,y+dy))
                    thickness.append(dmap[y+dy][x+dx])
                    skelPixels.append((x+dx,y+dy))
                    if not visited[y+dy][x+dx]:
                        return ((x+dx,y+dy), 'j', l + np.hypot(dx,dy), thickness, skelPixels) 
                                              # length estimate l + final step
                    else:
                        print("   loop detected! (%d, %d)"%(x+dx,y+dy))
                        return ((x+dx,y+dy), 'l', l + np.hypot(dx,dy), thickness, skelPixels) 
                                             # length estimate l + final step
                if (x+dx,y+dy) in terminals:
                    # we found a terminal
                    #print ("terminal", (x+dx,y+dy))
                    thickness.append(dmap[y+dy][x+dx])
                    skelPixels.append((x+dx,y+dy))
                    return ((x+dx,y+dy), 't', l + np.hypot(dx,dy), thickness, skelPixels)
                if not visited[y+dy][x+dx]:
                    c += 1
                    nxt = (x+dx, y+dy)
                    dl = np.hypot(dx,dy)

        # estimate the length of this branch
        l += dl

        if c > 1:
            # There is more than one unvisited neighbor, without this pixel or a neighbor being
            # a junction - i.e. not in the junctions list.
            # This can happen where two diagonal lines cross so that they form a 2x2 square.
            print("Found %d neighbors without finding a junction at (%d, %d)"%(c,x,y))
            return ((x,y), 's', l, thickness, skelPixels)

        if c == 0:
            # There is nowhere to go, but this isn't a terminal.
            # This happens when a junction has an unnecessary neighbor pixel, 
            # which can apppear as a very short branch.
            if l > 0: #no warning if the branch is short
                print ("Found no neighbors without finding a terminal at (%d, %d) length %f"%(x,y,l))
            return (None, 'e', l, thickness, skelPixels)

        #there is one unvisited neighbor, go to it
        x,y = nxt
        

# Construct a tree from the skeleton graph, starting at p, which is
# assumed to be a junction or a terminal.
#
# G is a NetworkX graph, to which junctions and terminals are added as they are discovered.
# Nodes in G are identified by their coordinate tuple.
def buildTree(skel, visited, dmap, p, junctions, terminals, G):
    #print ("buildTree at (%d,%d)"%p)

    # find unvisited neighbors of this pixel
    neighbors = []
    x,y = p
    visited[y][x] = True
    for dx,dy in NeighborOffsets:
        if skel[y+dy][x+dx] and not visited[y+dy][x+dx]:
            if (x+dx, y+dy) in junctions:
                #our starting junction has a junction as a neighbor, this can happen.
                #we add *it's* neighbors to the neighbor list
                #some pixels will be there twice, which is OK
                #if 3 junctions can be neighbors, this must be done in more passes.
                neighbors += [(x+dx+ox,y+dy+oy) for (ox,oy) in NeighborOffsets]
            else:
                #visited[y+dy][x+dx] = True
                neighbors += [(x+dx,y+dy)]
            
    # start a skeletonWalk at every neighbor
    for n in neighbors:
        pp,t,l,thickness,skelPixels = skeletonWalk(skel, visited, n, p, junctions, terminals, dmap)
        l += np.hypot(n[0]-p[0],n[1]-p[1])  #add length of first step

        if t in ('j', 't', 'l', 's'):
            thickness = np.sqrt(thickness) # get real widths, from squared ones
            G.add_edge(p, pp)
            G[p][pp]['branchlength'] = l
            G[p][pp]['branchlength_e'] = np.sqrt(distsq(p, pp))  #Euclidian distance between end points

            G[p][pp]['skelPixels'] = skelPixels         #all the pixels visited during the skeleton walk
            # statistics of the branch width
            G[p][pp]['thicknessProfile'] = 2*thickness  #the full thickness profile along the skeleton
            G[p][pp]['W_min']  = 2*np.min(thickness) 
            G[p][pp]['W_max']  = 2*np.max(thickness)
            G[p][pp]['W_mean'] = 2*np.mean(thickness)
            G[p][pp]['W_std']  = 2*np.std(thickness)
        if t == 'j' or t == 's':
            buildTree(skel, visited, dmap, pp, junctions, terminals, G)
        #if t == 't':            
        #if t == 'l':  # this branch forms a loop. 
            # added to the networkX graph above in any case.


# plot all edges in the graph G
# but better to use networkx's builtin plots
def plotG(G):
    for e in G.edges():
        plt.plot((e[0][0], e[1][0]), (e[0][1], e[1][1]), 'c-')


# find the Strahler order of each node in the graph
# S = max(S_children) + delta
#  where delta = 1 if two (or more) children have the max value
#  for a leaf, S = 1
#
# we want Strahler order for edges as well. 
# for each edge, the number is the smaller of the Strahler numbers of 
# the two nodes.
# S(e) = min (S(e[0]), S(e[1]))
#
# also compute 'level':
# the level of a leaf is 1. The level of a node is max(children's levels) + 1
# property: a node has a higher level than any of it's children
# used to identify parent when measuring angles.

def StrahlerOrder(G, root):
    # initialize all nodes to 1
    for p in G:
        G.node[p]['Strahler'] = 1
        G.node[p]['level'] = 0      # longest path from this node to a leaf - leaves have level 1.
    
    nodes = nx.dfs_postorder_nodes(G, root)
    for p in nodes:
        neigh = nx.all_neighbors(G, p)
        delta = 0
        mx = 0
        level = 0
        for q in neigh:
            if G.node[q]['Strahler'] == mx: # we met the max value again
                delta = 1 
            if G.node[q]['Strahler'] > mx:
                mx = G.node[q]['Strahler']
                delta = 0
            level = max(level, G.node[q]['level'])
        G.node[p]['Strahler'] = mx + delta
        G.node[p]['level'] = level + 1 
        
    # assign Strahler orders to edges too
    for e in G.edges():
        s = min(G.node[e[0]]['Strahler'], G.node[e[1]]['Strahler'])
        G[e[0]] [e[1]] ['Strahler'] = s
        s = min(G.node[e[0]]['level'], G.node[e[1]]['level'])
        G[e[0]] [e[1]] ['level'] = s

    # determine the parent of each node, based on the node 'level'
    for n in G:
        m = 0
        p = None
        for n2 in G[n]:
            if G.node[n2]['level'] > m and G.node[n2]['level'] > G.node[n]['level']:
                m = G.node[n2]['level']
                p = n2
        G.node[n]['parent'] = p
    return

# return the diameter of the minimal disk centered at p,
# which touches the background
# look at all pixels within a square with side 2 R, if no 
# background is found, double R and try again
#
# This function is not in use -  mahotas.distance used instead, faster and gives identical results.
def maxDisk(im, p):
    
    R = 10
    rminsq = 100000000

    X,Y = p
    Ysize, Xsize =  im.shape
    while True:
        for y in np.arange(max(0, Y-R), min(Ysize, Y+R+1)):
            for x in np.arange(max(0, X-R), min(Xsize, X+R+1)):
                if not im[y][x]:
                    rminsq = min(rminsq, (x-X)**2+(y-Y)**2)

        if rminsq < R**2: 
            # rminsq can be trusted only if all pixels at this distance are in the square
            return 2*np.sqrt(rminsq)
        R *= 2


# find the point on path which crosses the circle centered on P with radius R
# used for branch angles.
# For short branches, it's possible that the branch is completely inside the disk.
# then, return the more distant end point of the branch.
def crossing(P, R, path):
    for p in path:
        r = np.sqrt(distsq(p, P))
        if np.abs(r - R) < 1.0:
            return p

    r0 = np.sqrt(distsq(path[0], P))
    if r0 > r:
        return path[0]
    else:
        return p;


# determine the angle between the vectors ab and ac, given points a, b, and c as tuples    
def angle(a, b, c):
    dot = (b[0]-a[0]) * (c[0]-a[0]) + (b[1]-a[1]) * (c[1]-a[1])
    Dab = dist(a, b)
    Dac = dist(a, c)
    return  np.arccos(dot/(Dab*Dac))
    
# measure node diameters
def measureDia(G, dmap):
    for n in G:
        # d1 = maxDisk(im, n) # find disks manually 
        x,y = n
        d = 2*np.sqrt(dmap[y][x]) 
        G.node[n]['dia'] = d

def measureApicalDist(G):
    for n in G:
        # measure apical distance
        # we need it only for terminal nodes, but we measure for all nodes
        # so that all nodes have the same keys
        r = 1000000
        for n2 in G:
            if G.node[n2]['level'] == 1 and n2 != n:
                r = min(r, dist(n, n2))
        G.node[n]['apicaldist'] = r
                    
        
# Measure branch angles in the graph. Must be updated if the root moves.
# measureDia and StrahlerOrder should be called before this
def measureAngles(G):
    # for each edge, create a dictionary mapping end point to theta angle
    # e[2] is the data dictionary associated with each edge
    for e in G.edges(data=True):
        e[2]['theta'] = {}
    
    for n in G:
        # iterate over neighbors, determine angle to every one
        d = G.node[n]['dia']
        for n2 in G[n]:
            # p is the point where the skeleton crosses this node's disk
            p = crossing(n, d/2.0, G[n][n2]['skelPixels'] )
            dx = p[0] - n[0]  # vector from n to p
            dy = p[1] - n[1]
            dist = np.hypot(dx, dy)
            theta = np.arctan2(-dy, dx) #angle for this vector, in radians, measured counterclockwise from the +x axis
                                        # -dy since y grows downwards here, and upwards traditionally
            #print (dist, theta, d/2.0)
            
            G[n][n2]['theta'][n] = theta

        # for each child branch, measure the angle relative to the parent branch, 'alpha'
        # also measure the angle between the straigh lines to the parent node and the child node, 'alpha_e'
        theta_p = None
        parent = G.node[n]['parent']
        if parent != None:
            theta_p = G[n][parent]['theta'][n]
        for n2 in G[n]: 
            if n2 != parent: # a child branch
                if theta_p != None:
                    theta = G[n][n2]['theta'][n]
                    alpha = abs(theta - theta_p)
                    if alpha > np.pi:
                        alpha = 2*np.pi - alpha
                    alpha_e = angle(n, n2, parent)
                else:  # no parent - this is the root
                    alpha   = np.pi  # these are excluded from the averages in report() based on adjacency to the root
                    alpha_e = np.pi
                G[n][n2]['alpha'] = alpha
                G[n][n2]['alpha_e'] = alpha_e


# remove spurious brances
# a branch is spurious if:
#    it's a terminal branch AND
#       it's shorter than minLength
#       OR
#       it has a disk diameter < minDia AND branch length < maxRemoveLen
# The root cannot be deleted.
#
# return a list of the removed nodes
def cleanup(G, root, minLength=5, minDia=8, maxRemoveLen = 20):
    removed = []
    for n in G:
        if G.degree(n) == 1 and n != root:
            p = list(G.neighbors(n))[0] #the parent (only neighbor)
            if G[n][p]['branchlength'] < minLength or (G.node[n]['dia'] < minDia and G[n][p]['branchlength'] < maxRemoveLen):
                removed.append(n)
    G.remove_nodes_from(removed)
    print('removed nodes', removed)
    return removed

# find the node in G closest to p and delete it
def deleteNode(G, p):
    print('deleteNode (%5.1f, %5.1f)'%p)
    p = findClosest(list(G.nodes()), p)
    print('closest node is (%5.1f, %5.1f)'%p)

    parent = G.node[p]['parent']
    for n in G.neighbors(p):
        if n!=parent:
            G.add_edge(n,parent)
            G[parent][n]['branchlength'] = G[parent][p]['branchlength']+G[p][n]['branchlength']
            G[parent][n]['branchlength_e'] = np.sqrt(distsq(parent,n))
            G[parent][n]['skelPixels'] = G[parent][p]['skelPixels'] + G[p][n]['skelPixels']
            G[parent][n]['W_min'] = min(G[parent][p]['W_min'], G[p][n]['W_min'])
            G[parent][n]['W_max'] = max(G[parent][p]['W_max'], G[p][n]['W_max'])
            G[parent][n]['W_mean'] = (G[parent][p]['W_mean']+G[p][n]['W_mean'])/2
            G[parent][n]['W_std'] = (G[parent][p]['W_std']+G[p][n]['W_std'])/2
            print("added edge")
            
    G.remove_node(p)


    
# plot the graph G with nodes and edges.
# keep handles to everything, so that they can be updated - i.e. removed and added again
def plot_graph(G):
    global nodes, edges, edge_labels, node_labels, rad
    # remove old elements from the plot if they exist
    if nodes:
        nodes.remove()
    if edges:
        edges.remove()
    if node_labels:
        for p in node_labels:
            node_labels[p].remove()

    if edge_labels:
        for e in edge_labels:
            edge_labels[e].remove()

            
    pos = {}
    nlabels = {}
    elabels = {}
    rad=[]

    # build a dictionary mapping nodes to their positions needed for
    # drawing the graph the nodes themselves are also their positions, but
    # there seems to be no way to use that directly.  Also make dictionary
    # of node labels

    for p in G.nodes:
        pos[p] = p
        #nlabels[p] = "%.1f"%G.node[p]['dia']
        nlabels[p] = "%d"%G.node[p]['Strahler']
        rad.append(G.node[p]['dia'] * .5)
        
    #for e in G.edges():
        #elabels[e] = G.edge[e[0]][e[1]]['Strahler']
        #elabels[e] = G.edge[e[0]][e[1]]['level']
        #if G.edge[e[0]][e[1]]['alpha'] == None:
        #    elabels[e] = '-'
        #else:
        #    elabels[e] = "%.0f"%(G.edge[e[0]][e[1]]['alpha']*180/np.pi) 
    
    nodes = nx.draw_networkx_nodes(G, pos, node_color = [G.node[p]['node_color'] for p in G])
    plt.setp(nodes, edgecolors=node_edge_color, picker=True)
    
    edges = nx.draw_networkx_edges(G, pos, width = edge_width, alpha=edge_alpha, edge_color=edge_color)
    node_labels = nx.draw_networkx_labels(G, pos, nlabels) 
    edge_labels = nx.draw_networkx_edge_labels(G, pos, elabels) 
    
# save CSV files with data for all edges and nodes.
# columns are saved in the order they are listed in the edge_keys and node_keys below
def saveTreeText(G, edgeName, nodeName):
    edge_keys = ['Strahler', 'level', 'branchlength', 'branchlength_e', 'alpha', 'alpha_e', 'W_mean', 'W_max', 'W_min', 'W_std']
    node_keys = ['Strahler', 'level', 'dia']

    px_mm_factor = 1
    if px_mm != None:
        px_mm_factor = px_mm

    # make a dictionary of conversion factors
    conversion = {}
    for k in length_keys:
        conversion[k] = 1.0/px_mm_factor
    for k in angle_keys:
        conversion[k] = 180/np.pi
             
    eheader = ''
    for k in edge_keys:
        eheader += k + ', ' 
    nheader = ''
    for k in node_keys:
        nheader += k + ', ' 
    
    output = [[e[2][k] * (conversion[k] if k in conversion else 1)
               for k in edge_keys] for e in G.edges(data=True)]
    output.sort() # sort on Strahler order
    np.savetxt(edgeName, output, fmt = '%3d,%3d' + (len(edge_keys)-2) * ', %8.3f', header=eheader, comments='#')

    output = [[n[1][k] * (conversion[k] if k in conversion else 1) 
               for k in node_keys] for n in G.nodes(data=True)]
    output.sort() # sort on Strahler order
    np.savetxt(nodeName, output, fmt = '%3d,%3d' + (len(node_keys)-2) * ', %8.3f', header=nheader, comments='#')


def report(G):
    # extract values from the data dictionaries of all the nodes and all the edges 

    # list the keys we are interested in
    node_keys = ['dia', 'apicaldist']
    edge_keys = ['branchlength', 'branchlength_e', 'alpha', 'alpha_e', 'W_mean', 'W_max']
    

    all_keys = node_keys + edge_keys
    
    data    = {} # all branches and nodes
    data_t  = {} # terminal branches and nodes
    data_nt = {} # non-terminal branches and nodes
    
    # extract the data into new dictionaries
    # EXCLUDE the root's branches from alpha and alpha_e ?
    for k in edge_keys:
        if k == 'alpha' or k == 'alpha_e':
            # Special treatment of the angles.
            # Want to exclude the root's child branches which do not have a well-defined angle.
            data[k]    = []
            data_t[k]  = []
            data_nt[k] = []
            
            for e in G.edges(data=True):
                v = e[2][k]
                if e[0] != root and e[1] != root:
                    # this edge does not go to the root
                    # add the data to the appropriate lists
                    data[k].append(v)
                    if e[2]['level'] == 1:
                        data_t[k].append(v)
                    else:
                        data_nt[k].append(v)
                # else:
                #     print ('Excluded ' + str(v))
            data[k] = np.array(data[k])
            data_t[k] = np.array(data_t[k])
            data_nt[k] = np.array(data_nt[k])
        else:
            #standard treatment, keep all edges
            data[k]    = np.array([e[2][k] for e in G.edges(data=True) ])
            data_t[k]  = np.array([e[2][k] for e in G.edges(data=True) if e[2]['level'] == 1])
            data_nt[k] = np.array([e[2][k] for e in G.edges(data=True) if e[2]['level'] != 1])
        
    for k in node_keys:
        data[k]    = np.array([n[1][k] for n in G.nodes(data=True) ])
        data_t[k]  = np.array([n[1][k] for n in G.nodes(data=True) if n[1]['level'] == 1])
        data_nt[k] = np.array([n[1][k] for n in G.nodes(data=True) if n[1]['level'] != 1])
                
        
    # convert to degrees from radians
    for d in [data, data_t, data_nt]:
        for k in d: 
            if k in angle_keys:
                d[k] *= 180/np.pi
            
    # convert to mm
    if px_mm != None:
        unit = 'mm'
        for d in [data, data_t, data_nt]:
            for k in d: 
                if k in length_keys:
                    d[k] /= px_mm
        px_mm_factor = px_mm
    else:
        unit = 'PX'        
        px_mm_factor=1.0 #local variable used when scaling for the report file        
    
    print()
    print('RESULTS in %s and degrees'%unit)
    print('-------------------------')
    print('%d nodes, %d branches '%(len(data['dia']), len(data['alpha'])) )
    print('%d terminal nodes, %d terminal branches '%(len(data_t ['dia']), len(data_t ['alpha'])) )
    print('%d internal nodes, %d internal branches '%(len(data_nt['dia']), len(data_nt['alpha'])) )
    print()

    avg    = {}
    std    = {}
    avg_t  = {}
    std_t  = {}
    avg_nt = {}
    std_nt = {}
    for k in all_keys:
        avg[k]    = np.average(data[k])
        std[k]    = np.std(data[k])
        avg_t[k]  = np.average(data_t[k])
        std_t[k]  = np.std(data_t[k])
        avg_nt[k] = np.average(data_nt[k])
        std_nt[k] = np.std(data_nt[k])
        
    print('                      all      |     terminals   |  non-terminals')
    print('      parameter   avg    std   |    avg    std   |    avg    std ')    

    for k in all_keys:
        print(k.rjust (15), "%6.2f"%avg[k], "%6.2f"%std[k], ' | ', "%6.2f"%avg_t[k], "%6.2f"%std_t[k],  ' | ', "%6.2f"%avg_nt[k], "%6.2f"%std_nt[k])

    header_text = "#name "

    for n in ['all', 'term', 'internal']:
        for k in all_keys:
            header_text += ", %s_%s_avg, %s_%s_dev"%(n,k,n,k)

    report_text = basename.rjust(15)
    for (a,d) in [(avg, std), (avg_t, std_t), (avg_nt, std_nt)]:
        for k in all_keys:
            report_text += ", %6.2f, %6.2f"%(a[k], d[k])

    report_text += ', %5.0f'%px_mm_factor
    header_text += ', px/mm'

    #print (header_text)
    #print (report_text)
    return report_text, header_text

# return the directory of the leafsnake.py file
# used to find the leaf database
# http://stackoverflow.com/questions/4934806/python-how-to-find-scripts-directory
def getScriptPath():
    return os.path.dirname(os.path.realpath(sys.argv[0]))
    

# find the distance between a point p and the line segment between l1 and l2.
# http://paulbourke.net/geometry/pointlineplane/
def distPointLine(l1, l2, p):
    x = p[0]
    y = p[1]
    l1x = l1[0]
    l1y = l1[1]
    l2x = l2[0]
    l2y = l2[1]
    u = (x - l1x)*(l2x - l1x) + (y - l1y)*(l2y - l1y)
    u /= distsq(l1, l2)

    # px, py is the closest poitn on the line.
    # px = l1x + u (l2x - l1x)
    # py = l1y + u (l2y - l1y)
    if 0 < u and u < 1: #the closest point is inside the segment
        px = l1x + u * (l2x - l1x)
        py = l1y + u * (l2y - l1y)
        return dist((px, py), p)
    elif u < 0: #the closest point is l1
        return dist(l1, p)
    else:
        return dist(l2, p)

def findClosestEdge(G, p):
    emin = None
    dmin = 1e100
    for e in G.edges():
        d = distPointLine(e[0], e[1], p) 
        if d < dmin:
            dmin = d
            emin = e
    return emin

def findClosestSkel(skel,p):
   # we store coordinates of pixels in the skeleton in a list
   listskel = []
   Y,X = skel.shape
   for y in np.arange(1,Y-1):
       for x in np.arange(1,X-1):
           if skel[y][x]:
               listskel += [(x,y)]
   return findClosest(listskel,p)


def addNodeSkel(G,p,skel,junctions,terminals,dmap):
    visited = np.zeros_like(skel)
    neighbors = []
    x,y = p
    visited[y][x] = True
    for dx,dy in NeighborOffsets:
        if skel[y+dy][x+dx] and not visited[y+dy][x+dx]:
            neighbors += [(x+dx,y+dy)]
            
    for n in neighbors:
        pp,char,l,thickness,skelPixels = skeletonWalk(skel, visited, n, p, junctions, terminals, dmap)
        
        l += np.hypot(n[0]-p[0],n[1]-p[1])  #add length of first step
        
        if char in ('j', 't', 'l', 's') and pp in G:
            thickness = np.sqrt(thickness) 
            G.add_edge(p, pp)
            G[p][pp]['branchlength'] = l
            G[p][pp]['branchlength_e'] = np.sqrt(distsq(p, pp)) 
            G[p][pp]['skelPixels'] = skelPixels      
            G[p][pp]['thicknessProfile'] = 2*thickness 
            G[p][pp]['W_min']  = 2*np.min(thickness) 
            G[p][pp]['W_max']  = 2*np.max(thickness)
            G[p][pp]['W_mean'] = 2*np.mean(thickness)
            G[p][pp]['W_std']  = 2*np.std(thickness)

    #add the diameter for the added node
    x,y = p
    d = 2*np.sqrt(dmap[y][x]) 
    G.node[p]['dia'] = d
    G.node[p]['node_color'] = node_color


def initializeColorNode(G):
    for p in G:
        G.node[p]['node_color'] = node_color
        G.node[p]['fertile'] = False
        G.node[p]['tag'] = []

###################################################################################
# Main code starts here
#

# take the file name from the command line, if given
if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = 'img/f.png'


# read database of leaves
leaves = {}

fileHandle = None
try:
    leavesFile = getScriptPath() + '/leaves.json'
    print('Database file: ' + leavesFile)
    fileHandle = open(leavesFile, 'r')
except:
    print('Could not load the leaf data base')
    
if fileHandle != None:
    # if the file exists, but the JSON is not correct, this will fail
    # This is on purpose, since if we continue and save the database at the end,
    # the data in it will be overwritten
    leaves = json.load(fileHandle)

    
path_basename = os.path.splitext(filename)[0]
print('path_basename', path_basename)
basename = os.path.splitext(os.path.basename(filename))[0]

leafData = {}
try:
    #find the current leaf, based on file name
    leafData = leaves[basename]
except:
    print('%s was not found in the leaf data base'%basename)
    # add an entry for this leaf to the ditionary
    leaves[basename] = leafData
    
print('Data for this leaf: ' + str(leafData))

px_mm = None
root = None

try:
    px_mm = leafData['px_mm']
except:
    print('No resolution found for leaf %s'%basename)

try:
    root = tuple(leafData['root']) 
except:
    print('No root found for leaf %s'%basename)
    
print('Resolution: ' + str(px_mm) + ' px/mm')
print('Root:       ' + str(root))
print()

print('Reading image %s'%filename)
img = mh.imread(filename)

# color image handling - works if the background is brighter than the object
if len(img.shape) == 3: #color image
    #convert to grayscale
    img = mh.colors.rgb2gray(img, dtype=np.uint8)

# thresholding
T_otsu = mh.otsu(img) # finds a numeric threshold
img = (img > T_otsu)  # make image binary

# invert the image (just because the test image is stored the other way)
img = ~img 

# close single-pixel holes. Also makes the skeletonization much more well-behaved,
# with less tiny branches close to terminals.
#
# This can create loops if two branches are separated by < 3 px of background
img = mh.close(img)

print('Thinning...')
# skeletonization from scikit image.
# Zhang-Suen algorithm (apparently with staircase removal)
skel = morphology.skeletonize(img)

# Try mahotas skeletonization. Makes many spikes.
# works better after mh.close(), but still splits many tips in two.
# Also gives staircases in the skeleton.
#skel = mh.thin(img)  
#skel = morphology.skeletonize(skel)  # one pass of the other skeletonization to remove staircases


################

# find terminals and junctions. t and j are in the format [(x1, y1), (x2, y2), ...]
print('Features...')
t, j = terminals(skel) 

if root == None:
    # unpack the tuples to separate lists of x:s and y:s, for plotting and root selection
    #tx = [x[0] for x in t]
    ty = [x[1] for x in t]

    # find the index of the lowest terminal. Use that as the root, for now.
    iroot = ty.index(max(ty))

    # find the lowest node, to use as root
    root = t[iroot]
else:
    if root not in t+j:
        newroot = findClosest(t+j, root)
        print('Moving the root to a node in the tree. New root is', str(newroot), 'old root was', str(root), 'distance', dist(newroot, root))     
        root = newroot


print('Plotting')
# make the skeleton image's background transparent
# create a copy of the boolean skeleton for further use
skel2 = skel
skel    = np.ma.masked_where(skel==0, skel)
#visited = np.ma.masked_where (visited==0, visited)

fig = plt.figure()
fig2 = fig

# Plot images. Inversion here is just for nicer coloring
# interpolation=nearest is to turn smoothing off.

width = skel.shape[1]
height= skel.shape[0]

plt.axis((0,width,height,0))

plt.imshow(~img, cmap=leaf_colors, interpolation="nearest")
#plt.imshow(np.sqrt(dmap), cmap=mpl.cm.jet_r, interpolation="nearest")

plt.imshow(skel, cmap=skel_colors,  interpolation="nearest")
#plt.imshow(visited, cmap=mpl.cm.cool,  interpolation="nearest")

# update measures that depend on graph structure or root placement
def updateMeasures(G, root):
    print('  Strahler order...')
    StrahlerOrder(G, root)
    print('  apical distances...')
    measureApicalDist(G)
    print('  branch angles...')
    measureAngles(G)
    
def buildGraph(img, skel):
    print('Building tree...')
    G = nx.Graph()
    visited = np.zeros_like(skel) # a new array of zeros, same dimensions as skel
    
    print('  distance transform...')
    dmap = mh.distance(img)

    print('  constructing tree...')
    buildTree(skel, visited, dmap, root, j, t, G)

    # measure node diameters
    measureDia(G, dmap)

    # automatically remove bad nodes and branches, e.g. too small ones
    removed = cleanup(G, root)
    # show (automatically) removed nodes
    for x,y in removed:
        plt.gca().add_patch(plt.Circle((x,y), radius=4, alpha=.4))

    updateMeasures(G, root)
    initializeColorNode(G)
    print('Done.')
    return G

# read in the graph from a previous run, if it exists
try:
    G=nx.read_gpickle(path_basename+'_graph.pkl')
    print('Loaded graph from ' + path_basename+'_graph.pkl')
except:
    # could not read the graph. Constructing it now
    G = buildGraph(img, skel)
    
# handles to plot elements
nodes = None
edges = None
node_labels = None
edge_labels = None

rad=[]
plot_graph(G)

# semi-transparent circles on the nodes
Cercle = {}
for p,r in zip(G, rad):
	Cercle[p] = plt.Circle(p, radius=r, alpha=terminal_disk_alpha, color=terminal_disk_color)
	plt.gca().add_patch(Cercle[p])

root_patch = plt.Circle(root, radius=40, alpha=root_disk_alpha, color=root_disk_color)
plt.gca().add_patch(root_patch)


def updateCircles(G,C):
   ''' for p in G:
	if p not in G.nodes():
		Cercle[p].remove()
	else:
		Cercle[p].remove()
		Cercle[p] = plt.Circle(p, radius=G.node[p]['dia'], alpha=terminal_disk_alpha, color=terminal_disk_color)
		plt.gca().add_patch(Cercle[p])
'''
   for p in C:
       if p not in G:
           C[p].set_visible(False)
       else:
           C[p].set_visible(False)
           C[p] = plt.Circle(p, radius=G.node[p]['dia']/2, alpha=terminal_disk_alpha, color=terminal_disk_color)
           plt.gca().add_patch(C[p])

updateCircles(G,Cercle)
            
def setRoot(root):
    root_patch.center = root;

    # save the new root in database
    # convert to int from numpy type, for JSON to work later
    leafData['root'] = (int(root[0]), int(root[1]))
    
    
# a function called when the user clicks a node    
def onpick(event):
    global root
    # make the clicked node the new root
   
    # for some reason it's difficult to get the coordinates of the clicked node
    # so we use mouse coordinates and search for the closest node.
    root = findClosest(list(G.nodes()), (event.mouseevent.xdata, event.mouseevent.ydata))
    print('New root: ' + str(root))
    setRoot(root)
    updateMeasures(G, root)
    report(G)
    
    plot_graph(G)        
    fig.canvas.draw()
   

undo_stack = []

# a function called on keypress events
def keypress(event):
    global nodes, edges, node_labels, G, root
    print('press', event.key)
    sys.stdout.flush()
    if event.key=='d': # delete closest node
        p = findClosest(list(G.nodes()), (event.xdata, event.ydata))
        print('closest node is (%5.1f, %5.1f)'%p)
        if p == root:
            print('Cannot remove the root.')
            return
        undo_stack.append((G.copy(), root, Cercle))
        deleteNode(G,p)
        Cercle[p].remove()
        #report(G)
        
        updateMeasures(G, root)
        plot_graph(G)        
        fig.canvas.draw()
        
    if event.key=='x': # delete closest branch
        e = findClosestEdge(G, (event.xdata, event.ydata))
        undo_stack.append((G.copy(), root,Cercle))
        G.remove_edge(*e)

        updateMeasures(G, root)
        plot_graph(G)        
        fig.canvas.draw()

    if event.key=='u': # undo
        if len(undo_stack) > 0:
            print('Undo')
            G,root,C = undo_stack.pop()
            setRoot(root)
            updateMeasures(G, root)
            updateCircles(G,C)
            plot_graph(G)        
            fig.canvas.draw()
        else:
            print('No further undo')

    if event.key=='r': #re-build the graph from skeleton
        print('Rebuilding tree')
        undo_stack.append((G.copy(), root,Cercle))
        G = buildGraph(img, skel)
        plot_graph(G)
        updateCircles(G,Cercle)
        fig.canvas.draw()

    if event.key=='a': #marking of fertile nodes
        p = findClosest(list(G.nodes()), (event.xdata, event.ydata))
        print('this node is now fertile')
        undo_stack.append((G.copy(),root,Cercle))
        # G.node[p]['node_color'] = node_color_fertile
        G.node[p]['fertile'] = True
        plot_graph(G)
        fig.canvas.draw()

    if event.key == 'ctrl+a': #display all fertile nodes
        p = findClosest(list(G.nodes()), (event.xdata, event.ydata))
        print('display of all fertile nodes')
        undo_stack.append((G.copy(),root,Cercle))
        for p in G:
            if G.node[p]['fertile'] == True:
                G.node[p]['node_color'] = node_color_fertile
                plot_graph(G)
                fig.canvas.draw()

    if event.key == 'alt+e': #hide any latest selection
        p = findClosest(list(G.nodes()), (event.xdata, event.ydata))
        print('hide caracteristics')
        undo_stack.append((G.copy(),root,Cercle))
        for p in G:
            if G.node[p]['node_color'] != node_color:
                G.node[p]['node_color'] = node_color
                plot_graph(G)
                fig.canvas.draw() 

    if event.key == 'c': #choose a node and display all of his caracteristics
        p = findClosest(list(G.nodes()), (event.xdata, event.ydata))
        print('Write the nametag wanted for node (%5.1f,%5.1f):'%p)
        ntag = input()
        undo_stack.append((G.copy(),root,Cercle))
        (G.node[p]['tag']).append(ntag) 
        plot_graph(G)
        fig.canvas.draw()

    if event.key == 'b': # diameter modification
        p = findClosest(list(G.nodes()), (event.xdata, event.ydata))
        print('closest node is (%5.1f, %5.1f)'%p)
        undo_stack.append((G.copy(),root,Cercle))
        G.node[p]['dia'] = 2*dist(p,(event.xdata,event.ydata))
        updateMeasures(G, root)
        Cercle[p].remove()
        Cercle[p] = plt.Circle(p,radius=G.node[p]['dia']/2,alpha=terminal_disk_alpha,color=terminal_disk_color)
        plt.gca().add_patch(Cercle[p])
        plot_graph(G)        
        fig.canvas.draw()   
        
    if event.key == 'n': #add node
        p = findClosestSkel(skel2,(event.xdata,event.ydata))
        print('closest node is (%5.1f, %5.1f)'%p)
        undo_stack.append((G.copy(),root,Cercle))
        print('adding node');
        try:
            addNodeSkel(G,p,skel2,j,t,dmap)
        except:
            print('  distance transform...')
            dmap = mh.distance(img)
            addNodeSkel(G,p,skel2,j,t,dmap)

        Cercle[p] = plt.Circle(p, radius=(G.node[p]['dia']/2), alpha=terminal_disk_alpha, color=terminal_disk_color)
        updateCircles(G,Cercle)
        updateMeasures(G, root)
        report(G)
        plot_graph(G)        
        fig.canvas.draw()



        
# register the event callback functions
fig.canvas.mpl_connect('pick_event', onpick)
fig.canvas.mpl_connect('key_press_event', keypress)

# plot disks for the apical distance, just for testing
# for n in G:
#    if G.node[n]['level'] == 1:
#        r = G.node[n]['apicaldist']
#        plt.gca().add_patch(plt.Circle(n, radius=r, alpha=.4))

        
report(G)
plt.tight_layout()

# show the plot - program pauses here for as long as the window is open
plt.show()

saveTreeText(G, path_basename+'_branches.txt', path_basename+'_nodes.txt');
nx.write_gpickle(G, path_basename+'_graph.pkl')

def saveReport(reportFile, report_text, header_text):
    # if the file does not exist already, write the header
    present = os.path.isfile(reportFile)

    print('Saving leaf report in file ' + reportFile)
    try:
        f = open(reportFile, 'at')
        if not present:
            f.write(header_text+'\n')
            f.write(report_text+'\n')
    except:
        print('Error when saving the report')

report_text,header_text = report(G)
reportFile = 'results.txt'
saveReport(reportFile, report_text, header_text)

print('Saving the database...')
# save the leaf data base
of = open(leavesFile, "wt")
json.dump(leaves, of, sort_keys=True, indent=2, separators=(',', ': '))
print('done.')







