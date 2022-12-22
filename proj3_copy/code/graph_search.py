from heapq import heappush, heappop,heapify  # Recommended.
import heapq
#import heapdict
import numpy as np

from flightsim.world import World
#from heapdict import heapdict

from .occupancy_map import OccupancyMap # Recommended.

neighbourhood8 = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
neighbourhood26 = np.array([[-1,-1,-1],[-1,-1,0],[-1,-1,1],[-1,0,-1],[-1,0,0],[-1,0,1],[-1,1,-1],[-1,1,0],[-1,1,1],
                 [0,-1,-1],[0,-1,0],[0,-1,1],[0,0,-1],[0,0,1],[0,1,-1],[0,1,0],[0,1,1],
                 [1,-1,-1],[1,-1,0],[1,-1,1],[1,0,-1],[1,0,0],[1,0,1],[1,1,-1],[1,1,0],[1,1,1]])
w1 = 1
w2 = 1
def djikstra(map,start_index,goal_index):
    # Returns path
    heap = []
    #heap2 = heapdict()

    dist = {}
    visited = set()
    parent = {}
    dist[start_index] = 0

    heappush(heap,[0,start_index])
    #heap2[start_index] = 0

    goal_found = False
    nodes_expanded = 0


    #curr_idx,_ = heap2.popitem() 
    [d,curr_idx] = heappop(heap)

    visited.add(curr_idx)
    nodes_expanded+=1
    if goal_index==curr_idx:
        goal_found = True
    else:
       for n in neighbourhood26:
            # if n.tolist() ==d1.tolist():
            #     continue
            idx = tuple((np.array(curr_idx)+np.array(n)).tolist())
            #print(idx)
            #d = np.sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])
            d = np.abs(n).sum() 

            if idx not in visited:
                if not map.is_occupied_index(idx):
                    if not dist.get(idx):
                        dist[idx] = d+dist[curr_idx]
                        heappush(heap,[dist[idx],idx])
                        #heap2[idx] = dist[idx]
                        parent[idx] = curr_idx
                    elif ((d+dist[curr_idx])<dist[idx]):
                        dist[idx] = d+dist[curr_idx]
                        parent[idx] = curr_idx
                        #heap2[idx] = dist[idx]
                        heappush(heap,[dist[idx],idx])


    while heap:
        
        #curr_idx,_ = heap2.popitem() 
        [d,curr_idx] = heappop(heap)

        
        if goal_index==curr_idx:
            goal_found = True
            break
        #print(visited)
        if(curr_idx in visited):
            continue

        visited.add(curr_idx)
        nodes_expanded+=1
        d1 = np.clip(np.array(curr_idx)-np.array(parent[curr_idx]),-1,1)
        k = np.abs(d1)
        idx = tuple((np.array(curr_idx)+np.array(d1)).tolist())

        # if k.sum()==1 and not map.is_occupied_index(idx):
        #     idx = tuple((np.array(curr_idx)+np.array(d1)).tolist())
        #     d = np.sqrt(d1[0]*d1[0]+d1[1]*d1[1]+d1[2]*d1[2])
        #     neighs = [d1]
        #     # for a in range(2):
        #     #     if (map.is_occupied_index((idx[0]+int(k[0]==1)*a,idx[1]+int(k[1]==1)*a,idx[2]+int(k[2]==1)*a)):
        #     #         neighs.append(np.array([int(k[0]==1)*ad1[0]]))
        #     # if k[0] == 1:
        #     #     pass
        #     # elif k[1] == 1:
        #     # else:

        #     for n in neighs:
        #         idx = tuple((np.array(curr_idx)+np.array(n)).tolist())
        #         d = np.sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])
        #         if idx not in visited:
        #             if not map.is_occupied_index(idx):
        #                 if not dist.get(idx) or ((d+dist[curr_idx])<dist[idx]):
        #                     dist[idx] = d+dist[curr_idx]
        #                     heappush(heap,[dist[idx],idx])
        #                     parent[idx] = curr_idx
        # else:
        for n in neighbourhood26:
            if np.dot(n,d1)<=0:
                continue;
            idx = tuple((np.array(curr_idx)+np.array(n)).tolist())
            #print(idx)
            #d = np.sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])
            d = np.abs(n).sum() 

            if idx not in visited:
                if not map.is_occupied_index(idx):
                    if not dist.get(idx) or ((d+dist[curr_idx])<dist[idx]):
                        dist[idx] = d+dist[curr_idx]
                        heappush(heap,[dist[idx],idx])
                        parent[idx] = curr_idx


    if goal_found is True:
        return parse_path(map,parent,goal_index,start_index),nodes_expanded
    else:
        return None,nodes_expanded

    pass

def parse_path(map,parent,goal_index,start_index):
    idx = goal_index
    path = [map.index_to_metric_center(idx)]
    if idx not in parent.keys():
        return None
    while(parent[idx] != start_index):
        idx = parent[idx]
        path = [map.index_to_metric_center(idx)] + path
    path = [map.index_to_metric_center(start_index)]+path

    return path


def heuristic(goal_index,start_index):
    return np.abs(np.array(start_index)-np.array(goal_index)).max()

def heuristic2(goal_index,start_index):
    return np.sqrt((start_index[0]-goal_index[0])*(start_index[0]-goal_index[0])+(start_index[1]-goal_index[1])*(start_index[1]-goal_index[1])\
        +(start_index[2]-goal_index[2])*(start_index[2]-goal_index[2]))


def heuristic3(goal_index,start_index):
    return np.abs(np.array(start_index)-np.array(goal_index)).sum()/1.02

def heuristic4(goal_index,start_index):
    return np.abs(np.array(start_index)-np.array(goal_index)).sum() + 2*np.abs(np.array(start_index)-np.array(goal_index)).max()

def heuristic5(goal_index,start_index):
    k = np.abs(np.array(start_index)-np.array(goal_index))
    return k.sum() + 2*k.max()\
     + 1.0*np.sqrt(k[0]*k[0]+k[1]*k[1]+k[2]*k[2])

def a_star(map,start_index,goal_index):
   # Returns path
    heap = []

    #heap2 = heapdict()
    dist = {}
    dist[start_index] = 0
    visited = set()
    parent = {}


    #heap2[start_index] = heuristic2(goal_index,start_index)
    heappush(heap,[heuristic4(goal_index,start_index),start_index])

    [d,curr_idx] = heappop(heap)
        #curr_idx,_ = heap2.popitem()
    nodes_expanded = 0

    visited.add(curr_idx)
    nodes_expanded+=1
    

    goal_found = False
    if goal_index==curr_idx:
        goal_found = True
    else:
        for n in neighbourhood26:
            # if n.tolist() ==d1.tolist():
            #     continue

            idx = tuple((np.array(curr_idx)+np.array(n)).tolist())
            #d = np.linalg.norm(n)
            #d = np.abs(n).sum() + np.linalg.norm(n)
            #d = 1*np.abs(n).sum() + 2*1 + 1.0*np.sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])
            d = 1*np.abs(n).sum() + 2*1 
            if idx not in visited:
                if not map.is_occupied_index(idx):
                    if not dist.get(idx) or ((d+dist[curr_idx])<(dist[idx])):
                        dist[idx] = d+dist[curr_idx]
                        heappush(heap,[dist[idx] + heuristic4(idx,goal_index),idx])
                        parent[idx] = curr_idx


    while heap:
        [d,curr_idx] = heappop(heap)
        #curr_idx,_ = heap2.popitem()
        
        
        if goal_index==curr_idx:
            goal_found = True
            break
        #print(visited)
        if(curr_idx in visited):
            continue

        visited.add(curr_idx)
        nodes_expanded+=1

        d1 = np.clip(np.array(curr_idx)-np.array(parent[curr_idx]),-1,1)
        # dots = np.dot(neighbourhood26,d1)
        # p = np.argsort()
        # print(p)
        k = np.abs(d1).sum()
        idx = tuple((np.array(curr_idx)+np.array(d1)).tolist())
        # if k==1 and not map.is_occupied_index(idx):
        #     #print(idx)
        #     #d = 1*np.abs(d1).sum() + 2*1 + 1.0*np.sqrt(d1[0]*d1[0]+d1[1]*d1[1]+d1[2]*d1[2])
        #     d = 1*np.abs(d1).sum() + 2*1 
        #     #d = 1

        #     if idx not in visited:
        #         if not map.is_occupied_index(idx):
        #             if not dist.get(idx) or ((d+dist[curr_idx])<(dist[idx])):
        #                 dist[idx] = d+dist[curr_idx]
        #                 heappush(heap,[dist[idx] + heuristic4(idx,goal_index),idx])
        #                 #heap2[idx] = dist[idx] + heuristic2(idx,goal_index)
        #                 parent[idx] = curr_idx
        # else:
        for n in neighbourhood26:
            if np.dot(n,d1)<=0 :
                continue;
            idx = tuple((np.array(curr_idx)+np.array(n)).tolist())
            #d = np.linalg.norm(n)
            #d = np.abs(n).sum() + np.linalg.norm(n)
            #d = 1*np.abs(n).sum() + 2*1 + 1.0*np.sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2])
            d = 1*np.abs(n).sum() + 2*1 
            #d = 1

            if idx not in visited:
                if not map.is_occupied_index(idx):
                    if not dist.get(idx) or ((d+dist[curr_idx])<(dist[idx])):
                        dist[idx] = d+dist[curr_idx]
                        heappush(heap,[dist[idx] + heuristic4(idx,goal_index),idx])
                        #heap2[idx] = dist[idx] + heuristic2(idx,goal_index)
                        parent[idx] = curr_idx


    if goal_found is True:
        return (parse_path(map,parent,goal_index,start_index),nodes_expanded)
    else:
        return (None,nodes_expanded)



def graph_search(world, resolution, margin, start, goal, astar):
    """
    Parameters:
        world,      World object representing the environment obstacles
        resolution, xyz resolution in meters for an occupancy map, shape=(3,)
        margin,     minimum allowed distance in meters from path to obstacles.
        start,      xyz position in meters, shape=(3,)
        goal,       xyz position in meters, shape=(3,)
        astar,      if True use A*, else use Dijkstra
    Output:
        return a tuple (path, nodes_expanded)
        path,       xyz position coordinates along the path in meters with
                    shape=(N,3). These are typically the centers of visited
                    voxels of an occupancy map. The first point must be the
                    start and the last point must be the goal. If no path
                    exists, return None.
        nodes_expanded, the number of nodes that have been expanded
    """

    # While not required, we have provided an occupancy map you may use or modify.
    occ_map = OccupancyMap(world, resolution, margin)
    # Retrieve the index in the occupancy grid matrix corresponding to a position in space.
    start_index = tuple(occ_map.metric_to_index(start))
    goal_index = tuple(occ_map.metric_to_index(goal))

    if astar:
        path,nodes_expanded = a_star(occ_map,start_index,goal_index)
    else:
        path,nodes_expanded = djikstra(occ_map,start_index,goal_index)
    # Return a tuple (path, nodes_expanded)
    #print(path)
    if path is not None:
        path = [start] + path[1:-1] +[goal]
        return np.array(path),nodes_expanded
    else:
        return None,nodes_expanded

