import numpy as np
from calculus import Calculus
import time
import copy
from functools import reduce

def save_times(filename,times):
    with open(filename,'w') as file:
        for line in times:
            for ele in line:
                file.write(str(ele)+",")
            file.write("\n")

def init_1(buffer,calc):
    # create spatial variables of specified number
    temp = buffer.pop(0)
    Vars = int(temp.split()[0].strip())
    TypeId = temp.split()[1].strip()
    DALL = calc.B_dict['DALL']
    # initialize constraints
    Id = calc.Id
    cons_set=[]
  
    ConMatrix = [None for v in range(Vars)]
    for v in range(Vars):
        j = {}   
        j[v] = Id
        ConMatrix[v] = j

    # parse spatial CSP
    for line in buffer:
        l = line.strip().replace('(','').replace(')','').split()
        if l == ['.']:
            break
        s = reduce(lambda x, y: x | y, [calc.translateR(i) for i in l[2:]])
        cons_set.append((int(l[0]),s,int(l[1])))
        ConMatrix[int(l[0])][int(l[1])] = DALL
        ConMatrix[int(l[1])][int(l[0])] = DALL
    
    np.random.shuffle(cons_set)
    
    return TypeId,ConMatrix,cons_set
def directional_path_consistency_incremental_general(con_matrix,neighbors,pair,new_constraints,calc):                                                                                                                                                                                   
    P = dict()
    a,a_ = pair
    Q = np.zeros(len(a),dtype=int)
    num_var = len(a_)
    DALL = calc.B_dict['DALL']
    for vi,rij_,vj in new_constraints:
        if vj not in neighbors[vi]:
            if rij_ != DALL:
                # no effect on relations but structure of graph
                con_matrix[vi][vj]=rij_
                con_matrix[vj][vi]=calc.inverse(rij_)
                neighbors[vi].add(vj)
                neighbors[vj].add(vi)
                # ensure vj has higher order
                if a[vi] > a[vj]:
                    temp = vi
                    vi = vj
                    vj = temp
                if vj not in P:
                    P[vj] = {vi}
                else:
                    P[vj].add(vi)
                Q[a[vj]] = 1
        else:
            rijo = con_matrix[vi][vj]

            if rij_ == rijo:
                continue
            if (rij_ & rijo) != rij_:
                nbrs = [vk_ for vk_ in neighbors[vi] if vk_ in neighbors[vj] and a[vk_] > a[vi] and a[vk_]>a[vj] and (con_matrix[vi][vk_] != DALL or con_matrix[vk_][vj] != DALL)]
                for vk_ in nbrs:
                    rij_ = rij_ & calc.comp(con_matrix[vi][vk_],con_matrix[vk_][vj])
                    if not rij_:
                        return None
            if rij_ != rijo:
                con_matrix[vi][vj]=rij_
                con_matrix[vj][vi]=calc.inverse(rij_)
                if a[vi] > a[vj]:
                    temp = vi
                    vi = vj
                    vj = temp
                if vj not in P:
                    P[vj] = {vi}
                else:
                    P[vj].add(vi)
                Q[a[vj]] = 1
    
    p = num_var - 1
    while p != -1:
        if Q[p] == 1:
            vk = a_[p]
            neighbs = [y for y in neighbors[vk] if a[y] < p and con_matrix[y][vk] != DALL]
            n_nbrs = len(neighbs)

            for idxi in range(n_nbrs):
                for idxj in range(idxi+1,n_nbrs):
                    vi = neighbs[idxi]
                    vj = neighbs[idxj]
                    if vi in P[vk] or vj in P[vk]:
                        if a[vi] > a[vj]:
                            temp = vi
                            vi = vj
                            vj = temp
                        if vj not in neighbors[vi]:
                            con_matrix[vi][vj] = DALL
                            con_matrix[vj][vi] = DALL
                            neighbors[vi].add(vj)
                            neighbors[vj].add(vi)
                        tij = con_matrix[vi][vj] & calc.comp(con_matrix[vi][vk],con_matrix[vk][vj])
                        if not tij:
                            return None
                        if tij != con_matrix[vi][vj]:
                            con_matrix[vi][vj]=tij
                            con_matrix[vj][vi]=calc.inverse(tij)
                        Q[a[vj]] = 1
                        if vj not in P:
                            P[vj] = {vi}
                        else:
                            P[vj].add(vi) 
        p = p - 1
    return con_matrix
def add_all_newconstraints_experiment(filename,data_name,calc):
    times = []
    with open(filename) as f:
        total = 0
        meter = 0
        totalMCS = 0
        buffer = []

        for i in f:
            buffer.append(i)

            if i.strip() == '.':

                start = elapsed = elapsedMCS = 0
                meter += 1

                TypeId, ConMatrix,cons_set = init_1(buffer,calc)
                buffer = []

                from triangulation import MCS, FIC
                
                edjes = set([])

                neighbors = tuple([set([]) for i in range(len(ConMatrix))])
                for n, i in enumerate(ConMatrix):
                    for j in i:
                        if j > n:
                            edjes.add((n,j))
                            neighbors[n].add(j)
                            neighbors[j].add(n)

                # variable ordering computation with MCS
                start = time.perf_counter()
                a,a_ = MCS(ConMatrix, neighbors)
                elapsedMCS = time.perf_counter() - start

                fill = FIC(ConMatrix, neighbors, (a,a_))
                DALL = calc.B_dict['DALL']

                for i,j in fill:
                    if (i,j) not in edjes:
                        ConMatrix[i][j] = DALL
                        ConMatrix[j][i] = DALL
                        neighbors[i].add(j)
                        neighbors[j].add(i)

                # test incremental dpc with new constraints, adding them all together
                ConMatrix_temp = copy.deepcopy(ConMatrix)
                neighbors_temp = copy.deepcopy(neighbors)
                a_temp = copy.deepcopy(a)
                a_temp_ = copy.deepcopy(a_)
                start = time.perf_counter()
                solution = directional_path_consistency_incremental_general(ConMatrix_temp,neighbors_temp,(a_temp,a_temp_),cons_set,calc)
                elapsed = time.perf_counter() - start
                total += elapsed
                if solution is not None:
                    sat = 1
                else:
                    sat = 0
                times.append(["dpcia",total,elapsedMCS,elapsed,sat,0])


    save_times("../data/results/"+data_name+"_times_dpcia_temp.csv",times)
def add_constraints_back_1_by_1_experiment(filename,data_name,calc):
    times = []
    with open(filename) as f:
        total = 0
        meter = 0
        totalMCS = 0
        buffer = []

        for i in f:
            buffer.append(i)

            if i.strip() == '.':

                start = elapsed = elapsedMCS = 0
                meter += 1

                TypeId, ConMatrix,cons_set = init_1(buffer,calc)
                
                
                buffer = []

                from triangulation import MCS, FIC
                
                edjes = set([])

                neighbors = tuple([set([]) for i in range(len(ConMatrix))])
                for n, i in enumerate(ConMatrix):
                    for j in i:
                        if j > n:
                            edjes.add((n,j))
                            neighbors[n].add(j)
                            neighbors[j].add(n)

                # variable ordering computation with MCS
                start = time.perf_counter()
                a,a_ = MCS(ConMatrix, neighbors)
                elapsedMCS = time.perf_counter() - start
                fill = FIC(ConMatrix, neighbors, (a,a_))
                
                DALL = calc.B_dict['DALL']

                for i,j in fill:
                    if (i,j) not in edjes:
                        ConMatrix[i][j] = DALL
                        ConMatrix[j][i] = DALL
                        neighbors[i].add(j)
                        neighbors[j].add(i)

                # test incremental dpc with new constraints, adding them one by one
                ConMatrix_temp = copy.deepcopy(ConMatrix)
                neighbors_temp = copy.deepcopy(neighbors)
                a_temp = copy.deepcopy(a)
                a_temp_ = copy.deepcopy(a_)
                start = time.perf_counter()
                for cons in cons_set:
                    new_constraints_temp = [cons]
                    solution = directional_path_consistency_incremental_general(ConMatrix_temp,neighbors_temp,(a_temp,a_temp_),new_constraints_temp,calc)
                    if solution is None:
                        break
                elapsed = time.perf_counter() - start
                total += elapsed
                if solution is not None:
                    sat = 1
                else:
                    sat = 0
                times.append(["dpci1",total,elapsedMCS,elapsed,sat,0])

    save_times("../data/results/"+data_name+"_times_dpci1_temp.csv",times)


if __name__ == '__main__':
    calc = Calculus('rcc8','full')
    filename = "../data/aij/Footprint_1.csp"
    data_name= filename.split("/")[-1].split(".")[0]
    print(data_name)
    experiment = "1"
    if experiment== "all":
        add_all_newconstraints_experiment(filename,data_name,calc)
    elif experiment=="1":
        add_constraints_back_1_by_1_experiment(filename,data_name,calc)