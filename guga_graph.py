from re import A
import numpy as np
import pandas as pd
from sympy import *

class ShavittGraph():
    """ Generates the Distinct row table for Shavitt GUGA graph
    follows convention - W. Dobrautz  -J. Chem. Phys. 151, 094104 (2019); doi: 10.1063/1.5108908 """
    def __init__(self, abc_triple:tuple):


        self._abc_final = abc_triple
        self._nspatorbs = self._abc_final[0] + self._abc_final[1] + self._abc_final[2]
        self._abc_node_dict = self._get_node_dict()
        self._ki_dict = self._get_downward_chaining_indices()
        self._li_dict = self._get_upward_chaining_indices()


    class CSF(): 

        def __init__(self,index) -> None:
            self.index = index
            self.arc = None
            self.arcs_top = None

            pass

    def _get_node_dict(self) -> dict:

        # (a,b,c, N):j

        # Horrible indexing problem No other way I think. reminds me of first year phd suffering

        abc_node_dict = {}

        #Head Node
        node = 1
        abc_node_dict[(self._abc_final[0],self._abc_final[1],self._abc_final[2],self._nspatorbs)] = node
        cmax = self._abc_final[2]

        # Each loop is one level at a time
        abc_left_above = self._abc_final
        for i,nspat in enumerate(reversed(range(self._nspatorbs))):
            i = i+1
            if abc_left_above[2] > 0:
                left_abc = (abc_left_above[0],abc_left_above[1],abc_left_above[2]-1)
            elif (abc_left_above[2] ==0) and (abc_left_above[1] == 0):
                left_abc = (abc_left_above[0]-1,abc_left_above[1],abc_left_above[2])
            elif abc_left_above[2] == 0:
                left_abc = (abc_left_above[0],abc_left_above[1]-1,abc_left_above[2])
            else:
                raise(ValueError('Broken node (a,b,c) index'))
            # Set new left most node in orbital level

            # a groups
            a_left_abc = left_abc
            abc_left_above = left_abc
            amax = left_abc[0]
            if amax == 0:
                amin=0
            else:
                if (amax - i) >= 0:
                    amin = amax - i
                else:
                    amin = 0

            for a in reversed(range(amin, amax+1)):
                bmax = a_left_abc[1] # left most node of a group
                bplusc = bmax + a_left_abc[2] # Constant for all nodes in in a agroup
                if (bmax - cmax) < 0:
                    bmin = 0
                else:
                    bmin = bmax - (cmax - a_left_abc[2]) #This was a hack check it
                #bnodes in a group
                for b in reversed(range(bmin,bmax+1)):
                    node +=1
                    c = bplusc - b
                    abc = (a,b,c,nspat)
                    abc_node_dict[abc] = node

                a_left_abc = (a - 1, a_left_abc[1]+1, a_left_abc[2])

        return abc_node_dict


    def _get_downward_chaining_indices(self):

        ki_dict = {}
        for j,j_abc in zip(self._abc_node_dict.values(),self._abc_node_dict.keys()):

            k0_abc = (j_abc[0], j_abc[1], j_abc[2]-1,j_abc[3]-1)
            k1_abc = (j_abc[0], j_abc[1]-1, j_abc[2],j_abc[3]-1)
            k2_abc = (j_abc[0]-1, j_abc[1]+1, j_abc[2]-1,j_abc[3]-1)
            k3_abc = (j_abc[0]-1, j_abc[1], j_abc[2],j_abc[3]-1)

            k0 = self._abc_node_dict.get(k0_abc, 'None')
            k1 = self._abc_node_dict.get(k1_abc, 'None')
            k2 = self._abc_node_dict.get(k2_abc, 'None')
            k3 = self._abc_node_dict.get(k3_abc, 'None')

            ki_dict[j]= [k0,k1,k2,k3]

        return ki_dict


    def _get_upward_chaining_indices(self):

        li_dict = {}
        for j,j_abc in zip(self._abc_node_dict.values(),self._abc_node_dict.keys()):

            l0_abc = (j_abc[0], j_abc[1], j_abc[2] + 1, j_abc[3] + 1)
            l1_abc = (j_abc[0], j_abc[1] + 1, j_abc[2], j_abc[3] + 1)
            l2_abc = (j_abc[0] + 1, j_abc[1] - 1, j_abc[2] + 1, j_abc[3] + 1)
            l3_abc = (j_abc[0] + 1, j_abc[1], j_abc[2] , j_abc[3] + 1)

            l0 = self._abc_node_dict.get(l0_abc, 'None')
            l1 = self._abc_node_dict.get(l1_abc, 'None')
            l2 = self._abc_node_dict.get(l2_abc, 'None')
            l3 = self._abc_node_dict.get(l3_abc, 'None')

            li_dict[j]= [l0,l1,l2,l3]

        return li_dict



    @property
    def distinct_row_table(self):

        uirrep = []
        for u in reversed(range(self._nspatorbs+1)):
            a=1
            for abcN in self._abc_node_dict.keys():
                if abcN[3] == u:
                    uirrep.append(a)
                    a = a + 1

        nodes = [i for i in self._abc_node_dict.values()]
        data = {
            'a': [i[0] for i in self._abc_node_dict.keys()],
            'b': [i[1] for i in self._abc_node_dict.keys()],
            'c': [i[2] for i in self._abc_node_dict.keys()],
            'u': [i[3] for i in self._abc_node_dict.keys()],
            'uirrep': uirrep,
            'k0': [i[0] for i in self._ki_dict.values()],  #downward chaining indices
            'k1': [i[1] for i in self._ki_dict.values()],
            'k2': [i[2] for i in self._ki_dict.values()],
            'k3': [i[3] for i in self._ki_dict.values()],
            'l0': [i[0] for i in self._li_dict.values()],  #upward chaining indices
            'l1': [i[1] for i in self._li_dict.values()],
            'l2': [i[2] for i in self._li_dict.values()],
            'l3': [i[3] for i in self._li_dict.values()],
        }

        return pd.DataFrame(data,index=nodes).rename_axis('j')

    def _get_csfs(self):

        # Do we actually need to ever run this code if we have the nodes and the connectivity from the DRT

        #Things we need
        # - CSF object for each walk
        # ibrnch - list stored branching point opn current path at orbtial level N
        # jstat - list of nodes of csf path
        # d_vec - walk list
        # d

        # Lexicographical ordering = 3 - 2 - 1 - 0

        # Do we need to use numpy arrays

        jstat = np.zeros(self._nspatorbs + 1, dtype=int) # plus 1 bulshit
        ibrnch = np.zeros(self._nspatorbs + 1, dtype=int)
        idstat = np.zeros(self._nspatorbs, dtype=int)

        jstat[0] = self._nnode

        nback_to = 0

        ibrnch[nback_to] = jstat[nback_to]
        # idstat[nback_to] = 1  #first branch direction = 1
        idstat[nback_to] = 3
        idirn = 3  #First branching direction 1

        def recurive_walk(n, idirn):  # Does a single step, recurively updates for n+1
            n = n + 1  # Next level from previous brnaching point
            idstat[n - 1] = idirn
            # if (idirn == 1):
            #     jstat[n] = self.distinct_row_table['j1k'][jstat[n - 1]]
            # else:
            #     jstat[n] = self.distinct_row_table['j0k'][jstat[n - 1]]

            if (idirn == 3):
                jstat[n] = self.distinct_row_table['l3'][jstat[n - 1]]
            elif (idirn == 2):
                jstat[n] = self.distinct_row_table['l2'][jstat[n - 1]]
            elif (idirn == 1):
                jstat[n] = self.distinct_row_table['l1'][jstat[n - 1]]
            elif (idirn == 0):
                jstat[n] = self.distinct_row_table['l0'][jstat[n - 1]]
            else:
                raise(ValueError('d error'))

            #branching decide direction of next step
            # This needs more complex logic now that we have 4 step diretcion
            if ((self.distinct_row_table['j1k'][jstat[n]] != 0) and (self.distinct_row_table['j0k'][jstat[n]] != 0)):
                ibrnch[n] = jstat[n]
                idirn = 1

            if all(v is not None for v in [A, B, C, D, E])

            branches = [self.distinct_row_table['l3'][jstat[n]],self.distinct_row_table['l2'][jstat[n]],self.distinct_row_table['l1'][jstat[n]],self.distinct_row_table['l0'][jstat[n]]]

            

            #NO BRANCHING
            else:
                ibrnch[n] = 0
                if (self.distinct_row_table['j1k'][jstat[n]] != 0):
                    idirn = 1
                else:
                    idirn = 0

            if (n != self._nspinorbs):
                return recurive_walk(n, idirn)
            else:  # rewind to the last brnach step
                for nback in range(self._nspinorbs - 1, nback_to - 1, -1):
                    if (ibrnch[nback] != 0):
                        ibrnch[nback] = 0  #Correctly updates the ibrach
                        idirn = 0  #Step direction changes only when ibrnch[nback]!=0
                        n = nback
                        # nstate=nstate+1
                        slater_determinants.append(idstat.tolist()) # Do this in CSF class
                        slater_determinants_arcs.append([(j, d)for j, d in zip(jstat[:self._nspinorbs].tolist(),idstat.tolist())])
                        slater_determinants_arcs_top.append([(j, d) for j, d in zip(jstat[1:].tolist(),idstat.tolist())])
                        return recurive_walk(n, idirn)
                slater_determinants.append(idstat.tolist())
                slater_determinants_arcs.append([(j, d) for j, d in zip(jstat[:self._nspinorbs].tolist(), idstat.tolist())])
                slater_determinants_arcs_top.append([(j, d) for j, d in zip(jstat[1:].tolist(), idstat.tolist()) #this happens in csf object])
                return slater_determinants, slater_determinants_arcs, slater_determinants_arcs_top

        slater_determinants, slater_determinants_arcs, slater_determinants_arcs_top = recurive_walk(n, idirn)

        return slater_determinants, slater_determinants_arcs, slater_determinants_arcs_top  #return DF here



    def _get_slater_determiant_paths(self):

        jstat = np.zeros(self._nspinorbs + 1, dtype=int)
        ibrnch = np.zeros(self._nspinorbs + 1, dtype=int)
        idstat = np.zeros(self._nspinorbs, dtype=int)
        slater_determinants = []
        slater_determinants_arcs = []
        slater_determinants_arcs_top = []

        jstat[0] = self._nnode
        nstate = 0

        n = 0  #Base of graph
        nback_to = n  #branching point at the base

        ibrnch[nback_to] = jstat[nback_to]
        idstat[nback_to] = 1  #first branch direction = 1
        idirn = 1  #First branching direction 1

        def recurive_walk(
                n, idirn):  # Does a single step, recurively updates for n+1
            n = n + 1  # Next level from previous brnaching point
            idstat[n - 1] = idirn
            if (idirn == 1):
                jstat[n] = self.distinct_row_table['j1k'][jstat[n - 1]]
            else:
                jstat[n] = self.distinct_row_table['j0k'][jstat[n - 1]]

            #branching decide direction of next step
            if ((self.distinct_row_table['j1k'][jstat[n]] != 0)
                    and (self.distinct_row_table['j0k'][jstat[n]] != 0)):
                ibrnch[n] = jstat[n]
                idirn = 1

            #NO BRANCHING
            else:
                ibrnch[n] = 0
                if (self.distinct_row_table['j1k'][jstat[n]] != 0):
                    idirn = 1
                else:
                    idirn = 0

            if (n != self._nspinorbs):
                return recurive_walk(n, idirn)
            else:  # rewind to the last brnach step
                for nback in range(self._nspinorbs - 1, nback_to - 1, -1):
                    if (ibrnch[nback] != 0):
                        ibrnch[nback] = 0  #Correctly updates the ibrach
                        idirn = 0  #Step direction changes only when ibrnch[nback]!=0
                        n = nback
                        # nstate=nstate+1
                        slater_determinants.append(idstat.tolist())
                        slater_determinants_arcs.append([
                            (j, d)
                            for j, d in zip(jstat[:self._nspinorbs].tolist(),
                                            idstat.tolist())
                        ])
                        slater_determinants_arcs_top.append([
                            (j, d) for j, d in zip(jstat[1:].tolist(),
                                                   idstat.tolist())
                        ])
                        return recurive_walk(n, idirn)
                slater_determinants.append(idstat.tolist())
                slater_determinants_arcs.append([(j, d) for j, d in zip(
                    jstat[:self._nspinorbs].tolist(), idstat.tolist())])
                slater_determinants_arcs_top.append([
                    (j, d) for j, d in zip(jstat[1:].tolist(), idstat.tolist()) #this happens in csf object
                ])
                return slater_determinants, slater_determinants_arcs, slater_determinants_arcs_top

        slater_determinants, slater_determinants_arcs, slater_determinants_arcs_top = recurive_walk(n, idirn)

        return slater_determinants, slater_determinants_arcs, slater_determinants_arcs_top  #return DF here
