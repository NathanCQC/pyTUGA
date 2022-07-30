from re import A
import numpy as np
import pandas as pd
from scipy.special import binom
import sys

sys.setrecursionlimit(10**8)
import quimb as qu
import quimb.tensor as qtn
from autoray import do


class ShavittGraph():
    """ Generates the Distinct row table for Shavitt GUGA graph
    follows convention - W. Dobrautz  -J. Chem. Phys. 151, 094104 (2019); doi: 10.1063/1.5108908 """
    def __init__(self, abc_triple: tuple):

        self._abc_final = abc_triple
        self._nspatorbs = self._abc_final[0] + self._abc_final[
            1] + self._abc_final[2]
        self._abc_node_dict = self._get_node_dict()
        self._ki_dict = self._get_downward_chaining_indices()
        self._li_dict = self._get_upward_chaining_indices()
        self._nnode = list(self._abc_node_dict.values())[-1]
        self._distinct_row_table = self._get_distinct_row_table()
        self._csf_list = []
        self._n_csfs = None

    def abc_dim(self):

        abc = np.array(self._abc_final)
        n = abc.sum()
        a = abc[0]
        b = abc[1]
        c = abc[2]
        dim = (b + 1) / (n + 1) * binom(n + 1, a) * binom(n + 1, c)
        return int(dim)

    class CSF():
        def __init__(self, idstat, jstat, index) -> None:
            self.index = index
            self.idstat = idstat.tolist()
            self.jstat = jstat.tolist()
            self.nspatorbs = len(self.idstat)
            self.arcs_upper = {
                j: d
                for j, d in zip(self.jstat[1:], self.idstat)
            }
            self.arcs_lower = {
                j: d
                for j, d in zip(self.jstat[0:self.nspatorbs], self.idstat)
            }

    def _get_node_dict(self) -> dict:

        # (a,b,c, N):j

        # Horrible indexing problem No other way I think. reminds me of first year phd suffering

        abc_node_dict = {}

        #Head Node
        node = 1
        abc_node_dict[(self._abc_final[0], self._abc_final[1],
                       self._abc_final[2], self._nspatorbs)] = node
        cmax = self._abc_final[2]

        # Each loop is one level at a time
        abc_left_above = self._abc_final
        for i, nspat in enumerate(reversed(range(self._nspatorbs))):
            i = i + 1
            if abc_left_above[2] > 0:
                left_abc = (abc_left_above[0], abc_left_above[1],
                            abc_left_above[2] - 1)
            elif (abc_left_above[2] == 0) and (abc_left_above[1] == 0):
                left_abc = (abc_left_above[0] - 1, abc_left_above[1],
                            abc_left_above[2])
            elif abc_left_above[2] == 0:
                left_abc = (abc_left_above[0], abc_left_above[1] - 1,
                            abc_left_above[2])
            else:
                raise (ValueError('Broken node (a,b,c) index'))
            # Set new left most node in orbital level

            # a groups
            a_left_abc = left_abc
            abc_left_above = left_abc
            amax = left_abc[0]
            if amax == 0:
                amin = 0
            else:
                if (amax - i) >= 0:
                    amin = amax - i
                else:
                    amin = 0

            for a in reversed(range(amin, amax + 1)):
                bmax = a_left_abc[1]  # left most node of a group
                bplusc = bmax + a_left_abc[
                    2]  # Constant for all nodes in in a agroup
                if (bmax - cmax) < 0:
                    bmin = 0
                else:
                    bmin = bmax - (cmax - a_left_abc[2]
                                   )  #This was a hack check it
                #bnodes in a group
                for b in reversed(range(bmin, bmax + 1)):
                    node += 1
                    c = bplusc - b
                    abc = (a, b, c, nspat)
                    abc_node_dict[abc] = node

                a_left_abc = (a - 1, a_left_abc[1] + 1, a_left_abc[2])

        return abc_node_dict

    def _get_downward_chaining_indices(self):

        ki_dict = {}
        for j, j_abc in zip(self._abc_node_dict.values(),
                            self._abc_node_dict.keys()):

            k0_abc = (j_abc[0], j_abc[1], j_abc[2] - 1, j_abc[3] - 1)
            k1_abc = (j_abc[0], j_abc[1] - 1, j_abc[2], j_abc[3] - 1)
            k2_abc = (j_abc[0] - 1, j_abc[1] + 1, j_abc[2] - 1, j_abc[3] - 1)
            k3_abc = (j_abc[0] - 1, j_abc[1], j_abc[2], j_abc[3] - 1)

            k0 = self._abc_node_dict.get(k0_abc, 'None')
            k1 = self._abc_node_dict.get(k1_abc, 'None')
            k2 = self._abc_node_dict.get(k2_abc, 'None')
            k3 = self._abc_node_dict.get(k3_abc, 'None')

            ki_dict[j] = [k0, k1, k2, k3]

        return ki_dict

    def _get_upward_chaining_indices(self):

        li_dict = {}
        for j, j_abc in zip(self._abc_node_dict.values(),
                            self._abc_node_dict.keys()):

            l0_abc = (j_abc[0], j_abc[1], j_abc[2] + 1, j_abc[3] + 1)
            l1_abc = (j_abc[0], j_abc[1] + 1, j_abc[2], j_abc[3] + 1)
            l2_abc = (j_abc[0] + 1, j_abc[1] - 1, j_abc[2] + 1, j_abc[3] + 1)
            l3_abc = (j_abc[0] + 1, j_abc[1], j_abc[2], j_abc[3] + 1)

            l0 = self._abc_node_dict.get(l0_abc, 'None')
            l1 = self._abc_node_dict.get(l1_abc, 'None')
            l2 = self._abc_node_dict.get(l2_abc, 'None')
            l3 = self._abc_node_dict.get(l3_abc, 'None')

            li_dict[j] = [l0, l1, l2, l3]

        return li_dict

    def _get_distinct_row_table(self):

        uirrep = []
        for u in reversed(range(self._nspatorbs + 1)):
            a = 1
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
            'k0':[i[0] for i in self._ki_dict.values()],  #downward chaining indices
            'k1': [i[1] for i in self._ki_dict.values()],
            'k2': [i[2] for i in self._ki_dict.values()],
            'k3': [i[3] for i in self._ki_dict.values()],
            'l0':[i[0] for i in self._li_dict.values()],  #upward chaining indices
            'l1': [i[1] for i in self._li_dict.values()],
            'l2': [i[2] for i in self._li_dict.values()],
            'l3': [i[3] for i in self._li_dict.values()],
        }

        return pd.DataFrame(data, index=nodes).rename_axis('j')

    @property
    def distinct_row_table(self):
        return self._distinct_row_table

    def _get_csfs(self):

        jstat = np.zeros(self._nspatorbs + 1, dtype=int)  # plus 1 bulshit
        # ibrnch = np.zeros(self._nspatorbs + 1, dtype=int)
        ibrnch = [
            None for i in range(self._nspatorbs)
        ]  #[(node, next d)] #records the branching at the progressive levels
        idstat = np.zeros(self._nspatorbs, dtype=int)
        nstate = 1

        jstat[0] = self._nnode

        nback_to = 0
        n = 0

        ibrnch[nback_to] = jstat[nback_to]

        if self.distinct_row_table['l3'][jstat[0]] is not None:
            idirn = 3
            idstat[0] = 3
            ibrnch[0] = 3
        elif self.distinct_row_table['l1'][jstat[0]] is not None:
            idirn = 1
            idstat[0] = 1
            ibrnch[0] = 1
        elif self.distinct_row_table['l0'][jstat[0]] is not None:
            idirn = 0
            idstat[0] = 0
            ibrnch[0] = None  # No branchig here.
        else:
            raise (ValueError('inital d error'))

        # Recursive walk Theory CSF

        # 1. input d direction, find next j with drt j3k/j2k/j1k/j0k

        # 2. If you 3 Not None take that or 2 , 1 , 0. Take the nexty lowest d after ibrach (j,d) that is not None

        # 3. If branching possible. More than 2 Not None. Take highest d, ibrnch = (j,d)

        # - How know if there is brnaching and which brnach to take

        # 4. If n equal to norbs repeat
        # 5. if n orbs =nspinobs for back to j branch and set i dir 0
        # 5. Go to last branch that was not 0

        # print('n',n,'idstat',idirn,'jstat',jstat[n])

        def recurive_walk(
                n, idirn, ibrnch, jstat,
                nstate):  # Does a single step, recurively updates for n+1
            # Needed to populate the
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
                raise (ValueError('d error'))

            def next_branch_d(ibrnch, n):

                branches = [
                    self.distinct_row_table['l0'][jstat[n]],
                    self.distinct_row_table['l1'][jstat[n]],
                    self.distinct_row_table['l2'][jstat[n]],
                    self.distinct_row_table['l3'][jstat[n]]
                ]

                if ibrnch[n] == (None or False):
                    ibrnch[n] = 4  #Hacky
                for d, ld in reversed(list(enumerate(branches[0:ibrnch[n]]))):
                    if ld != 'None':
                        if d != 0:
                            for d1, ld1 in reversed(
                                    list(enumerate(branches[0:d]))):
                                if ld1 != 'None':  # this checks for not None
                                    ibrnch[n] = d
                                    idirn = d
                                    return ibrnch, idirn
                                else:
                                    ibrnch[n] = False
                            idirn = d
                            return ibrnch, idirn
                        else:
                            ibrnch[n] = False
                        idirn = d
                        return ibrnch, idirn
                    # else: # THIS IS POSSIBLY BUGGY
                    #     raise ValueError('invalid branching direction')

            if (n != self._nspatorbs):
                ibrnch, idirn = next_branch_d(ibrnch, n)
                return recurive_walk(n, idirn, ibrnch, jstat, nstate)
            else:  # rewind to the last brnach step
                csf = self.CSF(idstat, jstat, nstate)
                self._csf_list.append(csf)
                for nback in range(
                        self._nspatorbs - 1, nback_to - 1, -1
                ):  #This is going all the way back to 0 when it shouldt really. Should go to last branch. Seems ot work aslong as put to 0
                    if (ibrnch[nback] != False):
                        idirn = ibrnch[
                            nback]  #Step direction changes only when ibrnch[nback]!=0
                        n = nback
                        ibrnch, idirn = next_branch_d(ibrnch, n)
                        nstate = nstate + 1
                        return recurive_walk(n, idirn, ibrnch, jstat, nstate)
                return  #End of function

        recurive_walk(n, idirn, ibrnch, jstat,
                      nstate)  #Entry point to the function

        self._nstates = nstate

        return

    def tensor_network_state(self, type):
        #k are downward chaining indices
        j_tensors = []

        # Build the node tensors
        for j in self.distinct_row_table.index.tolist():
            # tagged [ j , u , uirrep] Maybe this could just be linked to the row of drt.
            # If ever in doubt just link index j to drt
            u = self.distinct_row_table['u'][j]
            u = f'U{u}'
            uirrep = self.distinct_row_table['uirrep'][j]
            uirrep = f'UIRREP{uirrep}'

            j_tensor = qtn.Tensor(
                tags={f'J{j}{type}', u, uirrep,
                      type.upper()})

            if type == 'bra':
                j_tensor.new_ind(f'{j}bra', size=4)
            elif type == 'ket':
                j_tensor.new_ind(f'{j}ket', size=4)
            else:
                raise ValueError('incorrect type keyword must be bra or ket')

            for l in range(4):
                if self.distinct_row_table[f'l{l}'][j] != 'None':
                    l_node = self.distinct_row_table[f'l{l}'][j]
                    j_tensor.new_ind(f'{l_node}_{j}', size=4)

            for k in range(4):
                if self.distinct_row_table[f'k{k}'][j] != 'None':
                    k_node = self.distinct_row_table[f'k{k}'][j]
                    j_tensor.new_ind(f'{j}_{k_node}', size=4)

            # Need to add in the paraeterisatio here
            # Need to pilt the J tensor to be lambda - d Jlambda'

            print(j_tensor)

            j_tensors.append(j_tensor)

        # Build tensor network
        tensor_network_state = j_tensors[0]
        for j in j_tensors[1:]:
            tensor_network_state = j & tensor_network_state

        if type == 'bra':
            tensor_network_state = tensor_network_state.H

        return tensor_network_state


    def ptensor_network_state(self): # Gate ride of type as we hope Tensor.H Handles this
        #k are downward chaining indices
        j_tensors = []

        def elements(params):
            def matrix_slice(data):
                c = do('cos', data[0] / 2)
                s = do('sin', data[0] / 2)
                e = do('exp', 1j*data[1])
                d = [[1, 0, 0, 0], [0, c, -s, 0], [0, s, c, 0], [0, 0, 0, e]]
                return d
            data = [matrix_slice(i) for i in params]
            return do('array', data, like=params)

        # Build the node tensors
        for j in self.distinct_row_table.index.tolist():
            # tagged [ j , u , uirrep] Maybe this could just be linked to the row of drt.
            # If ever in doubt just link index j to drt
            u = self.distinct_row_table['u'][j]
            u = f'U{u}'
            uirrep = self.distinct_row_table['uirrep'][j]
            uirrep = f'UIRREP{uirrep}'

            inds_l = ['l']
            for l in range(4):
                if self.distinct_row_table[f'l{l}'][j] != 'None':
                    l_node = self.distinct_row_table[f'l{l}'][j]
                    # j_tensor.new_ind(f'{l_node}_{j}', size=4)
                    inds_l.append(f'{l_node}_{j}')

            j_l_tensor = qtn.COPY_tensor(4,inds=inds_l) # I thinks this is 4 not len(inds_l)-1
            print(j_l_tensor.shape)

            # Upper Walks
            inds_k = ['k']
            for k in range(4):
                if self.distinct_row_table[f'k{k}'][j] != 'None':
                    k_node = self.distinct_row_table[f'k{k}'][j]
                    # j_tensor.new_ind(f'{j}_{k_node}', size=4)
                    inds_k.append(f'{k_node}_{j}')

            j_k_tensor =  qtn.COPY_tensor(4,inds=inds_k)
            print(j_k_tensor.shape)

            #Parameterised tensor
            j_ind = f'{j}'

            angles = [(0,0) for i in range(len(inds_k))] #This are the initial parameters
            if (len(inds_k) != 0) and (len(inds_l) != 0):
                print('angles',angles) 
                d_tensor = qtn.PTensor(elements, angles, inds=['l','k',j_ind])
            elif len(inds_l) == 0:
                d_tensor = qtn.PTensor(elements, angles, inds=['k',j_ind])
            elif len(inds_k) == 0:
                d_tensor = qtn.PTensor(elements, angles, inds=['l',j_ind])
            else:
                raise(ValueError('BUG ALERT'))

            j_tensor = j_l_tensor & d_tensor & j_k_tensor

            j_tensors.append(j_tensor)

        # Build tensor network
        tensor_network_state = j_tensors[0]
        for j in j_tensors[1:]:
            tensor_network_state = j & tensor_network_state

        if type == 'bra':
            tensor_network_state = tensor_network_state.H

        return tensor_network_state

    def tensor_network_operator(
            self):  #Need to figure our what goes in the tensor
        #k are downward chaining indices
        j_tensors = []

        # Build the node tensors
        for j in self.distinct_row_table.index.tolist():
            # tagged [ j , u , uirrep] Maybe this could just be linked to the row of drt.
            # If ever in doubt just link index j to drt
            u = self.distinct_row_table['u'][j]
            u = f'U{u}'
            uirrep = self.distinct_row_table['uirrep'][j]
            uirrep = f'UIRREP{uirrep}'

            j_tensor = qtn.Tensor(tags={f'J{j}', u, uirrep, 'OPERATOR'})

            j_tensor.new_ind(f'{j}bra', size=4)
            j_tensor.new_ind(f'{j}ket', size=4)

            for l in range(4):
                if self.distinct_row_table[f'l{l}'][j] != 'None':
                    l_node = self.distinct_row_table[f'l{l}'][j]
                    j_tensor.new_ind(f'{l_node}_{j}', size=4)

            for k in range(4):
                if self.distinct_row_table[f'k{k}'][j] != 'None':
                    k_node = self.distinct_row_table[f'k{k}'][j]
                    j_tensor.new_ind(f'{j}_{k_node}', size=4)

            j_tensors.append(j_tensor)

        # Build tensor network
        tensor_network_operator = j_tensors[0]
        for j in j_tensors[1:]:
            tensor_network_operator = j & tensor_network_operator

        return tensor_network_operator

    def tensor_network_expection(self):

        bra = self.tensor_network_state('bra')
        operator = self.tensor_network_operator()
        ket = self.tensor_network_state('ket')

        tensor_expectation = bra & operator & ket

        return tensor_expectation



class ShavittGraphAB():
    """ Generates the Distinct row table for Shavitt GUGA graph
    follows convention - W. Dobrautz  -J. Chem. Phys. 151, 094104 (2019); doi: 10.1063/1.5108908 """
    def __init__(self, abc_triple: tuple):

        if abc_triple[0] != abc_triple[2] and abc_triple[2] == 0:
            raise ValueError('abc triple must be (n,0,n)')

        self._abc_final = abc_triple
        self._nspatorbs = self._abc_final[0] + self._abc_final[
            1] + self._abc_final[2]
        self._abc_node_dict = self._get_node_dict()
        self._ki_dict = self._get_downward_chaining_indices()
        self._li_dict = self._get_upward_chaining_indices()
        self._nnode = list(self._abc_node_dict.values())[-1]
        self._distinct_row_table = self._get_distinct_row_table()
        self._csf_list = []
        self._n_csfs = None

    def abc_dim(self):

        abc = np.array(self._abc_final)
        n = abc.sum()
        a = abc[0]
        b = abc[1]
        c = abc[2]
        dim = (b + 1) / (n + 1) * binom(n + 1, a) 
        return int(dim)

    class CSF():
        def __init__(self, idstat, jstat, index) -> None:
            self.index = index
            self.idstat = idstat.tolist()
            self.jstat = jstat.tolist()
            self.nspatorbs = len(self.idstat)
            self.arcs_upper = {
                j: d
                for j, d in zip(self.jstat[1:], self.idstat)
            }
            self.arcs_lower = {
                j: d
                for j, d in zip(self.jstat[0:self.nspatorbs], self.idstat)
            }

    def _get_node_dict(self) -> dict:

        # This will only work for (n,0,n)

        # (a,b,c, N):j

        abc_node_dict = {}

        #Head Node
        node = 1
        abc_node_dict[(self._abc_final[0], self._abc_final[1],
                       self._abc_final[2], self._nspatorbs)] = node
        cmax = self._abc_final[2]

        # Each loop is one level at a time
        abc_left_above = self._abc_final
        for i, nspat in enumerate(reversed(range(self._nspatorbs))):
            # i = i + 1
            if abc_left_above[1] == 1:
                left_abc = (abc_left_above[0], 0, abc_left_above[0])
            elif abc_left_above[1] == 0:
                left_abc = (abc_left_above[0] - 1, 1, (abc_left_above[0] - 1))
            else:
                raise (ValueError('Broken node (a,b,c) index'))
            # Set new left most node in orbital level


            # a groups
            a_left_abc = left_abc
            abc_left_above = left_abc
            amax = left_abc[0]
            if amax == 0:
                amin = 0
            else:
                if (amax - i) >= 0:
                    amin = amax - i
                else:
                    amin = 0
            print('oi',nspat,i, amax,amin)

            c = left_abc[2]
            for a in reversed(range(amin, amax + 1)):
                print('a',a)
                node += 1
                b = nspat - a - c
                print(a,b,c,nspat)
                abc = (a, b ,c, nspat)
                abc_node_dict[abc] = node
                c = c - 1
            
        return abc_node_dict

    def _get_downward_chaining_indices(self):

        ki_dict = {}
        for j, j_abc in zip(self._abc_node_dict.values(),
                            self._abc_node_dict.keys()):

            # k0_abc = (j_abc[0], j_abc[1], j_abc[2] - 1, j_abc[3] - 1)
            k1_abc = (j_abc[0], j_abc[1] - 1, j_abc[2], j_abc[3] - 1)
            k2_abc = (j_abc[0] - 1, j_abc[1] + 1, j_abc[2] - 1, j_abc[3] - 1)
            # k3_abc = (j_abc[0] - 1, j_abc[1], j_abc[2], j_abc[3] - 1)

            # k0 = self._abc_node_dict.get(k0_abc, 'None')
            k1 = self._abc_node_dict.get(k1_abc, 'None')
            k2 = self._abc_node_dict.get(k2_abc, 'None')
            # k3 = self._abc_node_dict.get(k3_abc, 'None')

            ki_dict[j] = [k1, k2]

        return ki_dict

    def _get_upward_chaining_indices(self):

        li_dict = {}
        for j, j_abc in zip(self._abc_node_dict.values(),
                            self._abc_node_dict.keys()):

            # l0_abc = (j_abc[0], j_abc[1], j_abc[2] + 1, j_abc[3] + 1)
            l1_abc = (j_abc[0], j_abc[1] + 1, j_abc[2], j_abc[3] + 1)
            l2_abc = (j_abc[0] + 1, j_abc[1] - 1, j_abc[2] + 1, j_abc[3] + 1)
            # l3_abc = (j_abc[0] + 1, j_abc[1], j_abc[2], j_abc[3] + 1)

            # l0 = self._abc_node_dict.get(l0_abc, 'None')
            l1 = self._abc_node_dict.get(l1_abc, 'None')
            l2 = self._abc_node_dict.get(l2_abc, 'None')
            # l3 = self._abc_node_dict.get(l3_abc, 'None')

            li_dict[j] = [l1, l2]

        return li_dict

    def _get_distinct_row_table(self):

        uirrep = []
        for u in reversed(range(self._nspatorbs + 1)):
            a = 1
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
            # 'k0':[i[0] for i in self._ki_dict.values()],  #downward chaining indices
            'k1': [i[0] for i in self._ki_dict.values()],
            'k2': [i[1] for i in self._ki_dict.values()],
            # 'k3': [i[3] for i in self._ki_dict.values()],
            # 'l0':[i[0] for i in self._li_dict.values()],  #upward chaining indices
            'l1': [i[0] for i in self._li_dict.values()],
            'l2': [i[1] for i in self._li_dict.values()],
            # 'l3': [i[3] for i in self._li_dict.values()],
        }

        return pd.DataFrame(data, index=nodes).rename_axis('j')

    @property
    def distinct_row_table(self):
        return self._distinct_row_table

    def _get_csfs(self):

        jstat = np.zeros(self._nspatorbs + 1, dtype=int)  # plus 1 bulshit
        # ibrnch = np.zeros(self._nspatorbs + 1, dtype=int)
        ibrnch = [
            None for i in range(self._nspatorbs)
        ]  #[(node, next d)] #records the branching at the progressive levels
        idstat = np.zeros(self._nspatorbs, dtype=int)
        nstate = 1

        jstat[0] = self._nnode

        nback_to = 0
        n = 0

        ibrnch[nback_to] = jstat[nback_to]

        # if self.distinct_row_table['l3'][jstat[0]] is not None:
        #     idirn = 3
        #     idstat[0] = 3
        #     ibrnch[0] = 3
        if self.distinct_row_table['l1'][jstat[0]] is not None:
            idirn = 1
            idstat[0] = 1
            ibrnch[0] = 1
        # elif self.distinct_row_table['l0'][jstat[0]] is not None:
        #     idirn = 0
        #     idstat[0] = 0
        #     ibrnch[0] = None  # No branchig here.
        else:
            raise (ValueError('inital d error'))

        # Recursive walk Theory CSF

        # 1. input d direction, find next j with drt j3k/j2k/j1k/j0k

        # 2. If you 3 Not None take that or 2 , 1 , 0. Take the nexty lowest d after ibrach (j,d) that is not None

        # 3. If branching possible. More than 2 Not None. Take highest d, ibrnch = (j,d)

        # - How know if there is brnaching and which brnach to take

        # 4. If n equal to norbs repeat
        # 5. if n orbs =nspinobs for back to j branch and set i dir 0
        # 5. Go to last branch that was not 0

        # print('n',n,'idstat',idirn,'jstat',jstat[n])

        def recurive_walk(
                n, idirn, ibrnch, jstat,
                nstate):  # Does a single step, recurively updates for n+1
            # Needed to populate the
            n = n + 1  # Next level from previous brnaching point
            idstat[n - 1] = idirn

            # if (idirn == 1):
            #     jstat[n] = self.distinct_row_table['j1k'][jstat[n - 1]]
            # else:
            #     jstat[n] = self.distinct_row_table['j0k'][jstat[n - 1]]

            # if (idirn == 3):
            #     jstat[n] = self.distinct_row_table['l3'][jstat[n - 1]]
            if (idirn == 2):
                jstat[n] = self.distinct_row_table['l2'][jstat[n - 1]]
            elif (idirn == 1):
                jstat[n] = self.distinct_row_table['l1'][jstat[n - 1]]
            # elif (idirn == 0):
            #     jstat[n] = self.distinct_row_table['l0'][jstat[n - 1]]
            else:
                raise (ValueError('d error'))

            def next_branch_d(ibrnch, n):

                branches = [
                    # self.distinct_row_table['l0'][jstat[n]],
                    self.distinct_row_table['l1'][jstat[n]],
                    self.distinct_row_table['l2'][jstat[n]]
                    # self.distinct_row_table['l3'][jstat[n]]
                ]

                if ibrnch[n] == (None or False):
                    ibrnch[n] = 4  #Hacky
                for d, ld in reversed(list(enumerate(branches[0:ibrnch[n]]))):
                    if ld != 'None':
                        if d != 0:
                            for d1, ld1 in reversed(
                                    list(enumerate(branches[0:d]))):
                                if ld1 != 'None':  # this checks for not None
                                    ibrnch[n] = d
                                    idirn = d
                                    return ibrnch, idirn
                                else:
                                    ibrnch[n] = False
                            idirn = d
                            return ibrnch, idirn
                        else:
                            ibrnch[n] = False
                        idirn = d
                        return ibrnch, idirn
                    # else: # THIS IS POSSIBLY BUGGY
                    #     raise ValueError('invalid branching direction')

            if (n != self._nspatorbs):
                ibrnch, idirn = next_branch_d(ibrnch, n)
                return recurive_walk(n, idirn, ibrnch, jstat, nstate)
            else:  # rewind to the last brnach step
                csf = self.CSF(idstat, jstat, nstate)
                self._csf_list.append(csf)
                for nback in range(
                        self._nspatorbs - 1, nback_to - 1, -1
                ):  #This is going all the way back to 0 when it shouldt really. Should go to last branch. Seems ot work aslong as put to 0
                    if (ibrnch[nback] != False):
                        idirn = ibrnch[
                            nback]  #Step direction changes only when ibrnch[nback]!=0
                        n = nback
                        ibrnch, idirn = next_branch_d(ibrnch, n)
                        nstate = nstate + 1
                        return recurive_walk(n, idirn, ibrnch, jstat, nstate)
                return  #End of function

        recurive_walk(n, idirn, ibrnch, jstat,
                      nstate)  #Entry point to the function

        self._n_csfs  = nstate

        return

    def tensor_network_state(self, type):
        #k are downward chaining indices
        j_tensors = []

        # Build the node tensors
        for j in self.distinct_row_table.index.tolist():
            # tagged [ j , u , uirrep] Maybe this could just be linked to the row of drt.
            # If ever in doubt just link index j to drt
            u = self.distinct_row_table['u'][j]
            u = f'U{u}'
            uirrep = self.distinct_row_table['uirrep'][j]
            uirrep = f'UIRREP{uirrep}'

            j_tensor = qtn.Tensor(
                tags={f'J{j}{type}', u, uirrep,
                      type.upper()})

            if type == 'bra':
                j_tensor.new_ind(f'{j}bra', size=2)
            elif type == 'ket':
                j_tensor.new_ind(f'{j}ket', size=2)
            else:
                raise ValueError('incorrect type keyword must be bra or ket')

            for l in [1,2]:
                if self.distinct_row_table[f'l{l}'][j] != 'None':
                    l_node = self.distinct_row_table[f'l{l}'][j]
                    j_tensor.new_ind(f'{l_node}_{j}', size=2)

            for k in [1,2]:
                if self.distinct_row_table[f'k{k}'][j] != 'None':
                    k_node = self.distinct_row_table[f'k{k}'][j]
                    j_tensor.new_ind(f'{j}_{k_node}', size=2)

            # Need to add in the paraeterisatio here
            # Need to pilt the J tensor to be lambda - d Jlambda'

            print(j_tensor)

            # def j_ptensor(params, in_dim, out_dim):
            #     phi,psi = params
            #     c = do('cos', phi / 2)
            #     s = do('sin', phi / 2)
            #     e = do('exp', 1j*psi )
            #     data = [[1,0,0,0],[0,c, -s,0], [0,s, c,0],[0,0,0,e]]
            #     return do('array', data, like=params)

            # def j_ptensor(params. in_dim, out_dim):

            # # to do this the tensors need to be sqwashed to rank 3

            # for l in range(4):
            #     if self.distinct_row_table[f'l{l}'][j] != 'None':
            #         for k in range(4):
            #             if self.distinct_row_table[f'k{k}'][j] != 'None':
            #                 k_node = self.distinct_row_table[f'k{k}'][j]
            #                 j_tensor.new_ind(f'{j}_{k_node}', size=4)

            #     for p in params:
            #         c = do('cos', p[0] / 2)
            #         s = do('sin', p[0] / 2)
            #         e = do('exp', 1j*p[1] )
            #         data = [[1,0,0,0],[0,c, -s,0], [0,s, c,0],[0,0,0,e]]
            #     return do('array', data, like=params)

            # qtn.PTensor(j_ptensor, [0,3.14], inds=(4, 4))

            j_tensors.append(j_tensor)

        # Build tensor network
        tensor_network_state = j_tensors[0]
        for j in j_tensors[1:]:
            tensor_network_state = j & tensor_network_state

        if type == 'bra':
            tensor_network_state = tensor_network_state.H

        return tensor_network_state


    def ptensor_network_state(self,type='bra'): # Gate ride of type as we hope Tensor.H Handles this
        #k are downward chaining indices
        j_tensors = []



        # Build the node tensors
        for j in reversed(self.distinct_row_table.index.tolist()):
            # tagged [ j , u , uirrep] Maybe this could just be linked to the row of drt.
            # If ever in doubt just link index j to drt
            u = self.distinct_row_table['u'][j]
            u = f'U{u}'
            uirrep = self.distinct_row_table['uirrep'][j]
            uirrep = f'UIRREP{uirrep}'

            # Upper Walks (Very confusing l is upper)
            inds_l = [] # This is quite sketchky think of better way
            if j != self.distinct_row_table.index.tolist()[0]:
                inds_l.append(f'{j}l')
                for l in [1,2]:
                    if self.distinct_row_table[f'l{l}'][j] != 'None':
                        print('oij')
                        l_node = self.distinct_row_table[f'l{l}'][j]
                        inds_l.append(f'{l_node}_{j}')
                j_l_tensor = qtn.COPY_tensor(2,inds=inds_l) # I thinks this is 4 nodt len(inds_l)-1

            # Lower Walks
            inds_k = []
            if j != self.distinct_row_table.index.tolist()[-1]:
                inds_k.append(f'{j}k')
                for k in [1,2]:
                    if self.distinct_row_table[f'k{k}'][j] != 'None':
                        k_node = self.distinct_row_table[f'k{k}'][j]
                        # j_tensor.new_ind(f'{j}_{k_node}', size=4)
                        inds_k.append(f'{j}_{k_node}')
                j_k_tensor =  qtn.COPY_tensor(2,inds=inds_k)

            dim_walk_ind_d = 0 # Must be a better way of doing this.
            for k in [1,2]:
                if self.distinct_row_table[f'k{k}'][j] != 'None':
                    dim_walk_ind_d += 1 

            #Parameterised tensor
            j_ind = f'{j}'

            def elements(params):
                def matrix_slice(theta):
                    c = do('cos', theta / 2)
                    s = do('sin', theta / 2)
                    d = ((c, -s,), (s, c))
                    return d
                data = [matrix_slice(i) for i in params]
                return do('array', data)


            def elements_row(params): # Do we need a parameter here?, Does it need to be normalised? Like in tensor networks code
                def matrix_slice(theta):
                    c = do('cos', theta / 2)
                    s = do('sin', theta / 2)
                    d = (c, -s,)
                    return d
                data = [matrix_slice(i) for i in params]
                return do('array', data)


             #This are the initial parameters
            if inds_k != [] and inds_l != []: # Middle nodes
                print('head')
                angles = [0,0]
                d_tensor = qtn.PTensor(elements, angles, inds=[f'{j}l',f'{j}k',j_ind])
                j_tensor = j_l_tensor & d_tensor & j_k_tensor
                # rank 2 array
            elif inds_l == []: # No upper walks. Head node. head node should be a column. But not sure if this is applicable in tensor networks as long as the indexing is correct
                # Row vector
                angles = [0]
                d_tensor = qtn.PTensor(elements_row, angles, inds=[f'{j}k',j_ind]) # Head node should be a column, BUG IT SEEMS TO BE NOT CONNECTING
                j_tensor = d_tensor & j_k_tensor
            elif inds_k == []: # No upper walks. Tail node
                # column vector
                angles = [0]
                d_tensor = qtn.PTensor(elements_row, angles, inds=[f'{j}l',j_ind])
                j_tensor = j_l_tensor & d_tensor 
            else:
                raise(ValueError('BUG ALERT'))

            print(j,j_tensor)

            j_tensors.append(j_tensor)

        # Build tensor network
        tensor_network_state = j_tensors[0]
        for j in j_tensors[1:]:
            tensor_network_state = j & tensor_network_state

        if type == 'bra':
            tensor_network_state = tensor_network_state.H

        return tensor_network_state

    def tensor_network_operator(
            self):  #Need to figure our what goes in the tensor
        #k are downward chaining indices
        j_tensors = []

        # Build the node tensors
        for j in self.distinct_row_table.index.tolist():
            # tagged [ j , u , uirrep] Maybe this could just be linked to the row of drt.
            # If ever in doubt just link index j to drt
            u = self.distinct_row_table['u'][j]
            u = f'U{u}'
            uirrep = self.distinct_row_table['uirrep'][j]
            uirrep = f'UIRREP{uirrep}'

            j_tensor = qtn.Tensor(tags={f'J{j}', u, uirrep, 'OPERATOR'})

            j_tensor.new_ind(f'{j}bra', size=4)
            j_tensor.new_ind(f'{j}ket', size=4)

            for l in [1,2]:
                if self.distinct_row_table[f'l{l}'][j] != 'None':
                    l_node = self.distinct_row_table[f'l{l}'][j]
                    j_tensor.new_ind(f'{l_node}_{j}', size=4)

            for k in [1,2]:
                if self.distinct_row_table[f'k{k}'][j] != 'None':
                    k_node = self.distinct_row_table[f'k{k}'][j]
                    j_tensor.new_ind(f'{j}_{k_node}', size=4)

            j_tensors.append(j_tensor)

        # Build tensor network
        tensor_network_operator = j_tensors[0]
        for j in j_tensors[1:]:
            tensor_network_operator = j & tensor_network_operator

        return tensor_network_operator

        # Build tensor network
        tensor_network_operator = j_tensors[0]
        for j in j_tensors[1:]:
            tensor_network_operator = j & tensor_network_operator

        return tensor_network_operator

    def tensor_network_expection(self):

        bra = self.tensor_network_state('bra')
        operator = self.tensor_network_operator()
        ket = self.tensor_network_state('ket')

        tensor_expectation = bra & operator & ket

        return tensor_expectation





# 2 Mai problems

# If we are going for the full contraction approach I dont know how to contai the integral factors int the operator
# If we are going to the loop driven contraction this is easy, but require more programming

if __name__ == "__main__":

    abc = (3, 0, 3)
    guga = ShavittGraph(abc)
    print(guga.distinct_row_table)