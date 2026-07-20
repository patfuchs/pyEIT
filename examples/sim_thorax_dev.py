# coding: utf-8
""" Creation of thorax pfix with different electrode configurations"""
# Author: Patricia Fuchs
from __future__ import absolute_import, division, print_function

import matplotlib.pyplot as plt
import numpy as np
import pyEITdev.pyeit.mesh as mesh
from pyEITdev.pyeit.mesh.shape import thorax
from pyEITdev.pyeit.mesh.plot_mesh import plot_mesh
from pyEITdev.pyeit.mesh.wrapper import PyEITAnomaly_Circle

# -------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------- #

def plot_mesh2(mesh_obj, edgenodes, figsize: tuple = (6, 4), newpoints=[], newedges=[], title: str = "mesh"):
    """
    Plot a PyEITMesh

    Parameters
    ----------
    mesh_obj : PyEITMesh

    """
    plt.style.use("default")
    pts = mesh_obj.node
    tri = mesh_obj.element
    x, y = pts[:, 0], pts[:, 1]
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111)
    ax.tripcolor(
        x,
        y,
        tri,
        np.real(mesh_obj.perm_array),
        edgecolors="k",
        shading="flat",
        alpha=0.5,
        cmap=plt.cm.viridis,
    )
    # draw electrodes
    ax.plot(x[mesh_obj.el_pos], y[mesh_obj.el_pos], "ro")
    for i, e in enumerate(mesh_obj.el_pos):
        ax.text(x[e], y[e], str(i + 1), size=12)
    ax.plot(x[edgenodes], y[edgenodes], "bo")
   # for k in range(len(edgenodes)):
   #     ax.text(x[edgenodes[k]], y[edgenodes[k]], str(edgenodes[k]), size=8)
    if len(newedges) > 0:
        ax.plot(x[newedges], y[newedges], "go")
    #   for k in range(len(newedges)):
    #      ax.text(x[newedges[k]], y[newedges[k]], str(newedges[k]), size=8)
    if len(newpoints) > 0:
        ax.plot(newpoints[:, 0], newpoints[:, 1], "go")
    ax.set_title(title)
    ax.set_aspect("equal")
    ax.set_ylim([-1.2, 1.2])
    ax.set_xlim([-1.2, 1.2])
    fig.set_size_inches(6, 6)
    plt.show()


# -------------------------------------------------------------------------------------------------------------------- #
# Work on elems
def get_all_connecting_elements(edgeidx, elements):
    elems = []
    for k in elements:
        if k[0] == edgeidx or k[1] == edgeidx or k[2] == edgeidx:
            elems.append(k)
    return elems

def get_all_neighbors(edgenodes, elems, bRemovenonedges=True):
    neighbors = []
    for i in range(len(edgenodes)):
        neigh_elems = np.unique(np.array(get_all_connecting_elements(edgenodes[i], elems)).flatten())
        neigh_elems = np.delete(neigh_elems, np.where(neigh_elems == edgenodes[i])[0])
        neighbors.append(neigh_elems)
    if bRemovenonedges:
        for i in range(len(neighbors)):
            mask = np.ones(len(neighbors[i]), dtype=bool)
            for nindex, n in enumerate(neighbors[i]):
                if n in edgenodes:
                    foundNeighbor = True
                else:
                    mask[nindex] = False
            neighbors[i] = neighbors[i][mask]
    return neighbors



#todo delete
def check_n(edges, n: list, group, group_idx, edge_idx):
    n_idx = 0
    if edge_idx in group:
        return group, group_idx
    while edges[n_idx] != edge_idx:
        n_idx += 1
    if n_idx in group_idx:
        return group, group_idx
    elif len(n[n_idx]) >= 3:
        group.append(edge_idx)
        group_idx.append(n_idx)
        for n2 in range(len(n[n_idx])):
            group, group_idx = check_n(edges, n, group, group_idx, n[n_idx][n2])
        return group, group_idx
    else:
        return group, group_idx

# -------------------------------------------------------------------------------------------------------------------- #
class MeshForm:
    def __init__(self, elem):
        self.elem = list(elem)
        self.iCorners = len(elem)
        self.edges = self._make_edges()
        self.edges_reverse = self._make_edges(bReverse=True)

    def add_corner(self, index, corner):

        self.elem = self.elem[:index] + [corner] + self.elem[index:]
        self.iCorners = len(self.elem)
        self.edges = self._make_edges()
        self.edges_reverse = self._make_edges(bReverse=True)

    def inside_mesh(self, other):
        for e in other.elem:
            if e not in self.elem:
                return False
        return True

    def shares_edge(self, other):
        shared = False
        insert_idx = 0
        for e in other.edges:
            if e in self.edges or e in self.edges_reverse:
                insert_idx = 0
                shared = True
                for i in range(len(self.edges)):
                    if self.edges[i] == e or self.edges_reverse[i] == e:
                        insert_idx = i + 1

                return shared, insert_idx
        return shared, insert_idx

    def get_corners_no_x(self, x):
        if type(x) != list and type(x) != np.array:
            x = [x]
        elem_copy = self.elem.copy()
        for k in x:
            if k in elem_copy:
                elem_copy.remove(k)
        return elem_copy

    def _make_edges(self, bReverse=False):
        edges = []
        direction = 1 if not bReverse else -1
        if not bReverse:
            for i in range(0, self.iCorners):
                edges.append((self.elem[i % self.iCorners], self.elem[(i + 1) % self.iCorners]))
        else:
            for i in range(0, self.iCorners):
                edges.append((self.elem[(i + 1) % self.iCorners], self.elem[(i) % self.iCorners]))

        return edges


# -------------------------------------------------------------------------------------------------------------------- #
class MeshMatcher:
    def __init__(self, elements, nodes):
        self.elements = list(elements)
        self.nodes = nodes
        self.iN = len(nodes)
        self.iE = len(elements)
        self.posEdges = []
        self.piEdgeNodes = []

    def get_possible_edge_nodes(self):
        pass

    def test_enclosed_elem(self, iEdgeIdx):
        connecting_elems = get_all_connecting_elements(iEdgeIdx, self.elements)
        conn_mesh_elems = [MeshForm(connecting_elems[i]) for i in range(len(connecting_elems))]
        corners = [iEdgeIdx] + conn_mesh_elems[0].get_corners_no_x(iEdgeIdx)
        corners = MeshForm(corners)
        done = False
        i = 1
        fullyclosed = False
        k = 0
        while not done:
            e = conn_mesh_elems[i]
            if corners.inside_mesh(e):
                fullyclosed = True
                break
            shared_edges, insert_idx = corners.shares_edge(e)
            if not shared_edges:  # not fitting piece
                conn_mesh_elems = conn_mesh_elems[:i] + conn_mesh_elems[i + 1:] + [conn_mesh_elems[i]]
            else:
                to_add = e.get_corners_no_x(corners.elem)
                for a in to_add:
                    corners.add_corner(insert_idx, a)
                i += 1
            if i >= len(conn_mesh_elems):
                done = True
            k += 1
            if k == 10:
                print("loop")
                break

        if fullyclosed:
            return False
        else:
            return True

    def test_all_enclosed_elem(self, edgenodes):
        true_edges = []
        for e in edgenodes:
            edge = self.test_enclosed_elem(e)
            if edge == True:
                true_edges.append(e)
        self.piEdgeNodes = true_edges
        return true_edges



# -------------------------------------------------------------------------------------------------------------------- #
def get_edgenodes_by_neighbors(nodes, elements, firstouter, deltamin, thneighbors, thelems):
    def get_smallest_dist(n1, nodelist):
        norm = np.linalg.norm(nodelist - n1, ord=2, axis=1)
        norm[norm == 0] = 10  # do not show distance to one self
        k = np.argmin(norm)
        mindist = norm[k]
        return mindist

    def get_neighbors(point, nodelist, distance):
        norm = np.linalg.norm(nodelist - point, ord=2, axis=1)
        norm[norm == 0] = 10  # do not show distance to one self
        index = np.zeros(len(nodelist))
        index[norm < distance] = 1
        idx = [i for i, x in enumerate(index) if x == 1]
        return np.sum(index), idx

    def count_neigboring_elements(nodeidx, elements):
        z = np.zeros(len(elements))
        for e, elem in enumerate(elements):
            z[e] = 1 if ( nodeidx == elem[0] or nodeidx == elem[1] or nodeidx == elem[2]) else 0
        return np.sum(z)


    edgenodes = []
    edgelen = (1 + deltamin) * get_smallest_dist(firstouter, nodes)

    for j in range(len(nodes)):
        nNeighbors, idx_neighbors = get_neighbors(nodes[j], nodes, edgelen)
        n_elems = count_neigboring_elements(nodes[j], elements)
        if nNeighbors < thneighbors and n_elems < thelems:
            edgenodes.append(j)
    return edgenodes

# -------------------------------------------------------------------------------------------------------------------- #
""" Main function for getting the indexes of all edge nodes of a meshgrid"""
def get_index_all_edge_nodes(nodes, elements, firstouter, thelems=5, deltamin=0.35, thneighbors=5, mesh=None):
   # Get initial guess on edgenodes based on the amount of neighbors and neighboring elements
    edgenodes = get_edgenodes_by_neighbors(nodes, elements, firstouter, deltamin, thneighbors, thelems)

    edgenodes = np.array(edgenodes)
    if mesh is not None:
        plot_mesh2(mesh, edgenodes, title="Initial Edge Node Guess")
    edgenodes = list(edgenodes)
    # Get actual neighbors
    neighbors = get_all_neighbors(edgenodes, elements,bRemovenonedges=True)

    # Remove points without neighbors
    to_del = []
    for i in range(len(edgenodes)):
        if len(neighbors[i]) == 0:
            to_del.append(i)
    for d in reversed(sorted(to_del)):
        del edgenodes[d]
        del neighbors[d]

    # Find sorrounding elements for edge detection
    m1 = MeshMatcher(elements, nodes)
    true_edges = m1.test_all_enclosed_elem(edgenodes)
    true_neighbors = get_all_neighbors(true_edges, elements, True)

    edgenodes = np.array(edgenodes)
    plot_mesh2(mesh, edgenodes, title="Final Edge Node Picks", newedges=np.array(true_edges))
    edgenodes = list(edgenodes)

    # Order in consecutive order
    new_order = [true_edges[0]]
    curridx = 0
    for i in range(1, len(true_edges)):
        k = 0
        while true_neighbors[curridx][k] not in true_edges or true_neighbors[curridx][k] in new_order:
            k += 1
            if k == len(true_neighbors[curridx]):
                raise ValueError()
        # normal case
        new_order.append(true_neighbors[curridx][k])
        z = 0
        while true_edges[z] != true_neighbors[curridx][k]:
            z += 1

        curridx = z
    true_edges = list(reversed(new_order))
    true_edges = [true_edges[-1]] + true_edges[:-1]
    return np.array(true_edges)


# -------------------------------------------------------------------------------------------------------------------- #
""" Create new pfix arrays """
def get_pfix_doubled_electrodes(edgenodesidx, p16, nodes, nel=16):
    new_pfix = []
    for i in range(nel):
        idx_curr, idx_foll = p16[i], p16[(i + 1) % nel]
        idx_new = find_new_idx(edgenodesidx, idx_curr)
        new_pfix.append(nodes[idx_curr])
        idx_new2 = find_new_idx(edgenodesidx, idx_foll)
        if idx_new2 == 0:
            idx_new = idx_new + 1
            idx_new2 = len(edgenodesidx) - 1
        dist = np.abs(idx_new2 - idx_new)
        if dist % 2 == 1:
            new_pfix.append(np.round(nodes[edgenodesidx[(idx_new + idx_new2) // 2]], 4))
        else:
            k = (idx_new + idx_new2) // 2
            pnew = (nodes[edgenodesidx[k]] + nodes[edgenodesidx[k + 1]]) / 2
            new_pfix.append(np.round(pnew, 4))
    return np.array(new_pfix)


def find_new_idx(l, value):
    i = 0
    while l[i] != value:
        i += 1
    return i


# -------------------------------------------------------------------------------------------------------------------- #
# -------------------------------------------------------------------------------------------------------------------- #
""" Build Mesh """
# Mesh shape is specified with fd parameter in the instantiation, e.g : fd=thorax , Default :fd=circle
n_el = 16  # nb of electrodes
mesh2 = mesh.create(16, fd=thorax, h0=0.04, shift=8)
plot_mesh(mesh2)
anmly = PyEITAnomaly_Circle([0.2, -0.1], 5, 0.08)
alung1 = PyEITAnomaly_Circle([-0.5, -0.1], 0.0004, 0.25)
alung2 = PyEITAnomaly_Circle([-0.5, 0.1], 0.0004, 0.25)
alung3 = PyEITAnomaly_Circle([0.5, -0.1], 0.0004, 0.25)
alung4 = PyEITAnomaly_Circle([0.5, 0.1], 0.0004, 0.25)
m = mesh.set_perm(mesh2, [anmly, alung1, alung2, alung3, alung4], background=0.7)
plot_mesh(m)

# -------------------------------------------------------------------------------------------------------------------- #
""" Get Pfix for Thorax16 with electrodes centered """
edgenodes = get_index_all_edge_nodes(m.node, m.element, m.node[m.el_pos[0]], mesh=m)
pfix16 = m.node[m.el_pos]
pfix16mid = get_pfix_doubled_electrodes(edgenodes, m.el_pos, m.node)
pfix16mid = pfix16mid[::2, :2]
#pfixnew = get_16_electrodes_mid(edgenodes, m.el_pos, m.node)
plot_mesh2(m, edgenodes, newpoints=pfix16mid)
#pfixnew = pfixnew[:, :2]
mesh3 = mesh.create(16, p_fix=pfix16mid, fd=thorax, h0=0.04)
plot_mesh(mesh3)

""" Get Pfix Thorax32 """
pfix32 = get_pfix_doubled_electrodes(edgenodes, m.el_pos, m.node)
pfix32 = pfix32[:, :2]

mesh4 = mesh.create(32, p_fix=pfix32, fd=thorax, h0=0.04)  # Create Mesh Grid for 32 Electrodes
plot_mesh(mesh4)

# ------------------------------------------------ #
""" Get Pfix Thorax64 """
edgenodesnew = get_index_all_edge_nodes(mesh4.node, mesh4.element, mesh4.node[mesh4.el_pos[0]], Th_elems=6,
                                        Thneighbors=6, deltamin=0.35, mesh=mesh4)
pfix32idx = mesh4.el_pos
pfix64 = get_pfix_doubled_electrodes(edgenodesnew, pfix32idx, mesh4.node, nel=32)
plot_mesh2(m, edgenodes, newpoints=pfix64)

# ------------------------------------------------ #
""" Save pfix Arrays"""
save_path = "../../pyEITdev/pyeit/mesh/thorax_pfix/"
bSave = True
if bSave:
    np.save(save_path + "pfix32_electrodes.npy", pfix32)
    np.save(save_path + "pfix16_electrodes.npy", pfix16)
    np.save(save_path + "pfix64_electrodes.npy", pfix64)
    np.save(save_path + "pfix16mid_electrodes.npy", pfix16mid)
