import numpy as np

"""
A clustering algorithm is used to identify the strips corresponding to charges that came from the
electron cloud. The algorithm first looks for a strip whose signal is above a defined threshold,
0.08 AU in this work, to find a strip that definitely contains a signal. Then, starting from this
seed strip, the algorithm scans in both direction to construct a cluster. The scan stops when it
finds a strip that collected a negative charge, which means just noise. So we are going to find
the positions n1 and n2. Which n1 represents the first strip of the cluster and n2 the last one
"""


class Cluster:
    def __init__(self, seed: float, threshold: float, Charges: list, Positions: list):
        self.seed = seed
        self.Q = Charges
        self.P = Positions
        self.zero = threshold

    def find_seed(self) -> int:
        i = 0
        Qi = self.Q[i]
        while Qi < self.seed:
            i = i + 1
            Qi = self.Q[i]
        seed_index = i
        return seed_index

    def scans_left(self, n1: int) -> int:
        Qleft = self.Q[n1]
        while Qleft > self.zero and n1 != 0:
            n1 = n1 - 1
            Qleft = self.Q[n1]

        if Qleft > self.zero:
            return n1
        else:
            n1 = n1 + 1
            return n1

    def scans_right(self, n2: int) -> int:
        Qright = self.Q[n2]
        last_index = len(self.Q) - 1
        while Qright > self.zero and n2 < last_index:
            n2 = n2 + 1
            Qright = self.Q[n2]
        if Qright > self.zero:
            return n2
        else:
            n2 = n2 - 1
            return n2

    def Find_Cluster(self) -> tuple:
        if len([k for k in self.Q if k > self.seed]) == 0:
            raise NameError("All strips collected less charge than the seed value")
        seed_index = self.find_seed()
        n1 = seed_index
        n2 = seed_index
        last_index = len(self.Q) - 1

        if seed_index > 0 and seed_index < last_index:
            n1 = self.scans_left(seed_index - 1)
            n2 = self.scans_right(seed_index + 1)
        elif seed_index == 0:
            n1 = 0
            n2 = self.scans_right(seed_index + 1)
        elif seed_index == last_index:
            n2 = last_index
            n1 = self.scans_left(n2)

        if n1 != n2:
            n2 = n2 + 1
            P_Cluster = self.P[n1:n2]
            Q_Cluster = self.Q[n1:n2]
        else:
            P_Cluster = np.zeros(1)
            Q_Cluster = np.zeros(1)
            P_Cluster[0] = self.P[n1]
            Q_Cluster[0] = self.Q[n1]

        return Q_Cluster, P_Cluster
