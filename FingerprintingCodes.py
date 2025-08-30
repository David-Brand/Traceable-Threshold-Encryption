from bitstring import BitArray
import random
import math


class FingerprintingCode:

    def __init__(self, n, e, d=0):
        if d == 0:
            self.d = math.floor(2 * math.pow(n, 2) * math.log((2 * n) / e))
        else:
            self.d = d
        self.n = n
        self.e = e
        self.l = self.d * (n - 1)
        self.permutation = list(range(self.l))
        random.shuffle(self.permutation)
        self.codeBook = []
        self.writeCodeBook()

    def writeCodeBook(self):
        for i in range(0, self.n):
            tmp = BitArray(self.l)
            tmp.set(1, range(i * self.d, self.l))

            word = BitArray(self.l)
            for j in range(self.l):
                word[self.permutation[j]] = tmp[j]

            self.codeBook.append(word)

    def trace(self, x):
        guilty = []
        if self.weight(x, range(self.d)) > 0:
            guilty.append(0)

        for i in range(1, self.n - 1):
            k = self.weight(x, range((i - 1) * self.d, (i + 1) * self.d))
            val = (k / 2) - math.sqrt((k / 2) * math.log((2 * self.n) / self.e))
            if self.weight(x, range((i - 1) * self.d, i * self.d)) < val:
                guilty.append(i)

        if self.weight(x, range(self.l - self.d, self.l)) < self.d:
            guilty.append(self.n - 1)
        return guilty

    def weight(self, x, r):
        tmp = 0
        for i in r:
            tmp += x[self.permutation[i]]
        return tmp

    def out(self):
        for word in self.codeBook:
            print(word.bin)

    def collude(self, i, j):
        c = BitArray(self.l)
        for k in range(self.l):
            c[k] = (self.codeBook[i][k] == 1 & self.codeBook[j][k] == 1)
        return c


if __name__ == '__main__':
    a = FingerprintingCode(4, 0.2)
    a.out()
    print(a.trace(a.collude(1, 2)))
