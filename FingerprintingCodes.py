from bitstring import BitArray
import random
import math


class FingerprintingCode:

    def __init__(self, n, e, d=0):
        if d == 0:
            self.d = math.ceil(2 * math.pow(n, 2) * math.log((2 * n) / e))
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


class LogLengthCodes:

    def __init__(self, N, c, e):
        self.N = N
        self.L = math.ceil(2 * c * math.log((2 * N) / e))
        self.c = c
        self.n = 2 * c
        self.e = e
        self.d = math.ceil(2 * math.pow(self.n, 2) * math.log((4 * self.n * self.L) / e))
        self.codes = [FingerprintingCode(self.n, e, self.d) for i in range(self.L)]
        self.hiddenCode = self.createHiddenCode()
        self.codeBook = self.writeCodeBook()

    def createHiddenCode(self):
        h = []
        for i in range(self.N):
            h.append([random.randint(0, self.n - 1) for j in range(self.L)])
        return h

    def writeCodeBook(self):
        cb = []
        for h in self.hiddenCode:
            tmp = BitArray()
            for i in range(self.L):
                tmp.append(self.codes[i].codeBook[h[i]])
            cb.append(tmp)
        return cb

    def trace(self, x):
        l = self.d*(self.n-1)
        y = []
        for i in range(self.L):
            y.append(self.codes[i].trace(BitArray(bin=x.bin[(i*l):((i+1)*l)]))[0])
        index = 0
        matches = 0
        for i in range(self.N):
            count = 0
            for j in range(self.L):
                if self.hiddenCode[i][j] == y[j]:
                    count += 1
            if count > matches:
                index = i
                matches = count
        return index

    def collude(self, i, j):
        l = self.L * self.d * (self.n-1)
        c = BitArray(l)
        for k in range(l):
            c[k] = (self.codeBook[i][k] == 1 & self.codeBook[j][k] == 1)
        return c


if __name__ == '__main__':
    a = LogLengthCodes(10, 2, 0.1)
    print(a.trace(a.collude(2, 4)))
