#!/usr/bin/env python
# -*- coding: utf-8 -*-

from .tools_karkkainen_sanders import direct_kark_sort
from array import array
import numpy as np


class Rstr_max:
    def __init__(self):
        self.array_str = []

    def add_str(self, seq):
        self.array_str.append(seq)

    def step1_sort_suffix(self):
        self.global_suffix = self.array_str[0]
        for i in range(1, len(self.array_str)):
            self.global_suffix += [chr(2)] + self.array_str[i]
        nbChars = len(self.global_suffix)
        init = [-1]*nbChars
        self.idxString = array('i', init)
        self.idxPos = array('i', init)
        self.endAt = array('i', init)
        k = idx = 0
        for mot in self.array_str:
            last = k + len(mot)
            for p in range(len(mot)):
                self.idxString[k] = idx
                self.idxPos[k] = p
                self.endAt[k] = last
                k += 1
            idx += 1
            k += 1

        self.res = direct_kark_sort(self.global_suffix)

    def step2_lcp(self):
        n = len(self.res)
        init = [0]*n
        rank = array('i', init)
        LCP = array('i', init)

        s = self.global_suffix
        suffix_array = self.res
        endAt = self.endAt

        for i in range(len(self.array_str), n):
            v = self.res[i]
            rank[v] = i
        l = 0
        for j in range(n):
            if(l > 0):
                l -= 1
            i = rank[j]
            j2 = suffix_array[i-1]
            if i:
                while l + j < endAt[j] and l + j2 < endAt[j2] and s[j+l] == s[j2+l]:
                    l += 1
                LCP[i-1] = l
            else:
                l = 0
        self.lcp = LCP

    def step3_rstr(self):
        prev_len = 0
        idx = 0
        results = {}
        len_lcp = len(self.lcp) - 1

        class Stack:
            pass
        stack = Stack()
        stack._top = 0
        stack.lst_max = []

        if len(self.res) == 0:
            return {}

        pos1 = self.res[0]
        for idx in range(len_lcp):
            current_len = self.lcp[idx]
            pos2 = self.res[idx+1]
            end_ = max(pos1, pos2) + current_len
            n = prev_len - current_len
            if n < 0:
                stack.lst_max.append([-n, idx, end_])
                stack._top += -n
            elif n > 0:
                self.removeMany(stack, results, n, idx)
            elif stack._top > 0 and end_ > stack.lst_max[-1][-1]:
                stack.lst_max[-1][-1] = end_

            prev_len = current_len
            pos1 = pos2

        if(stack._top > 0):
            self.removeMany(stack, results, stack._top, idx+1)

        return results

    def removeMany(self, stack, results, m, idxEnd):
        prevStart = -1
        while m > 0:
            n, idxStart, maxEnd = stack.lst_max.pop()
            if prevStart != idxStart:
                id_ = (maxEnd, idxEnd-idxStart+1)
                if id_ not in results or results[id_][0] < stack._top:
                    results[id_] = (stack._top, idxStart)
                prevStart = idxStart
            m -= n
            stack._top -= n
        if m < 0:
            stack.lst_max.append([-m, idxStart, maxEnd-n-m])
            stack._top -= m

    def go(self):
        self.step1_sort_suffix()
        self.step2_lcp()
        r = self.step3_rstr()
        return r


def GetMotifs(sequence):
    '''
    Given a sequence, return a list of maximal motifs with 
    length greater than 1
    Returns: [(motif length), [<start_indices>]]
    '''
    rstr = Rstr_max()
    rstr.add_str(sequence)
    r = rstr.go()
    result = []  # (word_index_start, word_index_end)
    for (_, nb), (l, start) in r.items():
        if l == 1: continue
        occurrences = [rstr.idxPos[rstr.res[o]]
                       for o in range(start, start+nb)]
        occurrences.sort()
        result.append((l, occurrences))
    return result


if (__name__ == '__main__'):
    str1 = [1, 1, 1, 2, 2, 1, 1]
    print(GetMotifs(str1))
