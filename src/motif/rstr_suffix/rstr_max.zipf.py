#!/usr/bin/env python
# -*- coding: utf-8 -*-

from tools_karkkainen_sanders import *
from rstr_max import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import math

str1 = "Otto von Bolschwing est mort sous le soleil de Sacramento (Californie), le 9 mars 1982. Il avait 72 ans et plusieurs vies derrière lui. Dans la première, ce rejeton de la noblesse allemande, membre des services de renseignement de la SS - l'ossature nazie du IIIe Reich - a exercé ses talents d'agent secret en Palestine. Devenu aide de camp du responsable des Affaires juives, Adolf Eichmann, il s'est attelé avec zèle à l'élaboration du programme d'extermination des juifs d'Europe. Après la guerre, le capitaine von Bolschwing a rejoint les rangs de la CIA, sous le nom de code de Grossbahn. En récompense de ses loyaux services, ses nouveaux maîtres lui ont offert une deuxième vie aux Etats-Unis, où l'ancien SS s'est taillé une belle carrière dans l'industrie. Jusqu'à devenir, en 1969, président de l'entreprise californienne de high-tech TCI."

str1 = open('002.art').read()
str1_unicode = unicode(str1,'utf-8','replace')
#str2_unicode = str1_unicode[::-1]

rstr = Rstr_max()
rstr.add_str(str1_unicode)
r = rstr.go()

#for ((idStr, end), nb), (l, start_plage, end_plage) in r.iteritems():
#  ss = rstr.array_str[idStr][end-l:end]
#  print ss
#  for o in range(start_plage, end_plage+1) :
#    offset, id_str = rstr.array_suffix[o]
#    print '   %d %d'%(offset, id_str)
#    sss = rstr.array_str[id_str][offset:offset+l]

l = r.keys()

def cmpval(x,y):
  return y[1] - x[1]

l.sort(cmpval)
list_x = [math.log(elt[1]) for elt in l]
list_y = [math.log(i) for i in range(1, len(list_x)+1)]

plt.plot(list_y, list_x, 'r-')
filename = 'zipf_rstr.png'
plt.xlabel('log(rank)')
plt.ylabel('log(freq.)')
plt.savefig(filename)
