#!/usr/bin/env python
# -*- coding: utf-8 -*-
import unittest
from rstr_max import Rstr_max

class Test_rstrmax:
  def setUp(self):
    self.list_s = self.getString()
    self.rstr = Rstr_max()
    for s in self.list_s :
      self.rstr.add_str(s)

  def test_rstr_max(self) :
    r = self.rstr.go()
    for (offset_end, nb), (l, start_plage) in r.iteritems():
      ss = self.rstr.global_suffix[offset_end-l:offset_end]
#      ss = self.rstr.array_str[idStr][end-l:end]
      offset_end -= 1
      id_chaine = self.rstr.idxString[offset_end]
      s = self.rstr.global_suffix
      idx = 0
      for i in xrange(nb):
        idx = s.index(ss, idx) + 1
#      except ValueError, e:
#        print "+++", ss, end, i, nb
#      try:
#      self.assertRaises(ValueError, s.index, ss, idx)
#        print "***", ss, end, i, nb
#      except ValueError, e:
#        pass

  def test_maximality(self) :
    r = self.rstr.go()
    for (offset_end, nb), (l, start_plage) in r.iteritems():
      ss = self.rstr.global_suffix[offset_end-l:offset_end]
      offset_end -= 1
      id_chaine = self.rstr.idxString[offset_end]
      s = self.rstr.array_str[id_chaine]

      set_left, set_right = set(), set()

      for o in range(start_plage, start_plage + nb) :
        offset_global = self.rstr.res[o]
        su = (self.rstr.idxPos[offset_global],self.rstr.idxString[offset_global])
        ls = len(self.rstr.array_str[su[1]])

        char_left = "START_STR%i"%(su[1]) if(su[0] == 0) else self.rstr.array_str[su[1]][su[0]-1]
        set_left.add(char_left)

        char_right = "END_STR%i"%(su[1]) if(su[0]+l == ls) else self.rstr.array_str[su[1]][su[0]+l]
        set_right.add(char_right)

      self.assertNotEqual(len(set_left), 1)
      self.assertNotEqual(len(set_right), 1)

  def utest_left_maximality(self) :
    r = self.rstr.go()
#    for (idStr, end, nb), (l, start_plage) in r.iteritems():
    for (offset_end, nb), (l, start_plage) in r.iteritems():
      ss = self.rstr.global_suffix[offset_end-l:offset_end]
#      ss = self.rstr.array_str[idStr][end-l:end]
      offset_end -= 1
      id_chaine = self.rstr.idxString[offset_end]
      s = self.rstr.array_str[id_chaine]
#      s = self.rstr.array_str[idStr]
      set_left_char = set()
      for o in range(start_plage, start_plage + nb) :
        offset_global = self.rstr.res[o]
        su = (self.rstr.idxPos[offset_global],self.rstr.idxString[offset_global])
#        su = self.rstr.array_suffix[o]
        if(su[0] == 0) :
          char_left = "START_STR"
        else :
          char_left = self.rstr.array_str[su[1]][su[0]-1]
        set_left_char.add((char_left,su[1]))
      if(len(set_left_char) == 1) :
        print
        print '*'*10
        print set_left_char
        print ss.encode('utf-8')
        print '*'*10
        print
      self.assertNotEqual(len(set_left_char), 1)


  def utest_right_maximality(self) :
    r = self.rstr.go()
    for (offset_end, nb), (l, start_plage) in r.iteritems():
      ss = self.rstr.global_suffix[offset_end-l:offset_end]
      offset_end -= 1
      id_chaine = self.rstr.idxString[offset_end]
      s = self.rstr.array_str[id_chaine]
      set_right_char = set()
      for o in range(start_plage, start_plage + nb) :
        offset_global = self.rstr.res[o]
        su = (self.rstr.idxPos[offset_global],self.rstr.idxString[offset_global])
        ls = len(self.rstr.array_str[su[1]])
        if(su[0]+l == ls) :
          char_right = "END_STR"
        else :
          char_right = self.rstr.array_str[su[1]][su[0]+l]
        set_right_char.add((char_right,su[1]))

      if(len(set_right_char) == 1) :
        print
        print '*'*10
        print set_right_char
        print ss.encode('utf-8')
        print '*'*10
        print
      self.assertNotEqual(len(set_right_char), 1)

class Test_rstrmax_test1(Test_rstrmax, unittest.TestCase):
  def getString(self):
    str1 = " u u"
    return unicode(str1,'utf-8','replace')

#class Test_rstrmax_otto(Test_rstrmax, unittest.TestCase):
#  def getString(self):
#    str1 = open('otto.txt','r').read()
#    return unicode(str1,'utf-8','replace')
   
class iTest_rstrmax_python(Test_rstrmax, unittest.TestCase):
  def getString(self):
    str1 = open('Python.htm','r').read()
    return unicode(str1,'utf-8','replace')[:1000]

class Test_rstrmax_a(Test_rstrmax, unittest.TestCase):
  def getString(self):
    return ['a'*50]

#class Test_rstrmax_tititoto(Test_rstrmax, unittest.TestCase):
#  def getString(self):
#    return ['tititoto']

#class Test_rstrmax_toto(Test_rstrmax, unittest.TestCase):
#  def getString(self):
#    return ['toto']

class Test_rstrmax_one(Test_rstrmax, unittest.TestCase):
  def getString(self):
    return ['']


class Test_rstrmax_one_two(Test_rstrmax, unittest.TestCase):
  def getString(self):
    return ['a','']

#class iTest_rstrmax_art(Test_rstrmax, unittest.TestCase):
#  def getString(self):
#    str1 = open('002.art','r').read()
#    return unicode(str1,'utf-8','replace')

if (__name__ == '__main__') :
  unittest.main()

