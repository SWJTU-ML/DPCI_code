def MCS(ConMatrix, neighbors):
    
   a = [None for i in range(len(ConMatrix))]
   a_ = [None for i in range(len(ConMatrix))]     

   Set = [None for i in range(len(ConMatrix))]
   Size = [None for i in range(len(ConMatrix))] 

   for i in range(len(ConMatrix)):
      Set[i] = set([])     
      Size[i] = 0
      Set[0].add(i)

   j = 0
   for i in range(len(ConMatrix)):
      v = Set[j].pop()
      a[v] = i
      a_[i] = v
      if j==(len(ConMatrix)-1):
         break
      Size[v] = -1
      for u in neighbors[v]:
         if Size[u] >= 0:
            Set[Size[u]].remove(u)
            Size[u] = Size[u] + 1
            Set[Size[u]].add(u)
      j = j + 1
      while j >= 0 and not Set[j]:
         j = j - 1

   return a, a_

def FIC(ConMatrix, neighbors, pair ):
   a = pair[0]
   a_ = pair[1]
   f = [None for i in range(len(ConMatrix))]
   index = [None for i in range(len(ConMatrix))]
   fill = list()

   for i in range(len(ConMatrix)-1,-1,-1):
      w = a_[i]
      f[w] = w
      index[w] = i
      for v in neighbors[w]:
         if a[v] > i:
            x = v
            while index[x] > i:
               index[x] = i
               fill.append((x,w) if x<w else (w,x))
               x = f[x]
            if f[x] == x: f[x] = w

   return fill
