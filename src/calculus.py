from functools import reduce


class Calculus(object):
    """
    docstring
    """
    def __init__(self,calc_name='rcc8',comp_method='full') -> None:
        self.calc_name = calc_name
        self.B_dict = {}
        self.B_dict_reverse = {}

        with open(self.calc_name+"/"+self.calc_name+".relations") as f:
            for i, j in enumerate(f):
                self.B_dict[j.strip()] = int(2**i)
                self.B_dict_reverse[int(2**i)] =  j.strip()
        
        self.size = i+1

        self.B_dict['DALL'] = 2**(self.size) - 1 

        self.B = sorted(self.B_dict.values())[:-1]

        # initialize and fill list for set based on base relations
        self.bsplit = [(len(self.bitdecoding(i+1)),self.bitdecoding(i+1)) for i in range((2**self.size)-1)]
        self.bsplit_dict = {i+1:self.bitdecoding(i+1) for i in range(2**self.size-1)}
        self._init_composition_table(comp_method)

        inverseDict = {}
        with open(self.calc_name+"/"+self.calc_name+".conv") as f:
            for i in f:
                rel, irel = i.split('::')
                inverseDict[self.B_dict[rel.strip()]] = self.B_dict[irel.strip()]
            
        self._inv = [(reduce(lambda x, y: x | y, [inverseDict[i] for i in self.bsplit[j][1]])) for j in range((2**self.size)-1)]
        self.inverse = self._get_inverse
        with open(self.calc_name+"/"+self.calc_name+".identity") as f:
            self.Id = self.B_dict[f.readline().strip()]

    # translate a base relation from integer to its string representation
    def translate(self,BR):
        return self.B_dict_reverse[BR]

    # translate a base relation from its string representation to integer
    def translateR(self,BR):
        return self.B_dict[BR]

    # split a relation into its base relation representation
    def bitdecoding(self,b):
        l = []

        if b in self.B: return [b]

        if b == self.B_dict['DALL']: return self.B[:]

        l = [i for i in self.B if b&i != 0 and i <=b]

        return l
    
    def _init_composition_table(self,method='full'):
        if method == 'full':
            DALL = self.B_dict['DALL']
            def start(i, fcomp):
                relations, composition = i.split('::')
                relA, relB = relations.strip().split(':')
                composition = composition.replace(')','').replace('(','').strip().split()
                fcomp[self.B_dict[relA.strip()]-1][self.B_dict[relB.strip()]-1] = reduce(lambda x, y: x | y, [self.B_dict[i] for i in composition])

            def complete(fcomp):
                for i in range((2**self.size)-1):
                    for j in range((2**self.size)-1):
                        comp = 0 
                        for m in self.bsplit[i][1]:
                            for n in self.bsplit[j][1]:
                                comp |= fcomp[m-1][n-1]
                                if comp == DALL: 
                                    break # stop calculating if the global relation
                            else:
                                continue  # executed if the loop ended normally (no break)
                            break  # executed if 'continue' was skipped (break)
                        fcomp[i][j] = comp

            # initialize matrix to hold compositions between all 256 possible relations
            with open(self.calc_name+"/"+self.calc_name+".comp") as f:
                from array import array
                if self.size <= 8:
                    self._fcomp = [array('B',[0 for i in range((2**self.size)-1)]) for j in range((2**self.size)-1)]
                elif self.size <= 16:
                    self._fcomp = [array('H',[0 for i in range((2**self.size)-1)]) for j in range((2**self.size)-1)]
                else:
                    self._fcomp = [array('I',[0 for i in range((2**self.size)-1)]) for j in range((2**self.size)-1)]
                
                for i in f:
                    start(i,self._fcomp)
                complete(self._fcomp)
                self.comp = self._get_composition_full
        elif method == 'basic':
            # self._b_to_idx = {r:i for i,r in enumerate(self.B)}
            self._bcomp = {r1:{r2:0 for r2 in self.B} for r1 in self.B}
            with open(self.calc_name+"/"+self.calc_name+".comp") as f:
                for i in f:
                    relations, composition = i.split('::')
                    relA, relB = relations.strip().split(':')
                    composition = composition.replace(')','').replace('(','').strip().split()
                    self._bcomp[self.B_dict[relA.strip()]][self.B_dict[relB.strip()]] = reduce(lambda x, y: x | y, [self.B_dict[i] for i in composition])
            self.comp = self._get_composition_basic
        else:
            raise Exception("Composistion method should either be 'basic' or 'full'!")
    
    def _get_composition_full(self,r1,r2):
        return self._fcomp[r1-1][r2-1]
    
    def _get_composition_basic(self,r1,r2):
        r1_split = self.bsplit_dict[r1]
        r2_split = self.bsplit_dict[r2]
        composition = 0
        for b1 in r1_split:
            for b2 in r2_split:
                composition |= self._bcomp[b1][b2]

        return composition
    
    def _get_inverse(self,r):
        return self._inv[r-1]




            



    