############################################################################################################################# 
#                                                                                                                           #
#   Description: A python program that decides, under certain conditions on h and g,                                        #
#                whether the word h(g^omega(a)) is additive k-power-free.                                                   #
#                                                                                                                           #
#   Conventions:  - Words are written as lists of integers from {0,1,2,...,A-1}, where A is the size of the alphabet        #
#                   E.g., the word 0101 is written [0,1,0,1] (Note that arbitrarily large alphabet sizes may be used)       #
#                 - While it is important that we write A-ary words with the alphabet {0,1,...,A-1},                        #
#                   one can assign arbitrary integer weights to these letters by means of a weight function, denoted WT     #
#                 - The integer -1 denotes the empty word                                                                   #
#                                                                                                                           # 
#   Main Variables:   k = size of additive power, entered as an integer > 1                                                 #
#                     t = additive k-template, e.g., t = [-1,-1,-1,[0,0]] is an additive 2-template whose instances         #
#                         are additive squares, and t = [0,-1,2,0,[-1,0],[1,2]] is an additive 3-template                   #
#                 f,g,h = morphisms, e.g., g = [[0,0,1],[0,1,2],[2,1,2]] represents the morphism                            #
#                         defined by g(0)=001, g(1)= 012, g(2) = 212                                                        #
#               A,Ag,Ah = alphabet sizes; Ag and Ah are sizes of output alphabets for g and h respectively                  #
#            WT,WTg,WTh = weightings for letters, e.g., WT = [0,1,3] indicates the weights 0, 1, and 3                      #
#                                for the letters 0, 1, and 2, respectively                                                  #
#                     a = starting letter, i.e., we consider the h(g^omega(a))                                              #
#                                                                                                                           #
#############################################################################################################################

import numpy as np # imports the NumPy library, designation np.

def IsFactor(u,w):  # Returns True if and only if u is a factor of w
    m=len(u)
    n=len(w)
    if m>n:
        return False
    for i in range(n-m+1):
        if w[i:i+m]==u:
            return True
    return False

def ContainsAll(U, w): # Returns True if and only if every word in the set U is a factor of w
    
    for u in U:

        if not IsFactor(u,w):
            return False
    return True

def Factors(s, n): # Returns a list containing all factors of length n of the word s
    factors = []
    for i in range(len(s) - n + 1):
        factor = s[i:i+n]
        if factor not in factors:
            factors.append(factor)
    return factors
    
def AdditiveKPower(w,k,WT): # Returns an additive k-power factor of w with letters weighted by WT if one exists, and returns None otherwise
    n = len(w)
    W = [0]
    T = 0

    for i in range(n):
        T += int(WT[int(w[i])])
        W.append(T)

    for r in range(1, n // k + 1):
        for i in range(n - k * r + 1):
            sums = [W[i + j * r] - W[i + (j - 1) * r] for j in range(1, k + 1)]
            if all(s == sums[0] for s in sums):
                return w[i : i + r*k]
    return None

def WeightedSum(w,WT): # Returns the sum of the letters of w with letters weighted by WT
    s = 0

    for i in range(len(w)): 
      s+= int(WT[w[i]])

    return s

def Sigma(w,WT): # Calculates the vector sigma of a word w with letters weighted by WT
    return np.array([(len(w)),(WeightedSum(w,WT))])

def Splits(f,A,WT): 

    # f is a morphism, and A is the size of the output alphabet of f
    # Returns a list of length A+1 whose ith element is the list of f-splits with the letter i at the center
    # (The last element in the output list, which can be called with Splits(f,A)[-1], 
    # is the list of f-splits with -1, i.e., the empty word, in the center.)
    # Each split is recorded in the form [sigma(p),a,sigma(s),A]

    splits = [[]]

    for i in range(A):
        splits.append([])

    for i in range(len(f)):
        w = f[i]

        for j in range(len(w)):
            center_element = int(w[j])
            
            if center_element < A:
                splits[center_element].append([Sigma(w[:j],WT), w[j], Sigma(w[j + 1:],WT), i])

    splits[-1].append([Sigma([],WT), -1, Sigma([],WT), -1])

    for i in range(len(f)):
        w = f[i]

        for j in range(len(w) + 1):
            splits[-1].append([Sigma(w[:j],WT), -1, Sigma(w[j:],WT), i])

    return splits

def Morphism(f, w): # Applies a morphism f to a word w
    W=[]
    for i in range(len(w)):
        W+=f[w[i]]

    return W

def IterMorphism(f, w, p): # Applies the morphism f p times to the word w
    for i in range(p):
        w = Morphism(f, w)
    return w

def OmegaFactorsPure(g, a, n): # Returns a list containing all factors of the word w = g^omega(a) of lenth n
     
    factors = []

    w=[a]

    while len(w) < n:
        w = Morphism(g, w)

    F = [w[0:n]]

    new_F = F

    while new_F:
        old_F = new_F
        new_F = []

        for u in old_F:
            new_factors = Factors(Morphism(g,u), n)
        
            for x in new_factors:
                if x not in factors:
                    factors.append(x)
                    new_F.append(x)

    return factors

def OmegaFactorsOuter(g,h,a,n): # Returns a list containing all factors of the word w = g^omega(a) of lenth n
    factors = []

    F_factors = OmegaFactorsPure(g, a, n)

    for u in F_factors:
        new_factors = Factors(Morphism(h,u), n)
        
        for w in new_factors:
            if w not in factors:
                factors.append(w)

    return factors

def IterateNeeded(f,a,n): 
    # Takes in a morphism f prolongable on a, and a positive integer n. 
    # Calculates the smallest integer j such that f^j(a) contains all factors of f^omega(a) of length n.

    factors=OmegaFactorsPure(f,a,n)

    w = [a]
    j = 0
    while True:
        j+=1
        w = Morphism(f,w)

        if ContainsAll(factors,w):
            return j

def MorphismMatrix(f,WT): # Returns the matrix M_f associated with the affine morphism f and the alphabet WT, or an error if f is not affine
    
    if all(len(f[i]) - len(f[i+1]) == len(f[i+1]) - len(f[i+2]) for i in range(len(f) - 2)):
        B = len(f[1]) - len(f[0])
        A = len(f[0]) - 2*f.index(f[0])
    else:
        raise ValueError("This Morphism doesn't have a linear relation between its lengths")

    sums =[]

    for k in range(len(f)):
        sums.append(WeightedSum(f[k],WT))
        
    if all(sums[j] - sums[j+1] == sums[j+1] - sums[j+2] for j in range(len(f) - 2)):
            D = sums[1]-sums[0]
            C = sums[0]-2*f.index(f[0])
    else:
        raise ValueError("This Morphism doesn't have a linear relation between the sums of digits")
    
    return [[A,B],[C,D]]

def Swap(matrix): # Swaps the diagonal and changes the sign of the anti-diagonal entries of a 2x2 matrix
    if len(matrix) != 2 or len(matrix[0]) != 2 or len(matrix[1]) != 2:
        raise ValueError("Input matrix must be a 2x2 matrix.")

    matrix[0][0], matrix[1][1] = matrix[1][1], matrix[0][0]

    matrix[0][1] = -matrix[0][1]
    matrix[1][0] = -matrix[1][0]

    return matrix

def Delta(t,k): # Returns Delta(t), where t is an additive k-template
    
    L =[]
    maximum = 0

    for m in range(k+1,len(t)):
        L.append(t[m])
    
    for j in range(len(L)):
        candidate = L[j][0]

        if candidate >= maximum:
            maximum = candidate

    return maximum  

def B(delta_t,k,f): # Returns B_f(t), where t is an additive k-template with Delta(t)=delta_t 
    W_f = 0

    for q in range(len(f)):
        W_f = max(len(f[q]),W_f)

    B_f = k + 2 + k*(W_f-2) + int((((k-1)*k)/2))*delta_t

    return B_f  

def Parents(t,k,f,A,WT): # Returns the set of all f-parents of the additive k-template t, where A is the size of the output alphabet of f, and WT is a weight function
   
    i = 0
    OldParents = [[[],[]]]
    splits = Splits(f,A,WT) 
    m = MorphismMatrix(f,WT)
    det_m = m[0][0]*m[1][1]-m[0][1]*m[1][0]
    m = Swap(m)

    while i <= k:
        NewParents = []
        
        for S in OldParents:
            
            for split in splits[t[i]]:
                
                newS = [S[0].copy(),S[1].copy()]
                
                newS[0].append(split)
                
                if i < 2:
                    NewParents.append(newS)                       
                else:
                    
                    b = np.subtract(np.add(newS[0][i-1][2],newS[0][i][0]),np.add(newS[0][i-2][2],newS[0][i-1][0]))
                    
                    x = np.subtract(t[k+1+i-2],b)
                    
                    D = np.matmul(m,x)
                    
                    if D[0] % det_m == 0 and D[1] % det_m == 0:
                        D[0] = D[0]/det_m
                        D[1] = D[1]/det_m
                        
                        newS[1].append(D)  
                        NewParents.append(newS)
                                       
        OldParents = NewParents
        i+=1
    
    parents_set = set()

    for S in OldParents:

        T =[]
        
        for s in S[0]:
            T.append(s[-1])
           
        T += S[1]  

        for i, element in enumerate(T):
            if isinstance(element, np.ndarray):
                T[i] = tuple(element)

        T_tuple = tuple(T)  # Convert the list to a tuple for hashing
        
        parents_set.add(T_tuple)  # Add the tuple to parents
    
    return parents_set

def AncestorsPure(T,k,f,A,WT): 
    
    # Computes the set of all f-ancestors of additive k-template T
    # Returns the set of ancestors and the number of generations i required to obtain all of them

    Ancestors = Parents(T,k,f,A,WT)

    prevAncestors = Ancestors

    print(f"The initial template has {len(Ancestors)} parents.")

    i=0

    while True:

        i+=1

        print(f"The total number of ancestors at the end of generation {i} is {len(Ancestors)}.")

        newAncestors = set()

        for T in prevAncestors:

            T = list(T)

            for j in range(len(T)):
                if isinstance(T[j], tuple):
                    T[j] = list(T[j])

            Tparent = Parents(T,k,f,A,WT)

            for t in Tparent:
                if t not in Ancestors:
                    newAncestors.add(t)

        if len(newAncestors) == 0:
            print(f"In generation {i+1}, we found no new ancestors.")
            print(f"The total number of ancestors is {len(Ancestors)}, and they were all found after {i} generation(s).")
            break   

        print(f"In generation {i+1} we found {len(newAncestors)} new ancestors.")

        for t in newAncestors:
            Ancestors.add(t)

        prevAncestors = newAncestors

    return (Ancestors,i)

def AncestorsOuter(T,k,g,Ag,WTg,h,Ah,WTh):

    # Computes the set of all g-ancestors of all h-parents of additive k-template T
    # Returns the set of ancestors and the number of generations i required to obtain all of them
    # Ag and Ah are the sizes of the output alphabets for g and h, respectively, and WTg and WTh are the corresponding weights

    Ancestors = Parents(T,k,h,Ah,WTh)

    prevAncestors = Ancestors

    print(f"The initial template has {len(Ancestors)} h-parents.  Now we calculate all g-ancestors of all of these h-parents.")

    i=0

    while True:

        print(f"The total number of ancestors at the end of generation {i} is {len(Ancestors)}.")

        newAncestors = set()

        for T in prevAncestors:

            T = list(T)

            for j in range(len(T)):
                if isinstance(T[j], tuple):
                    T[j] = list(T[j])

            Tparent = Parents(T,k,g,Ag,WTg)

            for t in Tparent:
                if t not in Ancestors:
                    newAncestors.add(t)

        if len(newAncestors) == 0:
            print(f"In generation {i+1}, we found no new ancestors.")
            print(f"The total number of ancestors is {len(Ancestors)}, and they were all found after {i} generation(s).")
            break   

        print(f"In generation {i+1} we found {len(newAncestors)} new ancestors.")

        for t in newAncestors:
            Ancestors.add(t)

        prevAncestors = newAncestors

        i+=1


    return (Ancestors,i)

def MainPure(k,g,A,WT,a): 

    # Returns True if g^omega(a) is additive k-power-free with letters weighted by WT, and False otherwise
    # A is the size of the output alphabet of g

    # Build the template t_0 for additive k-powers

    t=[]
    for i in range(k+1):
        t.append(-1)
    for i in range(k-1):
        t.append([0,0])

    print(f'We are deciding whether or not the word g^omega({a}) is additive {k}-power-free with letters weighted by {WT}, where \ng = {g}.')

    # Initial Check

    B_g = B(Delta(t,k),k,g)

    print(f'For the Initial Check, we have B_g={B_g}, so we examine all factors of g^omega({a}) of length {B_g-1}.')

    omega_factors = OmegaFactorsPure(g,a,B_g-1)

    for factor in omega_factors:
        power=AdditiveKPower(factor,k,WT)
        if power:
            print(f'The word h(g^omega({a})) contains an additive {k}-power, namely {power}.')
            return False

    # Calculating Ancestors

    print(f'Finding no additive {k}-power of length less than {B_g}, we proceed onto Calculating Ancestors.')

    Ancestors,i = AncestorsPure(t,k,g,A,WT)

    # Final Check

    max_delta = 0

    for template in Ancestors:
        delta = Delta(template,k)

        if delta > max_delta:
            max_delta = delta

    z = B(max_delta,k,g)

    print(f'We have M={max_delta} and B_M ={z}.')

    j=IterateNeeded(g,a,z-1)
    
    print(f'All factors of g^omega({a}) of length {z-1} are contained in g^{j}({a})')

    print(f'So for the Final Check, we examine the word g^{i+j}({a}).')

    prefix = IterMorphism(g,[a],i+j) 
    power = AdditiveKPower(prefix,k,WT)

    if  power == None:
        print(f'We find that the word g^{i+j}({a}) is additive {k}-power-free, hence so is g^omega({a}).')
        return True

    print(f'The word g^{i+j}({a}) has an additive {k}-power, namely {power}.')
    return False


def MainOuter(k,g,Ag,WTg,h,Ah,WTh,a):

    # Returns True if h(g^omega(a)) is additive k-power-free with letters weighted by WTh, and False otherwise
    # Ag and Ah are the sizes of the output alphabets of g and h, respectively, and letters of g^omega(a) are weighted by WTg

    # Build the template t_0 for additive k-powers

    t=[]
    for i in range(k+1):
        t.append(-1)
    for i in range(k-1):
        t.append([0,0])
    
    print(f'We are deciding whether or not the word h(g^omega({a})) is additive {k}-power-free with letters weighted by {WTh}, where \ng = {g} and \nh = {h}.')

    # Initial Check

    B_h = B(Delta(t,k),k,h)

    print(f'For the Initial Check, we have B_h={B_h}, so we examine all factors of h(g^omega({a})) of length {B_h-1}.')
    
    omega_factors = OmegaFactorsOuter(g,h,a,B_h-1)
    
    for factor in omega_factors:
        power=AdditiveKPower(factor,k,WTh)
        if power:
            print(f'The word h(g^omega({a})) contains an additive {k}-power, namely {power}.')
            return False
    

    # Calculating Ancestors

    print(f'Finding no additive {k}-power of length less than {B_h}, we proceed onto Calculating Ancestors.')

    Anc,i = AncestorsOuter(t,k,g,Ag,WTg,h,Ah,WTh)

    # Final Check

    max_delta = 0
    
    for template in Anc:

        delta = Delta(template,k)

        if delta > max_delta:
            max_delta = delta
    
    z = B(max_delta,k,g)

    print(f'We have M={max_delta} and B_M={z}.')

    j=IterateNeeded(g,a,z-1)

    print(f'All factors of g^omega({a}) of length {z-1} are contained in g^{j}({a}).')
    print(f'So for the Final Check, we examine the word h(g^{i+j}({a})).')

    prefix = Morphism(h,IterMorphism(g,[a],i+j))
    
    power=AdditiveKPower(prefix,k,WTh)

    if power==None:
        print(f'We find that the word h(g^{i+j}({a})) is additive {k}-power-free, hence so is h(g^omega({a})).')
        return True

    print(f'The word h(g^{i+j}({a})) has an additive {k}-power, namely {power}.' )
    return False 

print('------------------------------------------')

# Verifies Proposition 3.1 in J. Andrade and L. Mol, Avoiding additive powers in rich words, preprint, 2024.
MainPure(5,[[0,0,0,0,1],[0,1,1,0,1]],2,[0,1],0)

print('------------------------------------------')

# Verifies Proposition 4.1 in J. Andrade and L. Mol, Avoiding additive powers in rich words, preprint, 2024.
MainPure(4,[[1,0,0,0,1],[1,0,1,2,1,0,1],[1,0,1,2,2,2,1,0,1]],3,[0,1,2],1)

print('------------------------------------------')

# Verifies the additive 4-power-freeness of the word in J. Currie, L. Mol, N. Rampersad, and J. Shallit, Extending Dekking's construction of an infinite binary word avoiding abelian 4-powers, preprint, 2021.  Available at arxiv.org/abs/2111.07857v1
MainOuter(4,[[0,0,1],[0,1,2],[2,1,2]],3,[0,1,2],[[0,0,0,1,0,0,1,1,1,0,0,1,0,0,0,1,1,0,0,0,1,1],[0,0,0,1,0,0,1,1,1,0,0,1,1,1,0,1,1,0,0,0,1,1],[0,1,1,1,0,0,1,1,1,0,0,1,1,1,0,1,1,0,0,0,1,1]],2,[0,1],0)

print('------------------------------------------')

# Verifies Theorem 1 in F. M. Dekking, Strongly non-repetitive sequences and progression-free sets, J. Combin. Theory Ser. A 27 (1979), 181-185.
MainPure(4,[[0,0,0,1],[0,1,1]],2,[0,1],0)

print('------------------------------------------')

# Verifies Theorem 2.2 in J. D. Currie and A. Aberkane, A cyclic binary morphism avoiding Abelian fourth powers, Theoret. Comput. Sci. 410 (2009), 44-52.
MainPure(4,[[0,0,1,0,0,0,1,0,1,1,1,0,1,0,0,0,1,0,1,1,0,0,0,1,0],[1,1,0,1,1,1,0,1,0,0,0,1,0,1,1,1,0,1,0,0,1,1,1,0,1]],2,[0,1],0)

print('------------------------------------------')
