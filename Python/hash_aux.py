
#################################################################################
# Auxiliary functions for hashing in S_n (hash_sn.py)
#################################################################################  

import math
import factoradic
import random
import string

class ut:

    # Convert words into little endian format (swap bytes)
    def swap_bytes(words):
        little_words=''
        length=len(words)
        for i in range(0,length,8):
            little_words=little_words+words[length-(8+i):length-i]
        return little_words

    # Return random bits
    def random_bits(length):
        rand=''
        for i in range(length):
            rand=rand+str(random.randint(0,1))
        return rand

    # Return random letters
    def random_letters(length):
        rand=''
        for i in range(length):
            rand=rand + random.choice(string.ascii_letters)
        return rand

    # Returns number of elements that differ between permutations
    def hamming_distance(m1,m2):
        n=len(m1)
        m1,m2=([i for i in m1],[i for i in m2])
        count=0
        for i in range(n):
            if m1[i]==m2[i]:
                count=count+1
        return n-count

    # Inverse of permutation
    def per_inv(per):
        return [per.index(i) for i in range(len(per))]

    # Composition of permutations
    def comp(per_1,per_2):
        n=len(per_1)
        if n!=len(per_2):
            raise IndexError('Permutations are not same length!')
        return [per_1[per_2[i]] for i in range(n)]

##################################################################################
# Encodings and decodings: integer to/from factoradic to/from permutation
##################################################################################

class en:

    # Convert binary to integer
    def bits_to_int(bits):
        n=0
        exp=0
        for i in bits[::-1]:
            if i=='1':
                n=n+2**exp
            exp=exp+1
        return n
    
    # Convert bits to string
    def bits_to_string(bits):
        if len(bits) % 8 : raise Exception('Uncompleted byte!')
        string=''
        length=len(bits) // 8
        for i in range(length):
            string=string+chr(int(bits[i*8:(i+1)*8],2))
        return string

    # Convert string to bits
    def string_to_bits(string):
        bits=''
        for i in string:
            bits=bits+format(ord(i),'08b')
        return bits



    
    # Encode a number into a permutation using
    # the factoradic system and the Lehmer code
    def int_to_per(n,t):
        # If n>32!-1 exit
        if n>math.factorial(t)-1:
            raise OverflowError('Number exceeds '+str(t)+'!-1')
        # To factoradic system
        fact=factoradic.to_factoradic(n)
        fact.reverse()
        fact=[0 for i in range(t-len(fact))]+fact
        per,T=[],[i for i in range(t)]
        # To permutation
        for i in range(t):
            per.append(T[fact[i]])
            T.remove(T[fact[i]])
        return per

    # Integer to factoradic
    def int_to_fact(n):
        return factoradic.to_factoradic(n)

    # Decode a permutation into a number
    def per_to_int(per,t):
        # To factoradic system
        fact=[]
        n=len(per)
        for i in range(n):
            count=0
            for j in range(n-i):
                if per[i+j]<per[i]:
                    count=count+1
            fact.append(count)
        k=0
        # To integer
        for i in range(1,n):
            k=k+fact[-i-1]*math.factorial(i)
        return k

##################################################################################
# Constructions for (>>,<<). Shifts and rotations.
##################################################################################

class op:
    
    # a shifts b from the left
    def a_shift_b(a,b,t):

        if type(a) is int:
            a = en.int_to_per(a,t)
        if type(b) is int:
            b = en.int_to_per(b,t)
        if type(a) is str:
            a = en.bits_to_int(a)
            a = en.int_to_per(a,t)
        if type(b) is str:
            b = en.bits_to_int(b)
            b = en.int_to_per(b,t)
        
        for i in range(t):
            s=a[i]
            if s not in [0,t-1]:
                b1=b[:i]+b[i+1:]
                b1=[b1[(j-s)%(t-1)] for j in range(t-1)]
                b=b1[:i]+[b[i]]+b1[i:]
            else:
                if s==0:
                    if t%2==0:
                        b=[b[a.index((t-a[j])%t)] for j in range(t)]
                    else:
                        b=[b[a.index(t-1-a[j])] for j in range(t)]
                else:
                    b=[b[(j+1)%(t)] for j in range(t)]
        return b

    # a shifts b from the right
    def a_shift_b_inv(a,b,t):

        if type(a) is int:
            a = en.int_to_per(a,t)
        if type(b) is int:
            b = en.int_to_per(b,t)
        if type(a) is str:
            a = en.bits_to_int(a)
            a = en.int_to_per(a,t)
        if type(b) is str:
            b = en.bits_to_int(b)
            b = en.int_to_per(b,t)
        
        t=len(a)
        for i in range(t-1,-1,-1):
            s=a[i]
            if s not in [0,t-1]:
                b1=b[:i]+b[i+1:]
                b1=[b1[(j+s)%(t-1)] for j in range(t-1)]
                b=b1[:i]+[b[i]]+b1[i:]
            else:
                if s==0:
                    if t%2==0:
                        b=[b[a.index((t-a[j])%t)] for j in range(t)]
                    else:
                        b=[b[a.index(t-1-a[j])] for j in range(t)]
                else:
                    b=[b[(j-1)%(t)] for j in range(t)]
        return b


    # a rotates b from the left
    def a_rotate_b(a,b,t):

        if type(a) is int:
            a = en.int_to_per(a,t)
        if type(b) is int:
            b = en.int_to_per(b,t)
        if type(a) is str:
            a = en.bits_to_int(a)
            a = en.int_to_per(a,t)
        if type(b) is str:
            b = en.bits_to_int(b)
            b = en.int_to_per(b,t)
        
        t=len(a)
        index=[0 for i in range(t)]
        for i in range(t):
            index[a[i]]=i
        for i in range(t):
            b1=b.copy()
            if a[i]!=0:
                r=b[i]
                b=b[:i]+b[i+1:]
                if a[i]>t-1-i:
                    b=b[1:i-(t-a[i])+1]+[r]+b[i-(t-a[i])+1:]+[b[0]]
                else:
                    b=b[:i+a[i]]+[r]+b[i+a[i]:]
            else:
                if t%2==0:
                    b=[b[index[(t-a[j])%t]] for j in range(t)]
                else:
                    b=[b[index[t-1-a[j]]] for j in range(t)]
        return b

    # a rotates b from the right  
    def a_rotate_b_inv(a,b,t):

        if type(a) is int:
            a = en.int_to_per(a,t)
        if type(b) is int:
            b = en.int_to_per(b,t)
        if type(a) is str:
            a = en.bits_to_int(a)
            a = en.int_to_per(a,t)
        if type(b) is str:
            b = en.bits_to_int(b)
            b = en.int_to_per(b,t)
        
        t=len(a)
        a=[a[t-1-i] for i in range(t)]
        b=[b[t-1-i] for i in range(t)]
        for i in range(t):
            if a[i]!=0:
                k=(i+t-a[i])%t
                r=b[k]
                b=b[:k]+b[k+1:]
                if a[i]>i:
                    b=b[1:i+1]+[r]+b[i+1:]+[b[0]]
                else:
                    b=b[:a[i]+k]+[r]+b[a[i]+k:]
            else:
                if t%2==0:
                    b=[b[a.index((t-a[j])%t)] for j in range(t)]
                else:
                    b=[b[a.index(t-1-a[j])] for j in range(t)]
                    
        return [b[t-1-i] for i in range(t)]

