
#################################################################################
# Hashing based on operations (>>) and (<<) in the symmetric group
#################################################################################  

import math,factoradic
import numpy as np
import pdb,string,time
from random import randint
import matplotlib.pyplot as plt

def menu():
    
    q,L,path='',['0','1','2','3'],'C:\\Users\\Tyson\\Desktop\\hash_sn\\'
    
    # Ask for message option, manual or load file
    while q not in L:
        q=input("Input string (0), Input bits (1) Load string file (2) Load bit file (3)? ")
    if q==L[0]:
        m1=input('Input string message: ')
        return hash_sn(m1,0)
    elif q==L[1]:
        m2=input('Input bit message: ')
        return hash_sn(m2,1)
    elif q==L[2]:
        m3=input('Input file name: ')
        file=open(path+m3,'r',encoding="utf8").read()
        return hash_sn(file,0)
    else:
        m3=input('Input file name: ')
        file=open(path+m3,'r',encoding="utf8").read()
        return hash_sn(file,1)

def hash_sn(m,bit_conversion):

    # 2^128<35!, 2^256<58!, p=2^132-347 or 10101572818049573286907662744861646324507 from https://bigprimes.org/
    global fact,block_size,p
    
    (block_size,fact,p)=(128,35,2**132-347)

    #start=time.time()

    m,length=hash_sn_pad(m,bit_conversion,block_size)
    
    M,n,seq=hash_blocks(m,block_size)

    #print('Convert into blocks:',time.time()-start)

    # Initial value for avalanch effect -> message transformed into integer transformed into permutation in S_fact
    e=[i for i in range(fact)]
    a=a_rotate_b_inv(e,int_to_per((n+length)%math.factorial(fact),fact))

    # Start hashing
    for i in seq:
        i=int(i)
        if i:
            a=a_rotate_b_inv(a,int_to_per(M[i],fact))
        else:
            a=a_rotate_b(a,int_to_per(M[i],fact))
    
    #print('Processed blocks:',time.time()-start)
            
    return format(per_to_int(a,fact)%(2**block_size),'032x')

        
# Bit conversion and padding
def hash_sn_pad(m,bit_conversion,block_size):
    
    if bit_conversion==0:
        
        # Convert message to bits
        byte_array=bytearray(m,'utf8')
        m=''
        
        for byte in byte_array:
            m=m+format(byte,'08b')
       
    length=len(m)
    pad1=length%block_size
    half_block_size=int(block_size/2)
    pad2=m[-min(half_block_size,length):]
    
    if pad1==half_block_size: m=m+str(10**(block_size-1))
    elif pad1>half_block_size: m=m+str(10**((block_size-1)-(pad1-half_block_size)))
    else: m=m+str(10**((half_block_size-1)-pad1))

    pad2="".join(str(elem) for elem in [pad2[-i] for i in range(1,len(pad2)+1)])

    m=m+pad2.rjust(half_block_size,'0')
    
    return m,length

# block_size-bit block processing
def hash_blocks(m,block_size):
    
    blocknum=int(len(m)/block_size)
    M=[]
    n=j=1
    length_prime=len(format(p,'0b'))
    
    for i in range(1,blocknum+1):
        
        n1=block_size*(i-1)
        n2=block_size*i
        block=m[n1:n2]
        s=(i-1)%block_size
        block_to_int=bits_to_int(block[s:]+block[:s])
        M.append(block_to_int)
        n=(n*block_to_int)%p
        
        if i%block_size==block_size-1:
            
            upd=format(n,'0'+str(length_prime)+'b')
            n=bits_to_int(upd[j:]+upd[:j])
            j=(j+1)%length_prime

    # Sequence of operations
    seq=format(n,'0b')[:blocknum]
    length=len(seq)
    
    if length<blocknum:
        
        if (blocknum-length)%2==1: seq='01'*(int((blocknum-length)/2))+'0'+seq
        else: seq='01'*(int((blocknum-block_size)/2))+seq
        
    return M,n,seq

#################################################################################
# Ad hoc functions
#################################################################################  

# Convert binary to integer
def bits_to_int(bits):
    n=0
    exp=0
    for i in bits[::-1]:
        if i=='1':
            n=n+2**exp
        exp=exp+1
    return n

# Convert words into little endian format (swap bytes)
def swap_bytes(words):
    little_words=''
    length=len(words)
    for i in range(0,length,8):
        little_words=little_words+words[length-(8+i):length-i]
    return little_words

# Convert bits to string
def bits_to_string(bits):
    string=''
    length=int(len(bits)/8)
    for i in range(length):
        string=string+chr(int(bits[i*8:(i+1)*8],2))
    return string

# Return random bits
def random_bits(length):
    rand=''
    for i in range(length):
        rand=rand+str(randint(0,1))
    return rand

# Return random letters
def random_letters(length):
    rand=''
    for i in range(length):
        rand=rand+random.choice(string.ascii_letters)
    return rand


# Shifts with fixed element
##def a_rotate_b(a,b):
##    n=len(a)
##    for i in range(n):
##        if a[i]!=0:
##            if a[i]!=n-1:
##                b=[b[a.index(n-1-a[j])] for j in range(n)]
##            r=b[i]
##            b=b[:i]+b[i+1:]
##            if a[i]>n-1-i:
##                b=b[1:i-(n-a[i])+1]+[r]+b[i-(n-a[i])+1:]+[b[0]]
##            else:
##                #b=[b[a.index((n-a[j])%n)] for j in range(n)]
##                b=b[:i+a[i]]+[r]+b[i+a[i]:]  
##        else:
##            #for i in range()
##            a0,an=a.index(0),a.index(n-1)
##            b0,bn=b[a0],b[an]
##            b[a0],b[an]=bn,b0
##        if len(set(b))!=n:
##            print('error!')
##    return b


##def a_rotate_b(a,b):
##    n=len(a)
##    for i in range(n-1):
##        j=i+n
##        b=b+b+b
##        b1=b[j-a[i]:j+a[i]+1]
##        n1=a[i]*2+1
##        b1=[b1[n1-1-i] for i in range(n1)]
##        b=b[:j-a[i]]+b1+b[j+a[i]+1:]
##        print(b[n:2*n])


# Cycles
def a_rotate_b(a,b):
    n=len(a)
    for i in range(n):
        b1=b.copy()
        if a[i]!=0:
            r=b[i]
            b=b[:i]+b[i+1:]
            if a[i]>n-1-i:
                b=b[1:i-(n-a[i])+1]+[r]+b[i-(n-a[i])+1:]+[b[0]]
            else:
                b=b[:i+a[i]]+[r]+b[i+a[i]:]
        else:
            if n%2==0:
                b=[b[a.index((n-a[j])%n)] for j in range(n)]
            else:
                b=[b[a.index(n-1-a[j])] for j in range(n)]
        #print(b)
        #print(per_to_int(b1,n)-per_to_int(b,n))
    return b
  
def a_rotate_b_inv(a,b):
    n=len(a)
    a=[a[n-1-i] for i in range(n)]
    b=[b[n-1-i] for i in range(n)]
    for i in range(n):
        if a[i]!=0:
            k=(i+n-a[i])%n
            r=b[k]
            b=b[:k]+b[k+1:]
            if a[i]>i:
                b=b[1:i+1]+[r]+b[i+1:]+[b[0]]
            else:
                b=b[:a[i]+k]+[r]+b[a[i]+k:]
        else:
            if n%2==0:
                b=[b[a.index((n-a[j])%n)] for j in range(n)]
            else:
                b=[b[a.index(n-1-a[j])] for j in range(n)]
    return [b[n-1-i] for i in range(n)]
        
# Encode a number into a permutation using
# the factoradic system and the Lehmer code
def int_to_per(n,base):
    # If n>32!-1 exit
    if n>math.factorial(base)-1:
        return 'Error: number exceeds '+str(base)+'!-1'
    # To factoradic system
    fact=factoradic.to_factoradic(n)
    fact.reverse()
    fact=[0 for i in range(base-len(fact))]+fact
    #print(fact)
    per,T=[],[i for i in range(base)]
    # To permutation
    for i in range(base):
        per.append(T[fact[i]])
        T.remove(T[fact[i]])
    return per

# Decode a permutation into a number
def per_to_int(per,base):
    # To factoradic system
    fact=[]
    n=len(per)
    for i in range(n):
        count=0
        for j in range(n-i):
            if per[i+j]<per[i]:
                count=count+1
        fact.append(count)
    #print(fact)
    k=0
    # To integer
    for i in range(1,n):
        k=k+fact[-i-1]*math.factorial(i)
    return k

def hamming_distance(m1,m2):
    n=len(m1)
    m1,m2=([i for i in m1],[i for i in m2])
    count=0
    for i in range(n):
        if m1[i]==m2[i]:
            count=count+1
    return n-count

# Given b and c find a such that a>>b=c
def rot_inv(b,c,inv):
    n=len(b)
    if len(c)!=n:
        return 'error size'
    a=int_to_per(randint(1,math.factorial(n)-1),n)
    if inv==1:
        while a_rotate_b_inv(a,b)!=c:
            a=int_to_per(randint(1,math.factorial(n)-1),n)
    else:
        while a_rotate_b(a,b)!=c:
            a=int_to_per(randint(1,math.factorial(n)-1),n)
    return a

# Inverse of permutation
def per_inv(per):
    return [per.index(i) for i in range(len(per))]

# Composition of permutations
def comp(per_1,per_2):
    if len(per_1)!=len(per_2):
        return 'Error, pers are not same length'
    n=len(per_1)
    return [per_1[per_2[i]] for i in range(n)]


def string_to_int(s):
    bits=''
    for i in s:
        bits=bits+format(ord(i),'08b')
    return bits_to_int(bits)


#################################################################################
# Testing
#################################################################################  

# order of a is the number of distinct elements it generates
# when iteratively rotating b and is constant.
def order(a):
    n=len(a)
    a0=a.copy()
    a=a_rotate_b(a0,a)
    count=1
    while a!=a0:
        a=a_rotate_b(a0,a)
        count=count+1
    return count

def order1(a):
    n=len(a)
    e=[i for i in range(n)]
    a=a_rotate_b(a,e)
    a0=a.copy()
    count=1
    a=comp(a,a)
    while a!=a0:
        a=comp(a0,a)
        count=count+1
    return count

# Given c, there exists n! pairs a,b s.t a>>b=c
# Given c and b there exists small amount of a s.t a>>b=c
def test_collide(c):
    n=len(c)
    M=[]
    for j in range(math.factorial(n)):
        b=int_to_per(j,n)
        for i in range(math.factorial(n)):
            a=int_to_per(i,n)
            if a_rotate_b(a,b)==c:
                pdb.set_trace()
                M.append([a,b])
    print('number of collisions:',len(M)) # n=7, c=[1, 6, 2, 3, 0, 5, 4]->5040
    L=[]
    for i in range(len(M)):
        if i==0:
            L.append([list(M[i][1]),1])
        if list(M[i][1]) not in list(np.transpose(np.array(L,dtype='object'))[0]):
            L.append([list(M[i][1]),1])
        else:
            L[-1][1]=L[-1][1]+1
    print('number of rotated elements:',len(L)) # n=7->898
    N=[]
    for i in range(len(L)):
        N.append(L[i][1])
    maxi=max(N)
    print('max number of rotator elements (a)',maxi,' for ',N.count(maxi),' rotated elements (b)') # n=7-> 38 a for 1 b
    print(L[N.index(maxi)]) # n=7->[[4, 3, 0, 5, 2, 1, 6], 38]
    print('average number of rotators for single rotated:',np.mean(N))
    print('variance number of rotators for single rotated:',np.var(N))

# Return >> in the form of compositions in S_n
def a_rotate_b_comp(a,b,direction):
    n=len(a)
    if len(b)!=n:
        return 'error size'
    if direction=='r':
        for i in range(n):
            if a[i]<n-1:
                b1=b.copy()
                b2=b.copy()
                b1.remove(b1[i])
                b1=b1[n-1-a[i]:]+b1[:n-1-a[i]]
                b=b1[:i]+[b[i]]+b1[i:]
            else:
                b=b[1:]+b[:1]
            print('a'+str(i)+':',comp(b,per_inv(b2)))
        return b
    if direction=='l':
        for i in range(n):
            if a[i]<n-1:
                b1=b.copy()
                b2=b.copy()
                b1.remove(b1[i])
                b1=b1[a[i]:]+b1[:a[i]]
                b=b1[:i]+[b[i]]+b1[i:]
            else:
                b=b[n-1:]+b[:n-1]
            print('a'+str(i)+':',comp(b,per_inv(b2)))
        return b

# return representatives of mappings a:S_n-->S_n where a(b)=a>>b    
def classes(n):
    partitions=5
    threshold=int((math.factorial(n)*0.5)/partitions)

    d1,d2,d3,d4,d5,d6={},{},{},{},{},{}
    
    e=[j for j in range(n)]
    start=time.time()
    start1=start
    for i in range(math.factorial(n)):
        if i%1000000==0:
            length=len(d1)+len(d2)+len(d3)+len(d4)+len(d5)+len(d6)
            
            print(str(round((i/math.factorial(n))*100,5))+'%'.ljust(5,' '),length)
            start=time.time()
        a=int_to_per(i,n)
        c=per_to_int(a_rotate_b(a,e),n)
        if len(d1)<threshold:
            if d1.get(c):
                d1[c]=d1[c]+1
            if not d1.get(c):
                d1[c]=1
        elif len(d2)<threshold:
            if d1.get(c):
                d1[c]=d1[c]+1
            elif d2.get(c):
                d2[c]=d2[c]+1
            else:
                d2[c]=1
        elif len(d3)<threshold:
            if d1.get(c):
                d1[c]=d1[c]+1
            elif d2.get(c):
                d2[c]=d2[c]+1
            elif d3.get(c):
                d3[c]=d3[c]+1
            else:
                d3[c]=1
        elif len(d4)<threshold:
            if d1.get(c):
                d1[c]=d1[c]+1
            elif d2.get(c):
                d2[c]=d2[c]+1
            elif d3.get(c):
                d3[c]=d3[c]+1
            elif d4.get(c):
                d4[c]=d4[c]+1
            else:
                d4[c]=1
        elif len(d5)<threshold:
            if d1.get(c):
                d1[c]=d1[c]+1
            elif d2.get(c):
                d2[c]=d2[c]+1
            elif d3.get(c):
                d3[c]=d3[c]+1
            elif d4.get(c):
                d4[c]=d4[c]+1
            elif d5.get(c):
                d5[c]=d5[c]+1
            else:
                d5[c]=1
        else:
            if d6.get(c):
                d6[c]=d6[c]+1
            if not d6.get(c):
                d6[c]=1

    print(time.time()-start1)

    S1=list(d1.values())
    S2=list(d2.values())
    S3=list(d3.values())
    S4=list(d4.values())
    S5=list(d5.values())
    S6=list(d6.values())

    K=list(set(S1))
    K=list(set(K+S2))
    K=list(set(K+S3))
    K=list(set(K+S4))
    K=list(set(K+S5))
    K=list(set(K+S6))

    maxi=max(K)
    print('')
    print('Max:',maxi)
    print('Min:',min(K))

    K.sort()
    print('')
    print('k'.ljust(5,' '),'Classes containing k elements'.ljust(35,' '))
    print(''.ljust(45,'-'))

    m=0
    d={}
    for i in range(len(K)):

        d[K[i]]=S1.count(K[i])+S2.count(K[i])+S3.count(K[i])+S4.count(K[i])+S5.count(K[i])+S6.count(K[i])
        
        print(str(K[i]).ljust(5,' '), d[K[i]])
        m=m+d[K[i]]
        
    print('Number of mapping classes sigma_a:',m)

#####
#####    n=12
#####    d={1:64313031,2:57725732,3:38451932,4:21182580,5:10205298,6:4448803,7:1791110,8:676342,9:242495,10:83621,11:27720,12:8811,13:2815,
#####       14:827,15:256,16:74,17:28,18:5,19:5}
#####    K=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
#####
    
    plt.grid()
    plt.plot([K[i] for i in range(len(K))],list(d.values()))
    plt.title('n='+str(n))
    plt.xlabel('k')
    plt.ylabel('number of classes with k elements')
    plt.show()

def copies(n,overlay):
    plt.grid()
    plt.title('Number of copies of mappings n=5,6,7,8')#+str(n))
    plt.xlabel('a in Sn')
    plt.ylabel('copies of a')
    if overlay:
        threshold=5
        for k in range(threshold,n+1):
            e=[i for i in range(k)]
            L=[]
            for i in range(math.factorial(threshold)):
                r=a_rotate_b(int_to_per(i,k),e)
                count=0
                for j in range(math.factorial(k)):
                    if a_rotate_b(int_to_per(j,k),e)==r:
                        count+=1
                L.append(count)
            # Number of classes
            plt.scatter([i for i in range(math.factorial(threshold))],L[:math.factorial(threshold)])
        plt.show()
    else:
        threshold=math.factorial(n)
        e=[i for i in range(n)]
        L=[]
        for i in range(threshold):
            r=a_rotate_b(int_to_per(i,n),e)
            count=0
            for j in range(math.factorial(n)):
                if a_rotate_b(int_to_per(j,n),e)==r:
                    count+=1
            L.append(count)
        # Number of classes
        plt.plot([i for i in range(threshold)],L)
        plt.show()

        
    
def plots(plot):
    if plot==1:
        # Plot Number of mappings as n grows vs exp(n)
        L=[[9,153809],[10,1511923],[11,16843800]]
        L_array_tr=np.transpose(np.array(L))
        plt.grid()
        plt.title('Number of mappings as n grows vs n! and exp(n)')
        plt.xlabel('n')
        plt.ylabel('number of mappings in tens of millions')
        # Number of classes
        plt.plot(list(L_array_tr[0]),list(L_array_tr[1]),color='b')
        # Number of elements
        plt.plot(list(L_array_tr[0]),[math.factorial(L_array_tr[0][i]) for i in range(len(L_array_tr[0]))],color='g')

        partition=[9+i/100 for i in range(200)]
        plt.plot(partition,[math.exp(partition[i]) for i in range(len(partition))],color='r')
        plt.show()
        
    if plot==2:
        # number of mappings as n grows
        L=[[5,47],[6,301],[7,2149],[8,16884],[9,153809],[10,1511923],[11,16843800],[12,199161485]]

        x=list(np.transpose(L[2:])[0])
        y=[list(np.transpose(L[2:])[1])[i]/math.factorial(x[i]) for i in range(len(L[2:]))]
        a = np.vstack([x, np.ones(len(x))]).T
        R=list(np.dot(np.linalg.inv(np.dot(a.T, a)), np.dot(a.T, y)))
        
        L_array_tr=np.transpose(np.array(L))
        plt.grid()
        plt.title('ratio of classes')
        plt.xlabel('n')
        plt.ylabel('number of classes/n!')
        # ratio of classes vs n! ((number of classes/n!))
        plt.plot(list(L_array_tr[0]),[(L_array_tr[1][i]/math.factorial(L_array_tr[0][i])) for i in range(len(L))],color='b')

        # Regression least squares
        plt.plot([5,12],[R[0]*5+R[1],R[0]*12+R[1]],'r')

        # Number of elements
        #plt.plot(list(L_array_tr[0]),[math.factorial(L_array_tr[0][i]) for i in range(len(L_array_tr[0]))],color='g')
        plt.show()

def randomness(n,R):
    for i in range(R[0],R[1]):
        L=[]
        a=int_to_per(i,n)
        for j in range(math.factorial(n)):
            L.append(per_to_int(comp(a,int_to_per(j,n)),n))
            #L.append(per_to_int(a_rotate_b(a,int_to_per(j,n)),n))
        #plt.grid()
        plt.title('a='+str(a))
        plt.xlabel('b')
        plt.ylabel('a>>b')
        plt.scatter([k for k in range(math.factorial(n))],L)
        plt.savefig('randomness_'+str(i)+'_'+str(n)+'.png')
        #plt.show()
        #plt.close()
    
def one_block_domain(s,t):
    # s -> blocksize multiple of 8 (1 char -> 1 byte -> 8 bits)
    # t -> smallest factorial larger than 2**s
    # ex: one_block_domain(16,9)
    L=[]
    fact=math.factorial(t)
    blocksize=2**s
    length=1
    l=int_to_per(length*(2**(int(s/2))),t)
    for i in range(2**(int(s/2))):
        if i%10000==0:
            print(i)
        r=int_to_per(i*2**(int(s/2)),t)
        block=int_to_per(i,t)
        h0=a_rotate_b_inv(l,r)
        if h0[0]%2==0:
            result=per_to_int(a_rotate_b_inv(a_rotate_b(h0,block),block),t)%blocksize
        else:
            result=per_to_int(a_rotate_b_inv(a_rotate_b_inv(h0,block),block),t)%blocksize
        #if result not in L:
        L.append(result)
    length=2
    l=int_to_per(length*(2**(int(s/2))),t)
    for i in range(2**(int(s/2)),blocksize):
        if i%10000==0:
            print(i)
        r=int_to_per(bits_to_int(format(i,'016b')[-8:])*2**(int(s/2)),t)
        block=int_to_per(i,t)
        if h0[0]%2==0:
            result=per_to_int(a_rotate_b_inv(a_rotate_b(h0,block),block),t)%blocksize
        else:
            result=per_to_int(a_rotate_b_inv(a_rotate_b_inv(h0,block),block),t)%blocksize
        #if result not in L:
        L.append(result)
    print('Plotting...')
    plt.grid()
    plt.scatter([i for i in range(blocksize)],L)
    plt.show()

#################################################################################
# Hash tables using (>>,<<)
#################################################################################  


# Input normal distribution array
def find_opt_hash():

    # 1. Generate 100000 random keys
    L=[randint(0,2**128-1) for i in range(100000)]

    # 2. Normally distributed keys
    L=list(set(np.round(np.random.normal(500000,18800,1000000))))
        
    fact=35 # we assume keys are at most 128 bits
    H=[]
    M=[]
    for i in range(50):
        if i%10==0:
            print(i)
        n=len(L)
        #mod=int(n*1.3)
        mod=200003
        d={}
        #N=[]
        # pick random hash function
        a=int_to_per(randint(0,math.factorial(fact)-1),fact)
        #b=int_to_per(randint(0,math.factorial(fact)-1),fact)
        for j in range(n):
            #h=a_rotate_b_inv(b,int_to_per(round(L[j]),fact))
            #h=per_to_int(a_rotate_b_inv(a,h),fact)%mod
            h=per_to_int(a_rotate_b(a,int_to_per(round(L[j]),fact)),fact)%mod
            m=randint(0,mod-1)
            if d.get(h):
                d[h]+=1
            else:
                d[h]=1
            #N.append(h)
        #pdb.set_trace()
        #plt.title('Hashed normal keys, n=m*1.3, t='+str(fact))
        #plt.hist(set(N), bins = 1000)
        #print((len(d)/n)*100,max(d.values()))
        #plt.show()
        #plt.close()
        H.append([round((len(d)/n)*100,2),max(d.values())])
        H_tr=np.array(H).transpose()
    H_maxi=H_tr[0].max()
    H_mini=H_tr[0].min()
    print('mean:',H_tr[0].mean(),H_tr[1].mean())
    print('max:',H_maxi,H[list(H_tr[0]).index(H_maxi)][1])
    print('min:',H_mini,H[list(H_tr[0]).index(H_mini)][1])
    # SHA1
##    e={}
##    for j in range(n):
##        m=(hashlib.md5(str(j).encode())).hexdigest()
##        if e.get(m):
##            e[m]+=1
##        else:
##            e[m]=1
##    M.append([(len(e)/n)*100,max(e.values())])
    #return H

